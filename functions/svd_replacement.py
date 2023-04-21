import torch

class H_functions:
    """
    A class replacing the SVD of a matrix H, perhaps efficiently.
    All input vectors are of shape (Batch, ...).
    All output vectors are of shape (Batch, DataDimension).
    """

    def V(self, vec):
        """
        Multiplies the input vector by V
        """
        raise NotImplementedError()

    def Vt(self, vec):
        """
        Multiplies the input vector by V transposed
        """
        raise NotImplementedError()

    def U(self, vec):
        """
        Multiplies the input vector by U
        """
        raise NotImplementedError()

    def Ut(self, vec):
        """
        Multiplies the input vector by U transposed
        """
        raise NotImplementedError()

    def singulars(self):
        """
        Returns a vector containing the singular values. The shape of the vector should be the same as the smaller dimension (like U)
        """
        raise NotImplementedError()

    def add_zeros(self, vec):
        """
        Adds trailing zeros to turn a vector from the small dimension (U) to the big dimension (V)
        """
        raise NotImplementedError()
    
    def H(self, vec):
        """
        Multiplies the input vector by H
        """
        temp = self.Vt(vec)
        singulars = self.singulars()
        return self.U(singulars * temp[:, :singulars.shape[0]])
    
    def Ht(self, vec):
        """
        Multiplies the input vector by H transposed
        """
        temp = self.Ut(vec)
        singulars = self.singulars()
        return self.V(self.add_zeros(singulars * temp[:, :singulars.shape[0]]))
    
    def H_pinv(self, vec):
        """
        Multiplies the input vector by the pseudo inverse of H
        """
        temp = self.Ut(vec)
        singulars = self.singulars()
        temp[:, :singulars.shape[0]] = temp[:, :singulars.shape[0]] / singulars
        return self.V(self.add_zeros(temp))

#a memory inefficient implementation for any general degradation H

class ultrasound1(H_functions):
    def __init__(self, channels, lbd, V, device):
        self.channels = channels
        self._singulars = lbd.repeat(self.channels).reshape(1,self.channels,-1).permute(0,2,1).reshape(-1).to(device) #torch.ones(3 * 256 ** 2, device=device)
        self._V = V
        self._Vt = self._V.transpose(0, 1)

        # ZERO = 2.4#1e-3
        # self._singulars[self._singulars < ZERO] = 0
        # print(len([x.item() for x in self._singulars if x == 0]))

    def V(self, vec):
        #return self.mat_by_vec(self._V, vec.clone()).to(vec.device)
        return torch.matmul(self._V, vec.clone().reshape(vec.shape[0], -1, self.channels).to('cpu')).permute(0,2,1).reshape(vec.shape[0],-1).to(vec.device)

    def Vt(self, vec):
        # return self.mat_by_vec(self._Vt, vec.clone()).to(vec.device)
        return torch.matmul(self._Vt, vec.clone().reshape(vec.shape[0], self.channels, -1).permute(0,2,1).to('cpu')).reshape(vec.shape[0],-1).to(vec.device)

    def U(self, vec):
        # return self.mat_by_vec(self._U, vec.clone()).to(vec.device)
        return vec.clone().reshape(vec.shape[0], -1, self.channels).permute(0, 2, 1).reshape(vec.shape[0], -1)

    def Ut(self, vec):
        # return self.mat_by_vec(self._Ut, vec.clone()).to(vec.device)
        return vec.clone().reshape(vec.shape[0], self.channels, -1).permute(0, 2, 1).reshape(vec.shape[0], -1)

    def singulars(self):
        return self._singulars

    def add_zeros(self, vec):
        return vec.clone().reshape(vec.shape[0], -1)
        # out = torch.zeros(vec.shape[0], self._V.shape[0]*3, device=vec.device)
        # out[:, :self._U.shape[0]*3] = vec.clone().reshape(vec.shape[0], -1)
        # return out


class ultrasound0(H_functions):
    def __init__(self, channels, U, lbd, V, device):
        #self._U, self._singulars, self._V = torch.svd(H, some=False)
        self.channels = channels
        self._singulars = lbd.repeat(self.channels).reshape(1,self.channels,-1).permute(0,2,1).reshape(-1).to(device) #torch.ones(3 * 256 ** 2, device=device)
        self._U = U
        self._V = V
        self._Vt = self._V.transpose(0, 1)
        self._Ut = self._U.transpose(0, 1)

        # ZERO = 40#1e-3
        # self._singulars[self._singulars < ZERO] = 0
        # print(len([x.item() for x in self._singulars if x == 0]))

    def V(self, vec):
        #return self.mat_by_vec(self._V, vec.clone()).to(vec.device)
        return torch.matmul(self._V, vec.clone().reshape(vec.shape[0], -1, self.channels).to('cpu')).permute(0,2,1).reshape(vec.shape[0],-1).to(vec.device)

    def Vt(self, vec):
        # return self.mat_by_vec(self._Vt, vec.clone()).to(vec.device)
        return torch.matmul(self._Vt, vec.clone().reshape(vec.shape[0], self.channels, -1).permute(0,2,1).to('cpu')).reshape(vec.shape[0],-1).to(vec.device)

    def U(self, vec):
        # return self.mat_by_vec(self._U, vec.clone()).to(vec.device)
        return torch.matmul(self._U, vec.clone().reshape(vec.shape[0], -1, self.channels).to('cpu')).permute(0,2,1).reshape(vec.shape[0],-1).to(vec.device)

    def Ut(self, vec):
        # return self.mat_by_vec(self._Ut, vec.clone()).to(vec.device)
        return torch.matmul(self._Ut, vec.clone().reshape(vec.shape[0], self.channels, -1).permute(0, 2, 1).to('cpu')).reshape(vec.shape[0], -1).to(vec.device)

    def singulars(self):
        return self._singulars

    def add_zeros(self, vec):
        return vec.clone().reshape(vec.shape[0], -1)
        # out = torch.zeros(vec.shape[0], self._V.shape[0]*3, device=vec.device)
        # out[:, :self._U.shape[0]*3] = vec.clone().reshape(vec.shape[0], -1)
        # return out

#Inpainting
class Inpainting(H_functions):
    def __init__(self, channels, img_dim, missing_indices, device):
        self.channels = channels
        self.img_dim = img_dim
        self._singulars = torch.ones(channels * img_dim**2 - missing_indices.shape[0]).to(device)
        self.missing_indices = missing_indices
        self.kept_indices = torch.Tensor([i for i in range(channels * img_dim**2) if i not in missing_indices]).to(device).long()

    def V(self, vec):
        temp = vec.clone().reshape(vec.shape[0], -1)
        out = torch.zeros_like(temp)
        out[:, self.kept_indices] = temp[:, :self.kept_indices.shape[0]]
        out[:, self.missing_indices] = temp[:, self.kept_indices.shape[0]:]
        return out.reshape(vec.shape[0], -1, self.channels).permute(0, 2, 1).reshape(vec.shape[0], -1)

    def Vt(self, vec):
        temp = vec.clone().reshape(vec.shape[0], self.channels, -1).permute(0, 2, 1).reshape(vec.shape[0], -1)
        out = torch.zeros_like(temp)
        out[:, :self.kept_indices.shape[0]] = temp[:, self.kept_indices]
        out[:, self.kept_indices.shape[0]:] = temp[:, self.missing_indices]
        return out

    def U(self, vec):
        return vec.clone().reshape(vec.shape[0], -1)

    def Ut(self, vec):
        return vec.clone().reshape(vec.shape[0], -1)

    def singulars(self):
        return self._singulars

    def add_zeros(self, vec):
        temp = torch.zeros((vec.shape[0], self.channels * self.img_dim**2), device=vec.device)
        reshaped = vec.clone().reshape(vec.shape[0], -1)
        temp[:, :reshaped.shape[1]] = reshaped
        return temp

#Denoising
class Denoising(H_functions):
    def __init__(self, channels, img_dim, device):
        self._singulars = torch.ones(channels * img_dim**2, device=device)

    def V(self, vec):
        return vec.clone().reshape(vec.shape[0], -1)

    def Vt(self, vec):
        return vec.clone().reshape(vec.shape[0], -1)

    def U(self, vec):
        return vec.clone().reshape(vec.shape[0], -1)

    def Ut(self, vec):
        return vec.clone().reshape(vec.shape[0], -1)

    def singulars(self):
        return self._singulars

    def add_zeros(self, vec):
        return vec.clone().reshape(vec.shape[0], -1)
