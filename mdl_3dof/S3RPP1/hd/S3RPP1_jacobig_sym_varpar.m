% Geometrischen Jacobi-Matrix für beliebiges Segment von
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% 
% Output:
% Jg [6x3]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = S3RPP1_jacobig_sym_varpar(qJ, link_index, r_i_i_C, pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),uint8(0),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_jacobig_sym_varpar: qJ has to be [3x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S3RPP1_jacobig_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S3RPP1_jacobig_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_jacobig_sym_varpar: pkin has to be [3x1] (double)');
%% Function calls
if link_index == 0
	Jg=S3RPP1_jacobig_0_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 1
	Jg=S3RPP1_jacobig_1_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 2
	Jg=S3RPP1_jacobig_2_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 3
	Jg=S3RPP1_jacobig_3_sym_varpar(qJ, r_i_i_C, pkin);
else
	Jg=NaN(6,3);
end