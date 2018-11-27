% Geometrischen Jacobi-Matrix für beliebiges Segment von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% Jg [6x7]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg = S7RRRRRRR1_jacobig_floatb_twist_sym_varpar(qJ, link_index, r_i_i_C, pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(7,1),uint8(0),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobig_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobig_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobig_floatb_twist_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobig_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');
%% Function calls
if link_index == 0
	Jg=S7RRRRRRR1_jacobig_0_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 1
	Jg=S7RRRRRRR1_jacobig_1_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 2
	Jg=S7RRRRRRR1_jacobig_2_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 3
	Jg=S7RRRRRRR1_jacobig_3_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 4
	Jg=S7RRRRRRR1_jacobig_4_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 5
	Jg=S7RRRRRRR1_jacobig_5_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 6
	Jg=S7RRRRRRR1_jacobig_6_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 7
	Jg=S7RRRRRRR1_jacobig_7_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
else
	Jg=NaN(6,7);
end