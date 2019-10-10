% Geometrische Jacobi-Matrix für beliebiges Segment von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
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
%   Geometrische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 17:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = S7RRRRRRR1_jacobig_sym_varpar(qJ, link_index, r_i_i_C, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),uint8(0),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobig_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobig_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobig_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobig_sym_varpar: pkin has to be [4x1] (double)');

% Function calls
Ja_transl = S7RRRRRRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin);
Jg_rot = S7RRRRRRR1_jacobig_rot_sym_varpar(qJ, link_index, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
end