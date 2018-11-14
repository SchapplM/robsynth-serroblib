% Analytische Jacobi-Matrix für beliebiges Segment von
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4RRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% Ja [6x4]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja = S4RRPR1_jacobia_floatb_twist_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_jacobia_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR1_jacobia_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR1_jacobia_floatb_twist_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_jacobia_floatb_twist_sym_varpar: pkin has to be [7x1] (double)');
%% Function calls
if link_index == 0
	Ja=S4RRPR1_jacobia_0_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 1
	Ja=S4RRPR1_jacobia_1_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 2
	Ja=S4RRPR1_jacobia_2_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 3
	Ja=S4RRPR1_jacobia_3_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 4
	Ja=S4RRPR1_jacobia_4_floatb_twist_sym_varpar(qJ, r_i_i_C, pkin);
else
	Ja=NaN(6,4);
end