% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4RRPR1_jacobia_rot_floatb_twist_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_jacobia_rot_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR1_jacobia_rot_floatb_twist_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_jacobia_rot_floatb_twist_sym_varpar: pkin has to be [7x1] (double)');
%% Function calls
if link_index == 0
	Ja_rot=S4RRPR1_jacobia_rot_0_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 1
	Ja_rot=S4RRPR1_jacobia_rot_1_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 2
	Ja_rot=S4RRPR1_jacobia_rot_2_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 3
	Ja_rot=S4RRPR1_jacobia_rot_3_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 4
	Ja_rot=S4RRPR1_jacobia_rot_4_floatb_twist_sym_varpar(qJ, pkin);
else
	Ja_rot=NaN(3,4);
end