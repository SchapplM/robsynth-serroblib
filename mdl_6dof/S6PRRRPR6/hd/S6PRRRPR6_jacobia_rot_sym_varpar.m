% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PRRRPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR6_jacobia_rot_floatb_twist_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobia_rot_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR6_jacobia_rot_floatb_twist_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobia_rot_floatb_twist_sym_varpar: pkin has to be [11x1] (double)');
%% Function calls
if link_index == 0
	Ja_rot=S6PRRRPR6_jacobia_rot_0_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 1
	Ja_rot=S6PRRRPR6_jacobia_rot_1_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 2
	Ja_rot=S6PRRRPR6_jacobia_rot_2_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 3
	Ja_rot=S6PRRRPR6_jacobia_rot_3_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 4
	Ja_rot=S6PRRRPR6_jacobia_rot_4_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 5
	Ja_rot=S6PRRRPR6_jacobia_rot_5_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 6
	Ja_rot=S6PRRRPR6_jacobia_rot_6_floatb_twist_sym_varpar(qJ, pkin);
else
	Ja_rot=NaN(3,6);
end