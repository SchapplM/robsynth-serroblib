% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RPRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:57
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobia_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
	t37 = cos(qJ(3));
	t34 = sin(qJ(3));
	t35 = sin(qJ(1));
	t43 = t35 * t34;
	t28 = atan2(-t43, -t37);
	t26 = sin(t28);
	t27 = cos(t28);
	t19 = -t26 * t43 - t27 * t37;
	t18 = 0.1e1 / t19 ^ 2;
	t38 = cos(qJ(1));
	t48 = t18 * t38 ^ 2;
	t33 = sin(qJ(4));
	t36 = cos(qJ(4));
	t40 = t38 * t36;
	t25 = t35 * t33 + t37 * t40;
	t23 = 0.1e1 / t25 ^ 2;
	t41 = t38 * t33;
	t24 = -t35 * t36 + t37 * t41;
	t47 = t23 * t24;
	t30 = t34 ^ 2;
	t46 = t30 / t37 ^ 2;
	t45 = t34 * t38;
	t29 = 0.1e1 / (t35 ^ 2 * t46 + 0.1e1);
	t44 = t35 * t29;
	t42 = t35 * t37;
	t39 = t24 ^ 2 * t23 + 0.1e1;
	t31 = 0.1e1 / t37;
	t22 = 0.1e1 / t25;
	t21 = (0.1e1 + t46) * t44;
	t20 = 0.1e1 / t39;
	t17 = 0.1e1 / t19;
	t16 = 0.1e1 / (t30 * t48 + 0.1e1);
	t1 = [t31 * t29 * t45, 0, t21, 0, 0; (-t17 * t43 - (-t27 * t30 * t31 * t44 + (t29 - 0.1e1) * t34 * t26) * t34 * t48) * t16, 0, (t37 * t17 - (-t26 * t42 + t27 * t34 + (t26 * t37 - t27 * t43) * t21) * t34 * t18) * t38 * t16, 0, 0; ((-t33 * t42 - t40) * t22 - (-t36 * t42 + t41) * t47) * t20, 0, (-t22 * t33 + t36 * t47) * t20 * t45, t39 * t20, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (216->32), mult. (662->85), div. (108->11), fcn. (985->11), ass. (0->45)
	t59 = sin(qJ(4));
	t63 = cos(qJ(4));
	t65 = cos(qJ(1));
	t67 = t65 * t63;
	t61 = sin(qJ(1));
	t64 = cos(qJ(3));
	t69 = t61 * t64;
	t47 = t59 * t69 + t67;
	t60 = sin(qJ(3));
	t73 = t60 * t59;
	t46 = atan2(-t47, t73);
	t43 = sin(t46);
	t44 = cos(t46);
	t37 = -t43 * t47 + t44 * t73;
	t36 = 0.1e1 / t37 ^ 2;
	t68 = t65 * t59;
	t50 = -t61 * t63 + t64 * t68;
	t78 = t36 * t50;
	t51 = t61 * t59 + t64 * t67;
	t58 = sin(qJ(5));
	t62 = cos(qJ(5));
	t70 = t60 * t65;
	t42 = t51 * t62 + t58 * t70;
	t40 = 0.1e1 / t42 ^ 2;
	t41 = t51 * t58 - t62 * t70;
	t77 = t40 * t41;
	t76 = t44 * t47;
	t75 = t50 ^ 2 * t36;
	t54 = 0.1e1 / t59;
	t56 = 0.1e1 / t60;
	t74 = t54 * t56;
	t72 = t60 * t61;
	t71 = t60 * t63;
	t66 = t41 ^ 2 * t40 + 0.1e1;
	t57 = 0.1e1 / t60 ^ 2;
	t55 = 0.1e1 / t59 ^ 2;
	t49 = t63 * t69 - t68;
	t45 = 0.1e1 / (t47 ^ 2 * t57 * t55 + 0.1e1);
	t39 = 0.1e1 / t42;
	t38 = 0.1e1 / t66;
	t35 = 0.1e1 / t37;
	t34 = (t47 * t54 * t57 * t64 + t61) * t45;
	t33 = 0.1e1 / (0.1e1 + t75);
	t32 = (t47 * t55 * t63 - t49 * t54) * t56 * t45;
	t1 = [-t50 * t45 * t74, 0, t34, t32, 0; (-t47 * t35 - (-t43 + (t74 * t76 + t43) * t45) * t75) * t33, 0, (t34 * t76 * t78 + (-t35 * t70 - (t44 * t64 + (-t34 * t60 + t72) * t43) * t78) * t59) * t33, (t51 * t35 - (t44 * t71 - t43 * t49 + (-t43 * t73 - t76) * t32) * t78) * t33, 0; ((-t49 * t58 + t62 * t72) * t39 - (-t49 * t62 - t58 * t72) * t77) * t38, 0, ((-t58 * t71 - t62 * t64) * t39 - (t58 * t64 - t62 * t71) * t77) * t38 * t65, (-t58 * t39 + t62 * t77) * t50 * t38, t66 * t38;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end