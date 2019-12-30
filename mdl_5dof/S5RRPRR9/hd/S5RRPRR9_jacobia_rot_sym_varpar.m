% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRR9
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
%   Wie in S5RRPRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPRR9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:32
	% EndTime: 2019-12-29 19:10:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:32
	% EndTime: 2019-12-29 19:10:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:32
	% EndTime: 2019-12-29 19:10:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:27
	% EndTime: 2019-12-29 19:10:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:27
	% EndTime: 2019-12-29 19:10:27
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (238->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t46 = qJ(2) + pkin(9);
	t45 = cos(t46);
	t44 = sin(t46);
	t48 = sin(qJ(1));
	t56 = t48 * t44;
	t39 = atan2(-t56, -t45);
	t35 = sin(t39);
	t36 = cos(t39);
	t30 = -t35 * t56 - t36 * t45;
	t29 = 0.1e1 / t30 ^ 2;
	t50 = cos(qJ(1));
	t62 = t29 * t50 ^ 2;
	t49 = cos(qJ(4));
	t52 = t50 * t49;
	t47 = sin(qJ(4));
	t55 = t48 * t47;
	t38 = t45 * t52 + t55;
	t34 = 0.1e1 / t38 ^ 2;
	t53 = t50 * t47;
	t54 = t48 * t49;
	t37 = t45 * t53 - t54;
	t61 = t34 * t37;
	t60 = t35 * t45;
	t41 = t44 ^ 2;
	t59 = t41 / t45 ^ 2;
	t58 = t44 * t50;
	t40 = 0.1e1 / (t48 ^ 2 * t59 + 0.1e1);
	t57 = t48 * t40;
	t51 = t37 ^ 2 * t34 + 0.1e1;
	t42 = 0.1e1 / t45;
	t33 = 0.1e1 / t38;
	t32 = 0.1e1 / t51;
	t31 = (0.1e1 + t59) * t57;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / (t41 * t62 + 0.1e1);
	t1 = [t42 * t40 * t58, t31, 0, 0, 0; (-t28 * t56 - (-t36 * t41 * t42 * t57 + (t40 - 0.1e1) * t44 * t35) * t44 * t62) * t27, (t45 * t28 - (-t48 * t60 + t36 * t44 + (-t36 * t56 + t60) * t31) * t44 * t29) * t50 * t27, 0, 0, 0; ((-t45 * t55 - t52) * t33 - (-t45 * t54 + t53) * t61) * t32, (-t33 * t47 + t49 * t61) * t32 * t58, 0, t51 * t32, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:27
	% EndTime: 2019-12-29 19:10:27
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (317->22), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
	t64 = qJ(2) + pkin(9);
	t61 = cos(t64);
	t60 = sin(t64);
	t66 = sin(qJ(1));
	t73 = t66 * t60;
	t55 = atan2(-t73, -t61);
	t53 = sin(t55);
	t54 = cos(t55);
	t46 = -t53 * t73 - t54 * t61;
	t45 = 0.1e1 / t46 ^ 2;
	t67 = cos(qJ(1));
	t79 = t45 * t67 ^ 2;
	t65 = qJ(4) + qJ(5);
	t63 = cos(t65);
	t69 = t67 * t63;
	t62 = sin(t65);
	t72 = t66 * t62;
	t52 = t61 * t69 + t72;
	t50 = 0.1e1 / t52 ^ 2;
	t70 = t67 * t62;
	t71 = t66 * t63;
	t51 = t61 * t70 - t71;
	t78 = t50 * t51;
	t77 = t53 * t61;
	t57 = t60 ^ 2;
	t76 = t57 / t61 ^ 2;
	t75 = t60 * t67;
	t56 = 0.1e1 / (t66 ^ 2 * t76 + 0.1e1);
	t74 = t66 * t56;
	t68 = t51 ^ 2 * t50 + 0.1e1;
	t58 = 0.1e1 / t61;
	t49 = 0.1e1 / t52;
	t48 = 0.1e1 / t68;
	t47 = (0.1e1 + t76) * t74;
	t44 = 0.1e1 / t46;
	t43 = 0.1e1 / (t57 * t79 + 0.1e1);
	t42 = t68 * t48;
	t1 = [t58 * t56 * t75, t47, 0, 0, 0; (-t44 * t73 - (-t54 * t57 * t58 * t74 + (t56 - 0.1e1) * t60 * t53) * t60 * t79) * t43, (t61 * t44 - (-t66 * t77 + t54 * t60 + (-t54 * t73 + t77) * t47) * t60 * t45) * t67 * t43, 0, 0, 0; ((-t61 * t72 - t69) * t49 - (-t61 * t71 + t70) * t78) * t48, (-t49 * t62 + t63 * t78) * t48 * t75, 0, t42, t42;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end