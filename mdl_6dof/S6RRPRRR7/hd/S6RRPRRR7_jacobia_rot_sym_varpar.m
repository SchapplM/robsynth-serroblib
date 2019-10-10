% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR7
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
%   Wie in S6RRPRRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:59
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (82->16), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->26)
	t29 = cos(qJ(1));
	t30 = t29 ^ 2;
	t28 = cos(qJ(2));
	t26 = sin(qJ(2));
	t27 = sin(qJ(1));
	t31 = t27 * t26;
	t18 = atan2(-t31, -t28);
	t16 = sin(t18);
	t17 = cos(t18);
	t14 = -t16 * t31 - t17 * t28;
	t13 = 0.1e1 / t14 ^ 2;
	t36 = t13 * t26;
	t35 = t16 * t28;
	t21 = t26 ^ 2;
	t24 = 0.1e1 / t28 ^ 2;
	t34 = t21 * t24;
	t22 = t27 ^ 2;
	t33 = t22 / t30;
	t19 = 0.1e1 / (t22 * t34 + 0.1e1);
	t32 = t27 * t19;
	t23 = 0.1e1 / t28;
	t20 = 0.1e1 / (t24 * t33 + 0.1e1);
	t15 = (0.1e1 + t34) * t32;
	t12 = 0.1e1 / t14;
	t11 = 0.1e1 / (t30 * t21 * t13 + 0.1e1);
	t1 = [t29 * t26 * t23 * t19, t15, 0, 0, 0, 0; (-t12 * t31 - (-t17 * t21 * t23 * t32 + (t19 - 0.1e1) * t26 * t16) * t30 * t36) * t11, (t28 * t12 - (-t27 * t35 + t17 * t26 + (-t17 * t31 + t35) * t15) * t36) * t29 * t11, 0, 0, 0, 0; (-0.1e1 - t33) * t23 * t20, -0.1e1 / t29 * t24 * t20 * t31, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (307->23), mult. (876->55), div. (85->9), fcn. (1282->11), ass. (0->37)
	t68 = sin(qJ(4));
	t71 = cos(qJ(2));
	t87 = sin(qJ(2));
	t88 = cos(qJ(4));
	t61 = -t71 * t68 + t87 * t88;
	t69 = sin(qJ(1));
	t53 = t61 * t69;
	t60 = t87 * t68 + t71 * t88;
	t48 = atan2(t53, t60);
	t45 = sin(t48);
	t46 = cos(t48);
	t43 = t45 * t53 + t46 * t60;
	t42 = 0.1e1 / t43 ^ 2;
	t89 = cos(qJ(1));
	t57 = t61 * t89;
	t81 = t57 ^ 2 * t42;
	t40 = 0.1e1 / (0.1e1 + t81);
	t41 = 0.1e1 / t43;
	t55 = t60 * t69;
	t56 = t60 * t89;
	t84 = t46 * t53;
	t59 = 0.1e1 / t60 ^ 2;
	t47 = 0.1e1 / (t53 ^ 2 * t59 + 0.1e1);
	t58 = 0.1e1 / t60;
	t91 = (-t53 * t59 * t61 - t55 * t58) * t47;
	t94 = (t42 * t57 * (-t45 * t55 + t46 * t61 + (-t45 * t60 + t84) * t91) + t56 * t41) * t40;
	t67 = sin(qJ(5));
	t70 = cos(qJ(5));
	t52 = t56 * t70 - t69 * t67;
	t50 = 0.1e1 / t52 ^ 2;
	t51 = t56 * t67 + t69 * t70;
	t79 = t51 ^ 2 * t50 + 0.1e1;
	t44 = 0.1e1 / t79;
	t49 = 0.1e1 / t52;
	t83 = t50 * t51;
	t90 = (-t67 * t49 + t70 * t83) * t44 * t57;
	t1 = [t57 * t58 * t47, -t91, 0, t91, 0, 0; (t53 * t41 - (-t45 + (-t58 * t84 + t45) * t47) * t81) * t40, -t94, 0, t94, 0, 0; ((-t55 * t67 + t89 * t70) * t49 - (-t55 * t70 - t89 * t67) * t83) * t44, t90, 0, -t90, t79 * t44, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (402->24), mult. (924->55), div. (90->9), fcn. (1345->11), ass. (0->39)
	t102 = sin(qJ(2));
	t103 = cos(qJ(4));
	t84 = sin(qJ(4));
	t86 = cos(qJ(2));
	t75 = t102 * t103 - t86 * t84;
	t85 = sin(qJ(1));
	t67 = t75 * t85;
	t74 = t102 * t84 + t86 * t103;
	t73 = 0.1e1 / t74 ^ 2;
	t61 = 0.1e1 / (t67 ^ 2 * t73 + 0.1e1);
	t69 = t74 * t85;
	t72 = 0.1e1 / t74;
	t106 = (-t67 * t73 * t75 - t69 * t72) * t61;
	t64 = atan2(t67, t74);
	t59 = sin(t64);
	t60 = cos(t64);
	t57 = t59 * t67 + t60 * t74;
	t56 = 0.1e1 / t57 ^ 2;
	t104 = cos(qJ(1));
	t71 = t75 * t104;
	t96 = t71 ^ 2 * t56;
	t54 = 0.1e1 / (0.1e1 + t96);
	t55 = 0.1e1 / t57;
	t70 = t74 * t104;
	t99 = t60 * t67;
	t109 = ((t106 * (-t59 * t74 + t99) - t59 * t69 + t60 * t75) * t56 * t71 + t70 * t55) * t54;
	t83 = qJ(5) + qJ(6);
	t81 = sin(t83);
	t82 = cos(t83);
	t66 = t70 * t82 - t85 * t81;
	t63 = 0.1e1 / t66 ^ 2;
	t65 = t70 * t81 + t85 * t82;
	t94 = t65 ^ 2 * t63 + 0.1e1;
	t58 = 0.1e1 / t94;
	t62 = 0.1e1 / t66;
	t98 = t63 * t65;
	t105 = (-t81 * t62 + t82 * t98) * t58 * t71;
	t51 = t94 * t58;
	t1 = [t71 * t72 * t61, -t106, 0, t106, 0, 0; (t67 * t55 - (-t59 + (-t72 * t99 + t59) * t61) * t96) * t54, -t109, 0, t109, 0, 0; ((t104 * t82 - t69 * t81) * t62 - (-t104 * t81 - t69 * t82) * t98) * t58, t105, 0, -t105, t51, t51;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end