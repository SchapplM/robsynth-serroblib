% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR7
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
%   Wie in S6RRRPRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->13), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->25)
	t30 = cos(pkin(6));
	t29 = sin(pkin(6));
	t34 = cos(qJ(1));
	t38 = t34 * t29;
	t26 = atan2(t38, t30);
	t23 = sin(t26);
	t24 = cos(t26);
	t18 = t23 * t38 + t24 * t30;
	t32 = sin(qJ(1));
	t42 = 0.1e1 / t18 ^ 2 * t32 ^ 2;
	t27 = t29 ^ 2;
	t25 = 0.1e1 / (0.1e1 + t34 ^ 2 * t27 / t30 ^ 2);
	t41 = t25 / t30;
	t31 = sin(qJ(2));
	t40 = t32 * t31;
	t33 = cos(qJ(2));
	t39 = t32 * t33;
	t37 = t34 * t31;
	t36 = t34 * t33;
	t22 = -t30 * t40 + t36;
	t20 = 0.1e1 / t22 ^ 2;
	t21 = t30 * t39 + t37;
	t35 = t21 ^ 2 * t20 + 0.1e1;
	t19 = 0.1e1 / t35;
	t1 = [-t32 * t29 * t41, 0, 0, 0, 0, 0; (0.1e1 / t18 * t38 - (-t24 * t27 * t34 * t41 + (t25 - 0.1e1) * t29 * t23) * t29 * t42) / (t27 * t42 + 0.1e1), 0, 0, 0, 0, 0; ((t30 * t36 - t40) / t22 - (-t30 * t37 - t39) * t21 * t20) * t19, t35 * t19, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (173->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
	t53 = sin(qJ(2));
	t56 = cos(qJ(2));
	t57 = cos(qJ(1));
	t54 = sin(qJ(1));
	t61 = cos(pkin(6));
	t59 = t54 * t61;
	t46 = -t53 * t59 + t57 * t56;
	t52 = sin(qJ(3));
	t55 = cos(qJ(3));
	t51 = sin(pkin(6));
	t64 = t51 * t54;
	t37 = t46 * t55 + t52 * t64;
	t35 = 0.1e1 / t37 ^ 2;
	t36 = t46 * t52 - t55 * t64;
	t68 = t35 * t36;
	t58 = t57 * t61;
	t42 = t54 * t53 - t56 * t58;
	t63 = t51 * t56;
	t40 = atan2(-t42, -t63);
	t39 = cos(t40);
	t67 = t39 * t42;
	t38 = sin(t40);
	t32 = -t38 * t42 - t39 * t63;
	t31 = 0.1e1 / t32 ^ 2;
	t45 = t57 * t53 + t56 * t59;
	t66 = t45 ^ 2 * t31;
	t48 = 0.1e1 / t51;
	t49 = 0.1e1 / t56;
	t65 = t48 * t49;
	t62 = t51 * t57;
	t60 = t36 ^ 2 * t35 + 0.1e1;
	t50 = 0.1e1 / t56 ^ 2;
	t44 = t53 * t58 + t54 * t56;
	t41 = 0.1e1 / (0.1e1 + t42 ^ 2 / t51 ^ 2 * t50);
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t60;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (0.1e1 + t66);
	t28 = (t42 * t50 * t53 + t44 * t49) * t48 * t41;
	t1 = [t45 * t41 * t65, t28, 0, 0, 0, 0; (-t42 * t30 - (-t38 + (-t65 * t67 + t38) * t41) * t66) * t29, (t46 * t30 - (t39 * t51 * t53 - t38 * t44 + (t38 * t63 - t67) * t28) * t45 * t31) * t29, 0, 0, 0, 0; ((-t44 * t52 - t55 * t62) * t34 - (-t44 * t55 + t52 * t62) * t68) * t33, (-t52 * t34 + t55 * t68) * t45 * t33, t60 * t33, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (221->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
	t59 = sin(qJ(2));
	t61 = cos(qJ(2));
	t62 = cos(qJ(1));
	t60 = sin(qJ(1));
	t66 = cos(pkin(6));
	t64 = t60 * t66;
	t50 = -t59 * t64 + t62 * t61;
	t55 = qJ(3) + pkin(12);
	t52 = sin(t55);
	t53 = cos(t55);
	t58 = sin(pkin(6));
	t69 = t58 * t60;
	t41 = t50 * t53 + t52 * t69;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = t50 * t52 - t53 * t69;
	t73 = t39 * t40;
	t63 = t62 * t66;
	t46 = t60 * t59 - t61 * t63;
	t68 = t58 * t61;
	t44 = atan2(-t46, -t68);
	t43 = cos(t44);
	t72 = t43 * t46;
	t42 = sin(t44);
	t36 = -t42 * t46 - t43 * t68;
	t35 = 0.1e1 / t36 ^ 2;
	t49 = t62 * t59 + t61 * t64;
	t71 = t49 ^ 2 * t35;
	t54 = 0.1e1 / t58;
	t56 = 0.1e1 / t61;
	t70 = t54 * t56;
	t67 = t58 * t62;
	t65 = t40 ^ 2 * t39 + 0.1e1;
	t57 = 0.1e1 / t61 ^ 2;
	t48 = t59 * t63 + t60 * t61;
	t45 = 0.1e1 / (0.1e1 + t46 ^ 2 / t58 ^ 2 * t57);
	t38 = 0.1e1 / t41;
	t37 = 0.1e1 / t65;
	t34 = 0.1e1 / t36;
	t33 = 0.1e1 / (0.1e1 + t71);
	t32 = (t46 * t57 * t59 + t48 * t56) * t54 * t45;
	t1 = [t49 * t45 * t70, t32, 0, 0, 0, 0; (-t46 * t34 - (-t42 + (-t70 * t72 + t42) * t45) * t71) * t33, (t50 * t34 - (t43 * t58 * t59 - t42 * t48 + (t42 * t68 - t72) * t32) * t49 * t35) * t33, 0, 0, 0, 0; ((-t48 * t52 - t53 * t67) * t38 - (-t48 * t53 + t52 * t67) * t73) * t37, (-t52 * t38 + t53 * t73) * t49 * t37, t65 * t37, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (314->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
	t72 = sin(qJ(2));
	t74 = cos(qJ(2));
	t75 = cos(qJ(1));
	t73 = sin(qJ(1));
	t79 = cos(pkin(6));
	t77 = t73 * t79;
	t63 = -t72 * t77 + t75 * t74;
	t67 = qJ(3) + pkin(12) + qJ(5);
	t65 = sin(t67);
	t66 = cos(t67);
	t71 = sin(pkin(6));
	t82 = t71 * t73;
	t54 = t63 * t66 + t65 * t82;
	t52 = 0.1e1 / t54 ^ 2;
	t53 = t63 * t65 - t66 * t82;
	t86 = t52 * t53;
	t76 = t75 * t79;
	t59 = t73 * t72 - t74 * t76;
	t81 = t71 * t74;
	t57 = atan2(-t59, -t81);
	t56 = cos(t57);
	t85 = t56 * t59;
	t55 = sin(t57);
	t49 = -t55 * t59 - t56 * t81;
	t48 = 0.1e1 / t49 ^ 2;
	t62 = t75 * t72 + t74 * t77;
	t84 = t62 ^ 2 * t48;
	t68 = 0.1e1 / t71;
	t69 = 0.1e1 / t74;
	t83 = t68 * t69;
	t80 = t71 * t75;
	t78 = t53 ^ 2 * t52 + 0.1e1;
	t70 = 0.1e1 / t74 ^ 2;
	t61 = t72 * t76 + t73 * t74;
	t58 = 0.1e1 / (0.1e1 + t59 ^ 2 / t71 ^ 2 * t70);
	t51 = 0.1e1 / t54;
	t50 = 0.1e1 / t78;
	t47 = 0.1e1 / t49;
	t46 = 0.1e1 / (0.1e1 + t84);
	t45 = (t59 * t70 * t72 + t61 * t69) * t68 * t58;
	t44 = t78 * t50;
	t1 = [t62 * t58 * t83, t45, 0, 0, 0, 0; (-t59 * t47 - (-t55 + (-t83 * t85 + t55) * t58) * t84) * t46, (t63 * t47 - (t56 * t71 * t72 - t55 * t61 + (t55 * t81 - t85) * t45) * t62 * t48) * t46, 0, 0, 0, 0; ((-t61 * t65 - t66 * t80) * t51 - (-t61 * t66 + t65 * t80) * t86) * t50, (-t65 * t51 + t66 * t86) * t62 * t50, t44, 0, t44, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:04:12
	% EndTime: 2019-10-10 12:04:12
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (1756->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
	t102 = qJ(3) + pkin(12) + qJ(5);
	t100 = sin(t102);
	t101 = cos(t102);
	t103 = sin(pkin(6));
	t110 = cos(qJ(1));
	t117 = t103 * t110;
	t104 = cos(pkin(6));
	t106 = sin(qJ(2));
	t114 = t110 * t106;
	t107 = sin(qJ(1));
	t109 = cos(qJ(2));
	t115 = t107 * t109;
	t95 = t104 * t114 + t115;
	t84 = t95 * t100 + t101 * t117;
	t120 = t103 * t106;
	t92 = t100 * t120 - t104 * t101;
	t79 = atan2(-t84, t92);
	t76 = sin(t79);
	t77 = cos(t79);
	t74 = -t76 * t84 + t77 * t92;
	t73 = 0.1e1 / t74 ^ 2;
	t119 = t103 * t107;
	t113 = t110 * t109;
	t116 = t107 * t106;
	t97 = -t104 * t116 + t113;
	t88 = t97 * t100 - t101 * t119;
	t127 = t73 * t88;
	t126 = t77 * t84;
	t108 = cos(qJ(6));
	t105 = sin(qJ(6));
	t96 = t104 * t115 + t114;
	t122 = t96 * t105;
	t89 = t100 * t119 + t97 * t101;
	t83 = t89 * t108 + t122;
	t81 = 0.1e1 / t83 ^ 2;
	t121 = t96 * t108;
	t82 = t89 * t105 - t121;
	t125 = t81 * t82;
	t91 = 0.1e1 / t92 ^ 2;
	t124 = t84 * t91;
	t123 = t88 ^ 2 * t73;
	t118 = t103 * t109;
	t112 = t82 ^ 2 * t81 + 0.1e1;
	t86 = -t100 * t117 + t95 * t101;
	t111 = -t76 * t92 - t126;
	t94 = t104 * t113 - t116;
	t93 = t104 * t100 + t101 * t120;
	t90 = 0.1e1 / t92;
	t80 = 0.1e1 / t83;
	t78 = 0.1e1 / (t84 ^ 2 * t91 + 0.1e1);
	t75 = 0.1e1 / t112;
	t72 = 0.1e1 / t74;
	t71 = 0.1e1 / (0.1e1 + t123);
	t70 = (t118 * t124 - t90 * t94) * t78 * t100;
	t69 = (t93 * t124 - t86 * t90) * t78;
	t68 = (-t105 * t80 + t108 * t125) * t88 * t75;
	t67 = (t89 * t72 - (t111 * t69 - t76 * t86 + t77 * t93) * t127) * t71;
	t1 = [-t88 * t90 * t78, t70, t69, 0, t69, 0; (-t84 * t72 - (-t76 + (t90 * t126 + t76) * t78) * t123) * t71, (-t96 * t100 * t72 - (t111 * t70 + (t77 * t118 - t76 * t94) * t100) * t127) * t71, t67, 0, t67, 0; ((-t105 * t86 - t94 * t108) * t80 - (t94 * t105 - t108 * t86) * t125) * t75, ((-t101 * t122 - t97 * t108) * t80 - (-t101 * t121 + t97 * t105) * t125) * t75, t68, 0, t68, t112 * t75;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end