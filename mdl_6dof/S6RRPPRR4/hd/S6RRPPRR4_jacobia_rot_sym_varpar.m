% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR4
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
%   Wie in S6RRPPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:41
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
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
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (73->16), mult. (212->41), div. (26->9), fcn. (311->11), ass. (0->28)
	t48 = cos(pkin(6));
	t45 = sin(pkin(11));
	t47 = cos(pkin(11));
	t49 = sin(qJ(2));
	t51 = cos(qJ(2));
	t54 = t51 * t45 + t49 * t47;
	t35 = t54 * t48;
	t36 = t49 * t45 - t51 * t47;
	t50 = sin(qJ(1));
	t52 = cos(qJ(1));
	t31 = -t50 * t35 - t52 * t36;
	t28 = 0.1e1 / t31 ^ 2;
	t53 = t36 * t48;
	t29 = t50 * t53 - t52 * t54;
	t58 = t29 ^ 2 * t28;
	t46 = sin(pkin(6));
	t55 = t52 * t46;
	t41 = atan2(t55, t48);
	t38 = sin(t41);
	t39 = cos(t41);
	t33 = t38 * t55 + t39 * t48;
	t57 = 0.1e1 / t33 ^ 2 * t50 ^ 2;
	t43 = t46 ^ 2;
	t40 = 0.1e1 / (0.1e1 + t52 ^ 2 * t43 / t48 ^ 2);
	t56 = t40 / t48;
	t27 = 0.1e1 / t31;
	t25 = 0.1e1 / (0.1e1 + t58);
	t1 = [-t50 * t46 * t56, 0, 0, 0, 0, 0; (0.1e1 / t33 * t55 - (-t39 * t43 * t52 * t56 + (t40 - 0.1e1) * t46 * t38) * t46 * t57) / (t43 * t57 + 0.1e1), 0, 0, 0, 0, 0; ((-t50 * t54 - t52 * t53) * t27 + (-t52 * t35 + t50 * t36) * t29 * t28) * t25, (t31 * t27 + t58) * t25, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (289->21), mult. (821->56), div. (53->11), fcn. (1181->11), ass. (0->37)
	t64 = sin(qJ(1));
	t66 = cos(qJ(1));
	t60 = sin(pkin(11));
	t62 = cos(pkin(11));
	t65 = cos(qJ(2));
	t71 = cos(pkin(6));
	t69 = t65 * t71;
	t63 = sin(qJ(2));
	t70 = t63 * t71;
	t67 = -t60 * t70 + t62 * t69;
	t68 = t65 * t60 + t63 * t62;
	t44 = -t64 * t68 + t66 * t67;
	t54 = t63 * t60 - t65 * t62;
	t61 = sin(pkin(6));
	t51 = t54 * t61;
	t41 = atan2(t44, t51);
	t39 = cos(t41);
	t74 = t39 * t44;
	t53 = t60 * t69 + t62 * t70;
	t47 = -t64 * t53 - t66 * t54;
	t59 = 0.1e1 / t64 ^ 2;
	t73 = 0.1e1 / (0.1e1 + t47 ^ 2 * t59 / t61 ^ 2) / t61;
	t38 = sin(t41);
	t37 = t38 * t44 + t39 * t51;
	t36 = 0.1e1 / t37 ^ 2;
	t45 = -t64 * t67 - t66 * t68;
	t72 = t45 ^ 2 * t36;
	t43 = -t66 * t53 + t64 * t54;
	t58 = 0.1e1 / t64;
	t52 = t68 * t61;
	t50 = 0.1e1 / t51 ^ 2;
	t49 = 0.1e1 / t51;
	t40 = 0.1e1 / (t44 ^ 2 * t50 + 0.1e1);
	t35 = 0.1e1 / t37;
	t34 = 0.1e1 / (0.1e1 + t72);
	t33 = (-t44 * t50 * t52 + t43 * t49) * t40;
	t1 = [t45 * t49 * t40, t33, 0, 0, 0, 0; (t44 * t35 + (t38 + (t49 * t74 - t38) * t40) * t72) * t34, (t47 * t35 + (t38 * t43 + t39 * t52 + (-t38 * t51 + t74) * t33) * t45 * t36) * t34, 0, 0, 0, 0; (-t66 * t47 * t59 + t43 * t58) * t73, t45 * t58 * t73, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (374->26), mult. (1051->65), div. (55->9), fcn. (1490->13), ass. (0->43)
	t72 = sin(pkin(11));
	t74 = cos(pkin(11));
	t77 = sin(qJ(2));
	t80 = cos(qJ(2));
	t69 = t77 * t72 - t80 * t74;
	t78 = sin(qJ(1));
	t81 = cos(qJ(1));
	t75 = cos(pkin(6));
	t84 = t80 * t72 + t77 * t74;
	t83 = t84 * t75;
	t62 = -t81 * t69 - t78 * t83;
	t58 = -t78 * t69 + t81 * t83;
	t73 = sin(pkin(6));
	t68 = t84 * t73;
	t52 = atan2(-t58, t68);
	t50 = cos(t52);
	t91 = t50 * t58;
	t82 = t69 * t75;
	t61 = t78 * t82 - t81 * t84;
	t76 = sin(qJ(5));
	t79 = cos(qJ(5));
	t88 = t73 * t78;
	t56 = -t61 * t76 + t79 * t88;
	t54 = 0.1e1 / t56 ^ 2;
	t55 = t61 * t79 + t76 * t88;
	t90 = t54 * t55;
	t49 = sin(t52);
	t47 = -t49 * t58 + t50 * t68;
	t46 = 0.1e1 / t47 ^ 2;
	t89 = t62 ^ 2 * t46;
	t87 = t73 * t81;
	t85 = t55 ^ 2 * t54 + 0.1e1;
	t67 = t69 * t73;
	t66 = 0.1e1 / t68 ^ 2;
	t65 = 0.1e1 / t68;
	t57 = -t78 * t84 - t81 * t82;
	t53 = 0.1e1 / t56;
	t51 = 0.1e1 / (t58 ^ 2 * t66 + 0.1e1);
	t48 = 0.1e1 / t85;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (0.1e1 + t89);
	t43 = (-t58 * t66 * t67 - t57 * t65) * t51;
	t1 = [-t62 * t65 * t51, t43, 0, 0, 0, 0; (-t58 * t45 - (-t49 + (t65 * t91 + t49) * t51) * t89) * t44, (t61 * t45 - (-t49 * t57 - t50 * t67 + (-t49 * t68 - t91) * t43) * t62 * t46) * t44, 0, 0, 0, 0; ((-t57 * t79 + t76 * t87) * t53 - (t57 * t76 + t79 * t87) * t90) * t48, (-t79 * t53 - t76 * t90) * t62 * t48, 0, 0, t85 * t48, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:41:13
	% EndTime: 2019-10-10 09:41:13
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (867->38), mult. (2363->97), div. (85->9), fcn. (3329->15), ass. (0->56)
	t105 = sin(qJ(5));
	t109 = cos(qJ(5));
	t101 = sin(pkin(6));
	t111 = cos(qJ(1));
	t119 = t101 * t111;
	t107 = sin(qJ(1));
	t103 = cos(pkin(6));
	t100 = sin(pkin(11));
	t102 = cos(pkin(11));
	t106 = sin(qJ(2));
	t110 = cos(qJ(2));
	t98 = t106 * t100 - t110 * t102;
	t113 = t98 * t103;
	t115 = t110 * t100 + t106 * t102;
	t86 = -t107 * t115 - t111 * t113;
	t81 = t105 * t119 - t86 * t109;
	t94 = t98 * t101;
	t92 = t103 * t105 - t94 * t109;
	t78 = atan2(t81, t92);
	t75 = sin(t78);
	t76 = cos(t78);
	t69 = t75 * t81 + t76 * t92;
	t68 = 0.1e1 / t69 ^ 2;
	t112 = t107 * t113 - t111 * t115;
	t120 = t101 * t107;
	t79 = t105 * t120 + t112 * t109;
	t126 = t68 * t79;
	t104 = sin(qJ(6));
	t108 = cos(qJ(6));
	t96 = t115 * t103;
	t116 = -t107 * t96 - t111 * t98;
	t80 = -t112 * t105 + t109 * t120;
	t74 = t104 * t116 + t80 * t108;
	t72 = 0.1e1 / t74 ^ 2;
	t73 = t80 * t104 - t108 * t116;
	t125 = t72 * t73;
	t124 = t76 * t81;
	t123 = t79 ^ 2 * t68;
	t91 = 0.1e1 / t92 ^ 2;
	t122 = t81 * t91;
	t121 = t116 * t105;
	t118 = t73 ^ 2 * t72 + 0.1e1;
	t117 = -t75 * t92 + t124;
	t87 = t107 * t98 - t111 * t96;
	t114 = t86 * t105 + t109 * t119;
	t95 = t115 * t101;
	t93 = t103 * t109 + t94 * t105;
	t90 = 0.1e1 / t92;
	t77 = 0.1e1 / (t81 ^ 2 * t91 + 0.1e1);
	t71 = 0.1e1 / t74;
	t70 = 0.1e1 / t118;
	t67 = 0.1e1 / t69;
	t66 = 0.1e1 / (0.1e1 + t123);
	t65 = (t95 * t122 - t87 * t90) * t77 * t109;
	t64 = (t114 * t90 - t93 * t122) * t77;
	t1 = [-t79 * t90 * t77, t65, 0, 0, t64, 0; (t81 * t67 - (-t75 + (-t90 * t124 + t75) * t77) * t123) * t66, (-t116 * t109 * t67 - (t117 * t65 + (-t75 * t87 - t76 * t95) * t109) * t126) * t66, 0, 0, (t80 * t67 - (t114 * t75 + t117 * t64 + t76 * t93) * t126) * t66, 0; ((t104 * t114 - t87 * t108) * t71 - (t87 * t104 + t108 * t114) * t125) * t70, ((t104 * t121 - t112 * t108) * t71 - (t112 * t104 + t108 * t121) * t125) * t70, 0, 0, (-t104 * t71 + t108 * t125) * t79 * t70, t118 * t70;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end