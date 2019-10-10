% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR5
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
%   Wie in S6RRPRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:57
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
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
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (73->16), mult. (212->41), div. (26->9), fcn. (311->11), ass. (0->28)
	t48 = cos(pkin(6));
	t45 = sin(pkin(12));
	t47 = cos(pkin(12));
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
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (374->26), mult. (1051->69), div. (55->9), fcn. (1490->13), ass. (0->45)
	t73 = sin(qJ(1));
	t76 = cos(qJ(1));
	t68 = sin(pkin(12));
	t70 = cos(pkin(12));
	t75 = cos(qJ(2));
	t83 = cos(pkin(6));
	t80 = t75 * t83;
	t72 = sin(qJ(2));
	t81 = t72 * t83;
	t77 = -t68 * t81 + t70 * t80;
	t78 = t75 * t68 + t72 * t70;
	t54 = -t73 * t78 + t76 * t77;
	t65 = t72 * t68 - t75 * t70;
	t69 = sin(pkin(6));
	t62 = t65 * t69;
	t48 = atan2(t54, t62);
	t46 = cos(t48);
	t88 = t46 * t54;
	t64 = t68 * t80 + t70 * t81;
	t58 = -t73 * t64 - t76 * t65;
	t71 = sin(qJ(4));
	t74 = cos(qJ(4));
	t85 = t69 * t73;
	t52 = t58 * t74 + t71 * t85;
	t50 = 0.1e1 / t52 ^ 2;
	t51 = t58 * t71 - t74 * t85;
	t87 = t50 * t51;
	t45 = sin(t48);
	t43 = t45 * t54 + t46 * t62;
	t42 = 0.1e1 / t43 ^ 2;
	t56 = -t73 * t77 - t76 * t78;
	t86 = t56 ^ 2 * t42;
	t84 = t69 * t76;
	t82 = t51 ^ 2 * t50 + 0.1e1;
	t79 = -t76 * t64 + t73 * t65;
	t63 = t78 * t69;
	t61 = 0.1e1 / t62 ^ 2;
	t60 = 0.1e1 / t62;
	t49 = 0.1e1 / t52;
	t47 = 0.1e1 / (t54 ^ 2 * t61 + 0.1e1);
	t44 = 0.1e1 / t82;
	t41 = 0.1e1 / t43;
	t40 = 0.1e1 / (0.1e1 + t86);
	t39 = (-t54 * t61 * t63 + t60 * t79) * t47;
	t1 = [t56 * t60 * t47, t39, 0, 0, 0, 0; (t54 * t41 + (t45 + (t60 * t88 - t45) * t47) * t86) * t40, (t58 * t41 + (t45 * t79 + t46 * t63 + (-t45 * t62 + t88) * t39) * t56 * t42) * t40, 0, 0, 0, 0; ((t71 * t79 - t74 * t84) * t49 - (t71 * t84 + t74 * t79) * t87) * t44, (t71 * t49 - t74 * t87) * t56 * t44, 0, t82 * t44, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (867->39), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->57)
	t104 = sin(qJ(4));
	t108 = cos(qJ(4));
	t100 = sin(pkin(6));
	t110 = cos(qJ(1));
	t117 = t100 * t110;
	t106 = sin(qJ(1));
	t101 = cos(pkin(12));
	t105 = sin(qJ(2));
	t109 = cos(qJ(2));
	t99 = sin(pkin(12));
	t113 = t109 * t101 - t105 * t99;
	t102 = cos(pkin(6));
	t112 = t105 * t101 + t109 * t99;
	t93 = t112 * t102;
	t82 = t106 * t113 + t110 * t93;
	t75 = t82 * t104 + t108 * t117;
	t92 = t112 * t100;
	t88 = -t102 * t108 + t92 * t104;
	t74 = atan2(-t75, t88);
	t71 = sin(t74);
	t72 = cos(t74);
	t65 = -t71 * t75 + t72 * t88;
	t64 = 0.1e1 / t65 ^ 2;
	t114 = -t106 * t93 + t110 * t113;
	t118 = t100 * t106;
	t79 = t104 * t114 - t108 * t118;
	t125 = t64 * t79;
	t107 = cos(qJ(5));
	t103 = sin(qJ(5));
	t111 = t113 * t102;
	t84 = -t106 * t111 - t110 * t112;
	t120 = t84 * t103;
	t80 = t104 * t118 + t108 * t114;
	t70 = t80 * t107 - t120;
	t68 = 0.1e1 / t70 ^ 2;
	t119 = t84 * t107;
	t69 = t80 * t103 + t119;
	t124 = t68 * t69;
	t123 = t72 * t75;
	t87 = 0.1e1 / t88 ^ 2;
	t122 = t75 * t87;
	t121 = t79 ^ 2 * t64;
	t116 = t69 ^ 2 * t68 + 0.1e1;
	t77 = -t104 * t117 + t82 * t108;
	t115 = -t71 * t88 - t123;
	t91 = t113 * t100;
	t89 = t102 * t104 + t92 * t108;
	t86 = 0.1e1 / t88;
	t81 = -t106 * t112 + t110 * t111;
	t73 = 0.1e1 / (t75 ^ 2 * t87 + 0.1e1);
	t67 = 0.1e1 / t70;
	t66 = 0.1e1 / t116;
	t63 = 0.1e1 / t65;
	t62 = 0.1e1 / (0.1e1 + t121);
	t61 = (t91 * t122 - t81 * t86) * t73 * t104;
	t60 = (t89 * t122 - t77 * t86) * t73;
	t1 = [-t79 * t86 * t73, t61, 0, t60, 0, 0; (-t75 * t63 - (-t71 + (t86 * t123 + t71) * t73) * t121) * t62, (t84 * t104 * t63 - (t115 * t61 + (-t71 * t81 + t72 * t91) * t104) * t125) * t62, 0, (t80 * t63 - (t115 * t60 - t71 * t77 + t72 * t89) * t125) * t62, 0, 0; ((-t103 * t77 - t81 * t107) * t67 - (t81 * t103 - t107 * t77) * t124) * t66, ((-t107 * t114 + t108 * t120) * t67 - (t103 * t114 + t108 * t119) * t124) * t66, 0, (-t103 * t67 + t107 * t124) * t79 * t66, t116 * t66, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (1006->40), mult. (2516->96), div. (90->9), fcn. (3539->15), ass. (0->59)
	t115 = sin(pkin(6));
	t114 = sin(pkin(12));
	t116 = cos(pkin(12));
	t119 = sin(qJ(2));
	t122 = cos(qJ(2));
	t125 = t122 * t114 + t119 * t116;
	t104 = t125 * t115;
	t117 = cos(pkin(6));
	t118 = sin(qJ(4));
	t121 = cos(qJ(4));
	t100 = t104 * t118 - t117 * t121;
	t123 = cos(qJ(1));
	t129 = t115 * t123;
	t105 = t125 * t117;
	t106 = t119 * t114 - t122 * t116;
	t120 = sin(qJ(1));
	t94 = t123 * t105 - t120 * t106;
	t87 = t94 * t118 + t121 * t129;
	t86 = atan2(-t87, t100);
	t83 = sin(t86);
	t84 = cos(t86);
	t77 = t84 * t100 - t83 * t87;
	t76 = 0.1e1 / t77 ^ 2;
	t126 = -t120 * t105 - t123 * t106;
	t130 = t115 * t120;
	t91 = t118 * t126 - t121 * t130;
	t137 = t76 * t91;
	t113 = qJ(5) + qJ(6);
	t112 = cos(t113);
	t111 = sin(t113);
	t124 = t106 * t117;
	t96 = t120 * t124 - t123 * t125;
	t132 = t96 * t111;
	t92 = t118 * t130 + t121 * t126;
	t82 = t92 * t112 - t132;
	t80 = 0.1e1 / t82 ^ 2;
	t131 = t96 * t112;
	t81 = t92 * t111 + t131;
	t136 = t80 * t81;
	t135 = t84 * t87;
	t99 = 0.1e1 / t100 ^ 2;
	t134 = t87 * t99;
	t133 = t91 ^ 2 * t76;
	t128 = t81 ^ 2 * t80 + 0.1e1;
	t89 = -t118 * t129 + t94 * t121;
	t127 = -t100 * t83 - t135;
	t103 = t106 * t115;
	t101 = t104 * t121 + t117 * t118;
	t98 = 0.1e1 / t100;
	t93 = -t120 * t125 - t123 * t124;
	t85 = 0.1e1 / (t87 ^ 2 * t99 + 0.1e1);
	t79 = 0.1e1 / t82;
	t78 = 0.1e1 / t128;
	t75 = 0.1e1 / t77;
	t74 = 0.1e1 / (0.1e1 + t133);
	t73 = (-t103 * t134 - t93 * t98) * t85 * t118;
	t72 = (t101 * t134 - t89 * t98) * t85;
	t71 = t128 * t78;
	t1 = [-t91 * t98 * t85, t73, 0, t72, 0, 0; (-t87 * t75 - (-t83 + (t98 * t135 + t83) * t85) * t133) * t74, (t96 * t118 * t75 - (t127 * t73 + (-t103 * t84 - t83 * t93) * t118) * t137) * t74, 0, (t92 * t75 - (t84 * t101 + t127 * t72 - t83 * t89) * t137) * t74, 0, 0; ((-t111 * t89 - t93 * t112) * t79 - (t93 * t111 - t112 * t89) * t136) * t78, ((-t112 * t126 + t121 * t132) * t79 - (t111 * t126 + t121 * t131) * t136) * t78, 0, (-t111 * t79 + t112 * t136) * t91 * t78, t71, t71;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end