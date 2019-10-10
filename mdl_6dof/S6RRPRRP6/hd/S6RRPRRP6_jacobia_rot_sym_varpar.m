% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP6
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
%   Wie in S6RRPRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:37
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:18
	% EndTime: 2019-10-10 10:37:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:18
	% EndTime: 2019-10-10 10:37:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:37:18
	% EndTime: 2019-10-10 10:37:19
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
	% StartTime: 2019-10-10 10:37:18
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (374->26), mult. (1051->69), div. (55->9), fcn. (1490->13), ass. (0->45)
	t73 = sin(qJ(1));
	t76 = cos(qJ(1));
	t68 = sin(pkin(11));
	t70 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (867->39), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->57)
	t104 = sin(qJ(4));
	t108 = cos(qJ(4));
	t100 = sin(pkin(6));
	t110 = cos(qJ(1));
	t117 = t100 * t110;
	t106 = sin(qJ(1));
	t101 = cos(pkin(11));
	t105 = sin(qJ(2));
	t109 = cos(qJ(2));
	t99 = sin(pkin(11));
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
	% StartTime: 2019-10-10 10:37:19
	% EndTime: 2019-10-10 10:37:19
	% DurationCPUTime: 0.65s
	% Computational Cost: add. (1741->46), mult. (4726->114), div. (107->9), fcn. (6600->15), ass. (0->65)
	t126 = cos(pkin(6));
	t123 = sin(pkin(11));
	t125 = cos(pkin(11));
	t129 = sin(qJ(2));
	t133 = cos(qJ(2));
	t137 = t133 * t123 + t129 * t125;
	t118 = t137 * t126;
	t121 = t129 * t123 - t133 * t125;
	t130 = sin(qJ(1));
	t134 = cos(qJ(1));
	t109 = t134 * t118 - t130 * t121;
	t128 = sin(qJ(4));
	t132 = cos(qJ(4));
	t124 = sin(pkin(6));
	t143 = t124 * t134;
	t103 = -t109 * t132 + t128 * t143;
	t127 = sin(qJ(5));
	t131 = cos(qJ(5));
	t135 = t121 * t126;
	t141 = -t130 * t137 - t134 * t135;
	t154 = t103 * t127 - t141 * t131;
	t140 = t141 * t127;
	t153 = t103 * t131 + t140;
	t117 = t137 * t124;
	t114 = t117 * t132 + t126 * t128;
	t116 = t121 * t124;
	t98 = t114 * t127 - t116 * t131;
	t86 = atan2(t154, t98);
	t83 = sin(t86);
	t84 = cos(t86);
	t82 = t154 * t83 + t84 * t98;
	t81 = 0.1e1 / t82 ^ 2;
	t138 = -t130 * t118 - t134 * t121;
	t144 = t124 * t130;
	t105 = t128 * t144 + t132 * t138;
	t111 = t130 * t135 - t134 * t137;
	t145 = t111 * t131;
	t93 = t105 * t127 + t145;
	t151 = t81 * t93;
	t150 = t84 * t154;
	t97 = 0.1e1 / t98 ^ 2;
	t149 = t154 * t97;
	t148 = t93 ^ 2 * t81;
	t104 = -t128 * t138 + t132 * t144;
	t94 = t105 * t131 - t111 * t127;
	t89 = 0.1e1 / t94 ^ 2;
	t147 = t104 ^ 2 * t89;
	t146 = t104 * t89;
	t142 = t127 * t132;
	t139 = -t83 * t98 + t150;
	t136 = t109 * t128 + t132 * t143;
	t113 = -t117 * t128 + t126 * t132;
	t106 = -t116 * t142 - t117 * t131;
	t99 = t114 * t131 + t116 * t127;
	t96 = 0.1e1 / t98;
	t95 = -t109 * t131 + t132 * t140;
	t88 = 0.1e1 / t94;
	t87 = 0.1e1 / (0.1e1 + t147);
	t85 = 0.1e1 / (t154 ^ 2 * t97 + 0.1e1);
	t80 = 0.1e1 / t82;
	t79 = 0.1e1 / (0.1e1 + t148);
	t78 = (-t113 * t149 + t136 * t96) * t85 * t127;
	t77 = (-t106 * t149 - t95 * t96) * t85;
	t76 = (-t99 * t149 + t153 * t96) * t85;
	t1 = [-t93 * t96 * t85, t77, 0, t78, t76, 0; (t154 * t80 - (-t83 + (-t96 * t150 + t83) * t85) * t148) * t79, ((t111 * t142 - t131 * t138) * t80 - (t84 * t106 + t139 * t77 - t83 * t95) * t151) * t79, 0, (t104 * t127 * t80 - (t139 * t78 + (t113 * t84 + t136 * t83) * t127) * t151) * t79, (t94 * t80 - (t139 * t76 + t153 * t83 + t84 * t99) * t151) * t79, 0; (t136 * t88 - t153 * t146) * t87, (-t111 * t128 * t88 - (t127 * t138 + t132 * t145) * t146) * t87, 0, (-t105 * t88 - t131 * t147) * t87, t93 * t87 * t146, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end