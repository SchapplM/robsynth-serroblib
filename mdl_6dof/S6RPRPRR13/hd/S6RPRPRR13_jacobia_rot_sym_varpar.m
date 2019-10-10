% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR13
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
%   Wie in S6RPRPRR13_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR13_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR13_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (25->13), mult. (89->32), div. (20->9), fcn. (140->9), ass. (0->23)
	t28 = cos(pkin(6));
	t26 = sin(pkin(6));
	t30 = cos(qJ(1));
	t32 = t30 * t26;
	t22 = atan2(t32, t28);
	t19 = sin(t22);
	t20 = cos(t22);
	t14 = t19 * t32 + t20 * t28;
	t29 = sin(qJ(1));
	t37 = 0.1e1 / t14 ^ 2 * t29 ^ 2;
	t23 = t26 ^ 2;
	t21 = 0.1e1 / (0.1e1 + t30 ^ 2 * t23 / t28 ^ 2);
	t36 = t21 / t28;
	t25 = sin(pkin(12));
	t35 = t29 * t25;
	t27 = cos(pkin(12));
	t34 = t29 * t27;
	t33 = t30 * t25;
	t31 = t30 * t27;
	t18 = -t28 * t35 + t31;
	t17 = t28 * t34 + t33;
	t16 = 0.1e1 / t18 ^ 2;
	t1 = [-t29 * t26 * t36, 0, 0, 0, 0, 0; (0.1e1 / t14 * t32 - (-t20 * t23 * t30 * t36 + (t21 - 0.1e1) * t26 * t19) * t26 * t37) / (t23 * t37 + 0.1e1), 0, 0, 0, 0, 0; ((t28 * t31 - t35) / t18 - (-t28 * t33 - t34) * t17 * t16) / (t17 ^ 2 * t16 + 0.1e1), 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (156->21), mult. (450->48), div. (25->9), fcn. (637->13), ass. (0->38)
	t61 = cos(pkin(6));
	t59 = cos(pkin(12));
	t65 = cos(qJ(1));
	t69 = t65 * t59;
	t56 = sin(pkin(12));
	t63 = sin(qJ(1));
	t72 = t63 * t56;
	t51 = -t61 * t69 + t72;
	t57 = sin(pkin(7));
	t60 = cos(pkin(7));
	t58 = sin(pkin(6));
	t73 = t58 * t65;
	t46 = -t51 * t57 + t60 * t73;
	t50 = -t58 * t59 * t57 + t61 * t60;
	t45 = atan2(t46, t50);
	t42 = sin(t45);
	t43 = cos(t45);
	t37 = t42 * t46 + t43 * t50;
	t70 = t65 * t56;
	t71 = t63 * t59;
	t53 = -t61 * t71 - t70;
	t74 = t58 * t63;
	t47 = t53 * t57 - t60 * t74;
	t75 = t47 ^ 2 / t37 ^ 2;
	t54 = -t61 * t72 + t69;
	t62 = sin(qJ(3));
	t64 = cos(qJ(3));
	t66 = t53 * t60 + t57 * t74;
	t41 = t54 * t64 + t66 * t62;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = t54 * t62 - t66 * t64;
	t68 = t40 ^ 2 * t39 + 0.1e1;
	t67 = t51 * t60 + t57 * t73;
	t52 = -t61 * t70 - t71;
	t49 = 0.1e1 / t50;
	t44 = 0.1e1 / (0.1e1 + t46 ^ 2 / t50 ^ 2);
	t38 = 0.1e1 / t68;
	t1 = [t47 * t49 * t44, 0, 0, 0, 0, 0; (t46 / t37 + (t42 + (t43 * t46 * t49 - t42) * t44) * t75) / (0.1e1 + t75), 0, 0, 0, 0, 0; ((t52 * t62 - t67 * t64) / t41 - (t52 * t64 + t67 * t62) * t40 * t39) * t38, 0, t68 * t38, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (423->35), mult. (1287->76), div. (47->9), fcn. (1778->13), ass. (0->50)
	t75 = sin(pkin(7));
	t78 = cos(pkin(7));
	t79 = cos(pkin(6));
	t74 = sin(pkin(12));
	t83 = cos(qJ(1));
	t90 = t83 * t74;
	t77 = cos(pkin(12));
	t81 = sin(qJ(1));
	t91 = t81 * t77;
	t86 = t79 * t91 + t90;
	t76 = sin(pkin(6));
	t96 = t76 * t81;
	t101 = -t75 * t96 + t86 * t78;
	t89 = t83 * t77;
	t92 = t81 * t74;
	t69 = -t79 * t89 + t92;
	t82 = cos(qJ(3));
	t95 = t76 * t83;
	t88 = t75 * t95;
	t93 = t78 * t82;
	t70 = t79 * t90 + t91;
	t80 = sin(qJ(3));
	t98 = t70 * t80;
	t55 = t69 * t93 + t82 * t88 + t98;
	t97 = t75 * t79;
	t62 = -t82 * t97 + (t74 * t80 - t77 * t93) * t76;
	t53 = atan2(-t55, t62);
	t51 = cos(t53);
	t100 = t51 * t55;
	t50 = sin(t53);
	t49 = -t50 * t55 + t51 * t62;
	t48 = 0.1e1 / t49 ^ 2;
	t71 = -t79 * t92 + t89;
	t58 = t101 * t82 + t71 * t80;
	t99 = t58 ^ 2 * t48;
	t94 = t78 * t80;
	t85 = t69 * t94 - t70 * t82 + t80 * t88;
	t66 = t86 * t75 + t78 * t96;
	t65 = 0.1e1 / t66 ^ 2;
	t64 = 0.1e1 / t66;
	t63 = t80 * t97 + (t74 * t82 + t77 * t94) * t76;
	t61 = 0.1e1 / t62 ^ 2;
	t60 = 0.1e1 / t62;
	t59 = -t101 * t80 + t71 * t82;
	t54 = 0.1e1 / (t59 ^ 2 * t65 + 0.1e1);
	t52 = 0.1e1 / (t55 ^ 2 * t61 + 0.1e1);
	t47 = 0.1e1 / t49;
	t46 = 0.1e1 / (0.1e1 + t99);
	t45 = (t55 * t61 * t63 + t60 * t85) * t52;
	t1 = [-t58 * t60 * t52, 0, t45, 0, 0, 0; ((-t98 + (-t69 * t78 - t88) * t82) * t47 - (-t50 + (t60 * t100 + t50) * t52) * t99) * t46, 0, (t59 * t47 - (t50 * t85 + t51 * t63 + (-t50 * t62 - t100) * t45) * t58 * t48) * t46, 0, 0, 0; (t85 * t64 - (-t69 * t75 + t78 * t95) * t59 * t65) * t54, 0, -t58 * t64 * t54, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:45
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (559->36), mult. (1669->85), div. (55->9), fcn. (2293->15), ass. (0->55)
	t90 = sin(pkin(6));
	t99 = cos(qJ(1));
	t110 = t90 * t99;
	t89 = sin(pkin(7));
	t104 = t89 * t110;
	t91 = cos(pkin(12));
	t105 = t99 * t91;
	t88 = sin(pkin(12));
	t96 = sin(qJ(1));
	t108 = t96 * t88;
	t93 = cos(pkin(6));
	t82 = -t105 * t93 + t108;
	t106 = t99 * t88;
	t107 = t96 * t91;
	t83 = t106 * t93 + t107;
	t92 = cos(pkin(7));
	t95 = sin(qJ(3));
	t98 = cos(qJ(3));
	t70 = (t82 * t92 + t104) * t98 + t83 * t95;
	t111 = t90 * t96;
	t84 = -t107 * t93 - t106;
	t101 = t111 * t89 + t84 * t92;
	t85 = -t108 * t93 + t105;
	t74 = -t101 * t98 + t85 * t95;
	t81 = t111 * t92 - t84 * t89;
	t94 = sin(qJ(5));
	t97 = cos(qJ(5));
	t65 = t74 * t94 + t81 * t97;
	t63 = 0.1e1 / t65 ^ 2;
	t64 = -t74 * t97 + t81 * t94;
	t116 = t63 * t64;
	t109 = t92 * t95;
	t100 = t95 * t104 + t109 * t82 - t83 * t98;
	t112 = t89 * t93;
	t79 = t95 * t112 + (t109 * t91 + t88 * t98) * t90;
	t69 = atan2(t100, t79);
	t67 = cos(t69);
	t115 = t67 * t100;
	t66 = sin(t69);
	t60 = t100 * t66 + t67 * t79;
	t59 = 0.1e1 / t60 ^ 2;
	t75 = t101 * t95 + t85 * t98;
	t114 = t75 ^ 2 * t59;
	t103 = t63 * t64 ^ 2 + 0.1e1;
	t80 = t110 * t92 - t82 * t89;
	t78 = t98 * t112 + (t91 * t92 * t98 - t88 * t95) * t90;
	t77 = 0.1e1 / t79 ^ 2;
	t76 = 0.1e1 / t79;
	t68 = 0.1e1 / (t100 ^ 2 * t77 + 0.1e1);
	t62 = 0.1e1 / t65;
	t61 = 0.1e1 / t103;
	t58 = 0.1e1 / t60;
	t57 = 0.1e1 / (0.1e1 + t114);
	t56 = (-t100 * t77 * t78 + t70 * t76) * t68;
	t1 = [-t75 * t76 * t68, 0, t56, 0, 0, 0; (t100 * t58 - (-t66 + (-t115 * t76 + t66) * t68) * t114) * t57, 0, (-t74 * t58 - (t66 * t70 + t67 * t78 + (-t66 * t79 + t115) * t56) * t75 * t59) * t57, 0, 0, 0; ((t70 * t97 + t80 * t94) * t62 - (-t70 * t94 + t80 * t97) * t116) * t61, 0, (-t116 * t94 - t62 * t97) * t75 * t61, 0, t103 * t61, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:45
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (1440->49), mult. (4121->117), div. (85->9), fcn. (5660->17), ass. (0->70)
	t114 = cos(pkin(6));
	t112 = cos(pkin(12));
	t122 = cos(qJ(1));
	t132 = t122 * t112;
	t109 = sin(pkin(12));
	t118 = sin(qJ(1));
	t135 = t118 * t109;
	t106 = -t114 * t132 + t135;
	t110 = sin(pkin(7));
	t113 = cos(pkin(7));
	t111 = sin(pkin(6));
	t137 = t111 * t122;
	t101 = -t106 * t110 + t113 * t137;
	t116 = sin(qJ(5));
	t120 = cos(qJ(5));
	t133 = t122 * t109;
	t134 = t118 * t112;
	t107 = t114 * t133 + t134;
	t117 = sin(qJ(3));
	t121 = cos(qJ(3));
	t128 = t106 * t113 + t110 * t137;
	t123 = t107 * t117 + t128 * t121;
	t86 = -t101 * t120 + t123 * t116;
	t152 = t101 * t116 + t123 * t120;
	t127 = -t114 * t134 - t133;
	t138 = t111 * t118;
	t149 = t110 * t138 + t127 * t113;
	t97 = -t107 * t121 + t128 * t117;
	t105 = -t111 * t112 * t110 + t114 * t113;
	t136 = t112 * t113;
	t139 = t110 * t114;
	t99 = -t121 * t139 + (t109 * t117 - t121 * t136) * t111;
	t93 = t105 * t116 - t99 * t120;
	t82 = atan2(t152, t93);
	t79 = sin(t82);
	t80 = cos(t82);
	t73 = t152 * t79 + t80 * t93;
	t72 = 0.1e1 / t73 ^ 2;
	t103 = -t127 * t110 + t113 * t138;
	t126 = -t114 * t135 + t132;
	t124 = t126 * t117 - t149 * t121;
	t87 = t103 * t116 - t124 * t120;
	t148 = t72 * t87;
	t119 = cos(qJ(6));
	t115 = sin(qJ(6));
	t98 = t149 * t117 + t126 * t121;
	t143 = t98 * t115;
	t88 = t103 * t120 + t124 * t116;
	t78 = t88 * t119 + t143;
	t76 = 0.1e1 / t78 ^ 2;
	t142 = t98 * t119;
	t77 = t88 * t115 - t142;
	t147 = t76 * t77;
	t146 = t80 * t152;
	t92 = 0.1e1 / t93 ^ 2;
	t145 = t152 * t92;
	t144 = t87 ^ 2 * t72;
	t131 = t77 ^ 2 * t76 + 0.1e1;
	t129 = -t79 * t93 + t146;
	t100 = t117 * t139 + (t109 * t121 + t117 * t136) * t111;
	t94 = t105 * t120 + t99 * t116;
	t91 = 0.1e1 / t93;
	t81 = 0.1e1 / (t152 ^ 2 * t92 + 0.1e1);
	t75 = 0.1e1 / t78;
	t74 = 0.1e1 / t131;
	t71 = 0.1e1 / t73;
	t70 = 0.1e1 / (0.1e1 + t144);
	t69 = (t100 * t145 - t91 * t97) * t81 * t120;
	t68 = (-t94 * t145 - t86 * t91) * t81;
	t1 = [-t87 * t91 * t81, 0, t69, 0, t68, 0; (t152 * t71 - (-t79 + (-t91 * t146 + t79) * t81) * t144) * t70, 0, (-t98 * t120 * t71 - (t129 * t69 + (-t100 * t80 - t79 * t97) * t120) * t148) * t70, 0, (t88 * t71 - (t129 * t68 - t79 * t86 + t80 * t94) * t148) * t70, 0; ((-t115 * t86 - t97 * t119) * t75 - (t97 * t115 - t119 * t86) * t147) * t74, 0, ((t116 * t143 + t124 * t119) * t75 - (-t124 * t115 + t116 * t142) * t147) * t74, 0, (-t115 * t75 + t119 * t147) * t87 * t74, t131 * t74;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end