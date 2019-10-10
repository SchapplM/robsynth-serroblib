% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR4
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
%   Wie in S6RRPRRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
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
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
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
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.22s
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
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (467->27), mult. (1134->69), div. (60->9), fcn. (1602->13), ass. (0->47)
	t89 = sin(qJ(1));
	t91 = cos(qJ(1));
	t85 = sin(pkin(12));
	t87 = cos(pkin(12));
	t90 = cos(qJ(2));
	t98 = cos(pkin(6));
	t95 = t90 * t98;
	t88 = sin(qJ(2));
	t96 = t88 * t98;
	t92 = -t85 * t96 + t87 * t95;
	t93 = t90 * t85 + t88 * t87;
	t68 = -t89 * t93 + t91 * t92;
	t79 = t88 * t85 - t90 * t87;
	t86 = sin(pkin(6));
	t76 = t79 * t86;
	t62 = atan2(t68, t76);
	t60 = cos(t62);
	t103 = t60 * t68;
	t100 = t86 * t89;
	t78 = t85 * t95 + t87 * t96;
	t72 = -t89 * t78 - t91 * t79;
	t84 = qJ(4) + qJ(5);
	t82 = sin(t84);
	t83 = cos(t84);
	t66 = t82 * t100 + t72 * t83;
	t64 = 0.1e1 / t66 ^ 2;
	t65 = -t83 * t100 + t72 * t82;
	t102 = t64 * t65;
	t59 = sin(t62);
	t57 = t59 * t68 + t60 * t76;
	t56 = 0.1e1 / t57 ^ 2;
	t70 = -t89 * t92 - t91 * t93;
	t101 = t70 ^ 2 * t56;
	t99 = t86 * t91;
	t97 = t65 ^ 2 * t64 + 0.1e1;
	t94 = -t91 * t78 + t89 * t79;
	t77 = t93 * t86;
	t75 = 0.1e1 / t76 ^ 2;
	t74 = 0.1e1 / t76;
	t63 = 0.1e1 / t66;
	t61 = 0.1e1 / (t68 ^ 2 * t75 + 0.1e1);
	t58 = 0.1e1 / t97;
	t55 = 0.1e1 / t57;
	t54 = 0.1e1 / (0.1e1 + t101);
	t53 = (-t68 * t75 * t77 + t74 * t94) * t61;
	t52 = t97 * t58;
	t1 = [t70 * t74 * t61, t53, 0, 0, 0, 0; (t68 * t55 + (t59 + (t74 * t103 - t59) * t61) * t101) * t54, (t72 * t55 + (t59 * t94 + t60 * t77 + (-t59 * t76 + t103) * t53) * t70 * t56) * t54, 0, 0, 0, 0; ((t82 * t94 - t83 * t99) * t63 - (t82 * t99 + t83 * t94) * t102) * t58, (-t83 * t102 + t82 * t63) * t70 * t58, 0, t52, t52, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:49
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (1726->40), mult. (3150->96), div. (115->9), fcn. (4438->15), ass. (0->60)
	t123 = qJ(4) + qJ(5);
	t121 = sin(t123);
	t122 = cos(t123);
	t127 = cos(pkin(6));
	t124 = sin(pkin(12));
	t126 = cos(pkin(12));
	t129 = sin(qJ(2));
	t132 = cos(qJ(2));
	t135 = t132 * t124 + t129 * t126;
	t115 = t135 * t127;
	t116 = t129 * t124 - t132 * t126;
	t130 = sin(qJ(1));
	t133 = cos(qJ(1));
	t136 = -t130 * t115 - t133 * t116;
	t125 = sin(pkin(6));
	t140 = t125 * t130;
	t102 = t121 * t140 + t122 * t136;
	t131 = cos(qJ(6));
	t134 = t116 * t127;
	t106 = t130 * t134 - t133 * t135;
	t128 = sin(qJ(6));
	t142 = t106 * t128;
	t92 = t102 * t131 - t142;
	t90 = 0.1e1 / t92 ^ 2;
	t141 = t106 * t131;
	t91 = t102 * t128 + t141;
	t147 = t90 * t91;
	t114 = t135 * t125;
	t110 = t114 * t121 - t127 * t122;
	t104 = t133 * t115 - t130 * t116;
	t139 = t125 * t133;
	t97 = t104 * t121 + t122 * t139;
	t96 = atan2(-t97, t110);
	t94 = cos(t96);
	t146 = t94 * t97;
	t101 = t121 * t136 - t122 * t140;
	t93 = sin(t96);
	t87 = t94 * t110 - t93 * t97;
	t86 = 0.1e1 / t87 ^ 2;
	t145 = t101 * t86;
	t144 = t101 ^ 2 * t86;
	t109 = 0.1e1 / t110 ^ 2;
	t143 = t109 * t97;
	t138 = t91 ^ 2 * t90 + 0.1e1;
	t99 = t104 * t122 - t121 * t139;
	t137 = -t110 * t93 - t146;
	t113 = t116 * t125;
	t111 = t114 * t122 + t127 * t121;
	t108 = 0.1e1 / t110;
	t103 = -t130 * t135 - t133 * t134;
	t95 = 0.1e1 / (t97 ^ 2 * t109 + 0.1e1);
	t89 = 0.1e1 / t92;
	t88 = 0.1e1 / t138;
	t85 = 0.1e1 / t87;
	t84 = 0.1e1 / (0.1e1 + t144);
	t83 = (-t103 * t108 - t113 * t143) * t95 * t121;
	t82 = (-t108 * t99 + t111 * t143) * t95;
	t81 = (-t128 * t89 + t131 * t147) * t88 * t101;
	t80 = (t102 * t85 - (t94 * t111 + t137 * t82 - t93 * t99) * t145) * t84;
	t1 = [-t101 * t108 * t95, t83, 0, t82, t82, 0; (-t97 * t85 - (-t93 + (t108 * t146 + t93) * t95) * t144) * t84, (t106 * t121 * t85 - (t137 * t83 + (-t103 * t93 - t113 * t94) * t121) * t145) * t84, 0, t80, t80, 0; ((-t103 * t131 - t128 * t99) * t89 - (t103 * t128 - t131 * t99) * t147) * t88, ((t122 * t142 - t131 * t136) * t89 - (t122 * t141 + t128 * t136) * t147) * t88, 0, t81, t81, t138 * t88;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end