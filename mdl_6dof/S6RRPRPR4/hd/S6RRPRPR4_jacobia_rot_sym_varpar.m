% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR4
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
%   Wie in S6RRPRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:09
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
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
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:42
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
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (422->27), mult. (1051->69), div. (55->9), fcn. (1490->13), ass. (0->46)
	t83 = sin(qJ(1));
	t85 = cos(qJ(1));
	t79 = sin(pkin(11));
	t81 = cos(pkin(11));
	t84 = cos(qJ(2));
	t92 = cos(pkin(6));
	t89 = t84 * t92;
	t82 = sin(qJ(2));
	t90 = t82 * t92;
	t86 = -t79 * t90 + t81 * t89;
	t87 = t84 * t79 + t82 * t81;
	t62 = -t83 * t87 + t85 * t86;
	t73 = t82 * t79 - t84 * t81;
	t80 = sin(pkin(6));
	t70 = t73 * t80;
	t56 = atan2(t62, t70);
	t54 = cos(t56);
	t97 = t54 * t62;
	t72 = t79 * t89 + t81 * t90;
	t66 = -t83 * t72 - t85 * t73;
	t78 = qJ(4) + pkin(12);
	t76 = sin(t78);
	t77 = cos(t78);
	t94 = t80 * t83;
	t60 = t66 * t77 + t76 * t94;
	t58 = 0.1e1 / t60 ^ 2;
	t59 = t66 * t76 - t77 * t94;
	t96 = t58 * t59;
	t53 = sin(t56);
	t51 = t53 * t62 + t54 * t70;
	t50 = 0.1e1 / t51 ^ 2;
	t64 = -t83 * t86 - t85 * t87;
	t95 = t64 ^ 2 * t50;
	t93 = t80 * t85;
	t91 = t59 ^ 2 * t58 + 0.1e1;
	t88 = -t85 * t72 + t83 * t73;
	t71 = t87 * t80;
	t69 = 0.1e1 / t70 ^ 2;
	t68 = 0.1e1 / t70;
	t57 = 0.1e1 / t60;
	t55 = 0.1e1 / (t62 ^ 2 * t69 + 0.1e1);
	t52 = 0.1e1 / t91;
	t49 = 0.1e1 / t51;
	t48 = 0.1e1 / (0.1e1 + t95);
	t47 = (-t62 * t69 * t71 + t68 * t88) * t55;
	t1 = [t64 * t68 * t55, t47, 0, 0, 0, 0; (t62 * t49 + (t53 + (t68 * t97 - t53) * t55) * t95) * t48, (t66 * t49 + (t53 * t88 + t54 * t71 + (-t53 * t70 + t97) * t47) * t64 * t50) * t48, 0, 0, 0, 0; ((t76 * t88 - t77 * t93) * t57 - (t76 * t93 + t77 * t88) * t96) * t52, (t76 * t57 - t77 * t96) * t64 * t52, 0, t91 * t52, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (1286->40), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->58)
	t107 = qJ(4) + pkin(12);
	t105 = sin(t107);
	t106 = cos(t107);
	t109 = sin(pkin(6));
	t117 = cos(qJ(1));
	t123 = t109 * t117;
	t108 = sin(pkin(11));
	t110 = cos(pkin(11));
	t113 = sin(qJ(2));
	t116 = cos(qJ(2));
	t100 = t113 * t108 - t116 * t110;
	t114 = sin(qJ(1));
	t111 = cos(pkin(6));
	t119 = t116 * t108 + t113 * t110;
	t99 = t119 * t111;
	t88 = -t114 * t100 + t117 * t99;
	t81 = t88 * t105 + t106 * t123;
	t98 = t119 * t109;
	t94 = t98 * t105 - t111 * t106;
	t80 = atan2(-t81, t94);
	t77 = sin(t80);
	t78 = cos(t80);
	t71 = -t77 * t81 + t78 * t94;
	t70 = 0.1e1 / t71 ^ 2;
	t120 = -t117 * t100 - t114 * t99;
	t124 = t109 * t114;
	t85 = t105 * t120 - t106 * t124;
	t131 = t70 * t85;
	t115 = cos(qJ(6));
	t112 = sin(qJ(6));
	t118 = t100 * t111;
	t90 = t114 * t118 - t117 * t119;
	t126 = t90 * t112;
	t86 = t105 * t124 + t106 * t120;
	t76 = t86 * t115 - t126;
	t74 = 0.1e1 / t76 ^ 2;
	t125 = t90 * t115;
	t75 = t86 * t112 + t125;
	t130 = t74 * t75;
	t129 = t78 * t81;
	t93 = 0.1e1 / t94 ^ 2;
	t128 = t81 * t93;
	t127 = t85 ^ 2 * t70;
	t122 = t75 ^ 2 * t74 + 0.1e1;
	t83 = -t105 * t123 + t88 * t106;
	t121 = -t77 * t94 - t129;
	t97 = t100 * t109;
	t95 = t111 * t105 + t98 * t106;
	t92 = 0.1e1 / t94;
	t87 = -t114 * t119 - t117 * t118;
	t79 = 0.1e1 / (t81 ^ 2 * t93 + 0.1e1);
	t73 = 0.1e1 / t76;
	t72 = 0.1e1 / t122;
	t69 = 0.1e1 / t71;
	t68 = 0.1e1 / (0.1e1 + t127);
	t67 = (-t97 * t128 - t87 * t92) * t79 * t105;
	t66 = (t95 * t128 - t83 * t92) * t79;
	t1 = [-t85 * t92 * t79, t67, 0, t66, 0, 0; (-t81 * t69 - (-t77 + (t92 * t129 + t77) * t79) * t127) * t68, (t90 * t105 * t69 - (t121 * t67 + (-t77 * t87 - t78 * t97) * t105) * t131) * t68, 0, (t86 * t69 - (t121 * t66 - t77 * t83 + t78 * t95) * t131) * t68, 0, 0; ((-t112 * t83 - t87 * t115) * t73 - (t87 * t112 - t115 * t83) * t130) * t72, ((t106 * t126 - t115 * t120) * t73 - (t106 * t125 + t112 * t120) * t130) * t72, 0, (-t112 * t73 + t115 * t130) * t85 * t72, 0, t122 * t72;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end