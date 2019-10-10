% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR6
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
%   Wie in S6RRPRPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:33
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (675->32), mult. (1865->82), div. (77->9), fcn. (2642->13), ass. (0->50)
	t85 = cos(pkin(6));
	t82 = sin(pkin(11));
	t84 = cos(pkin(11));
	t87 = sin(qJ(2));
	t90 = cos(qJ(2));
	t93 = t90 * t82 + t87 * t84;
	t76 = t93 * t85;
	t77 = t87 * t82 - t90 * t84;
	t88 = sin(qJ(1));
	t91 = cos(qJ(1));
	t65 = t91 * t76 - t88 * t77;
	t86 = sin(qJ(4));
	t89 = cos(qJ(4));
	t83 = sin(pkin(6));
	t95 = t83 * t91;
	t57 = t65 * t86 + t89 * t95;
	t75 = t93 * t83;
	t71 = t75 * t86 - t85 * t89;
	t56 = atan2(-t57, t71);
	t53 = sin(t56);
	t54 = cos(t56);
	t51 = -t53 * t57 + t54 * t71;
	t50 = 0.1e1 / t51 ^ 2;
	t68 = -t88 * t76 - t91 * t77;
	t96 = t83 * t88;
	t60 = t68 * t86 - t89 * t96;
	t101 = t50 * t60;
	t100 = t54 * t57;
	t70 = 0.1e1 / t71 ^ 2;
	t99 = t57 * t70;
	t98 = t60 ^ 2 * t50;
	t61 = t68 * t89 + t86 * t96;
	t92 = t77 * t85;
	t66 = t88 * t92 - t91 * t93;
	t63 = 0.1e1 / t66 ^ 2;
	t97 = t61 * t63;
	t59 = t65 * t89 - t86 * t95;
	t94 = -t53 * t71 - t100;
	t74 = t77 * t83;
	t72 = t75 * t89 + t85 * t86;
	t69 = 0.1e1 / t71;
	t64 = -t88 * t93 - t91 * t92;
	t62 = 0.1e1 / t66;
	t55 = 0.1e1 / (t57 ^ 2 * t70 + 0.1e1);
	t52 = 0.1e1 / (t61 ^ 2 * t63 + 0.1e1);
	t49 = 0.1e1 / t51;
	t48 = 0.1e1 / (0.1e1 + t98);
	t47 = (-t64 * t69 - t74 * t99) * t86 * t55;
	t46 = (-t59 * t69 + t72 * t99) * t55;
	t1 = [-t60 * t69 * t55, t47, 0, t46, 0, 0; (-t57 * t49 - (-t53 + (t69 * t100 + t53) * t55) * t98) * t48, (t66 * t86 * t49 - ((-t53 * t64 - t54 * t74) * t86 + t94 * t47) * t101) * t48, 0, (t61 * t49 - (t94 * t46 - t53 * t59 + t54 * t72) * t101) * t48, 0, 0; (t59 * t62 - t64 * t97) * t52, (-t66 * t89 * t62 - t68 * t97) * t52, 0, t60 * t62 * t52, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:33
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (867->39), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->57)
	t100 = sin(qJ(4));
	t104 = cos(qJ(4));
	t106 = cos(qJ(1));
	t96 = sin(pkin(6));
	t113 = t106 * t96;
	t102 = sin(qJ(1));
	t101 = sin(qJ(2));
	t105 = cos(qJ(2));
	t95 = sin(pkin(11));
	t97 = cos(pkin(11));
	t109 = t101 * t97 + t105 * t95;
	t98 = cos(pkin(6));
	t90 = t109 * t98;
	t91 = t101 * t95 - t105 * t97;
	t79 = -t102 * t91 + t106 * t90;
	t73 = -t100 * t113 + t104 * t79;
	t89 = t109 * t96;
	t86 = t100 * t98 + t104 * t89;
	t71 = atan2(-t73, t86);
	t68 = sin(t71);
	t69 = cos(t71);
	t62 = -t68 * t73 + t69 * t86;
	t61 = 0.1e1 / t62 ^ 2;
	t108 = -t102 * t90 - t106 * t91;
	t114 = t102 * t96;
	t77 = t100 * t114 + t104 * t108;
	t120 = t61 * t77;
	t103 = cos(qJ(6));
	t107 = t91 * t98;
	t81 = t102 * t107 - t106 * t109;
	t112 = t81 * t103;
	t76 = t100 * t108 - t104 * t114;
	t99 = sin(qJ(6));
	t67 = t76 * t99 - t112;
	t65 = 0.1e1 / t67 ^ 2;
	t115 = t81 * t99;
	t66 = -t103 * t76 - t115;
	t119 = t65 * t66;
	t118 = t69 * t73;
	t84 = 0.1e1 / t86 ^ 2;
	t117 = t73 * t84;
	t116 = t77 ^ 2 * t61;
	t111 = t65 * t66 ^ 2 + 0.1e1;
	t110 = -t68 * t86 - t118;
	t72 = t100 * t79 + t104 * t113;
	t88 = t91 * t96;
	t85 = -t100 * t89 + t104 * t98;
	t83 = 0.1e1 / t86;
	t78 = -t102 * t109 - t106 * t107;
	t70 = 0.1e1 / (t73 ^ 2 * t84 + 0.1e1);
	t64 = 0.1e1 / t67;
	t63 = 0.1e1 / t111;
	t60 = 0.1e1 / t62;
	t59 = 0.1e1 / (0.1e1 + t116);
	t58 = (-t117 * t88 - t78 * t83) * t70 * t104;
	t57 = (t117 * t85 + t72 * t83) * t70;
	t1 = [-t77 * t83 * t70, t58, 0, t57, 0, 0; (-t73 * t60 - (-t68 + (t118 * t83 + t68) * t70) * t116) * t59, (t81 * t104 * t60 - (t110 * t58 + (-t68 * t78 - t69 * t88) * t104) * t120) * t59, 0, (-t76 * t60 - (t110 * t57 + t68 * t72 + t69 * t85) * t120) * t59, 0, 0; ((t103 * t72 + t78 * t99) * t64 - (t103 * t78 - t72 * t99) * t119) * t63, ((-t100 * t112 + t108 * t99) * t64 - (t100 * t115 + t103 * t108) * t119) * t63, 0, (-t103 * t64 - t119 * t99) * t77 * t63, 0, t111 * t63;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end