% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR3
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
%   Wie in S6RRPPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:39
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
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
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:19
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
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (343->26), mult. (968->68), div. (50->9), fcn. (1378->13), ass. (0->44)
	t73 = sin(qJ(1));
	t75 = cos(qJ(1));
	t68 = sin(pkin(11));
	t71 = cos(pkin(11));
	t74 = cos(qJ(2));
	t81 = cos(pkin(6));
	t79 = t74 * t81;
	t72 = sin(qJ(2));
	t80 = t72 * t81;
	t76 = -t68 * t80 + t71 * t79;
	t77 = t74 * t68 + t72 * t71;
	t53 = -t73 * t77 + t75 * t76;
	t64 = t72 * t68 - t74 * t71;
	t69 = sin(pkin(6));
	t61 = t64 * t69;
	t47 = atan2(t53, t61);
	t45 = cos(t47);
	t86 = t45 * t53;
	t63 = t68 * t79 + t71 * t80;
	t57 = -t73 * t63 - t75 * t64;
	t67 = sin(pkin(12));
	t70 = cos(pkin(12));
	t83 = t69 * t73;
	t51 = t57 * t70 + t67 * t83;
	t49 = 0.1e1 / t51 ^ 2;
	t50 = t57 * t67 - t70 * t83;
	t85 = t49 * t50;
	t44 = sin(t47);
	t42 = t44 * t53 + t45 * t61;
	t41 = 0.1e1 / t42 ^ 2;
	t55 = -t73 * t76 - t75 * t77;
	t84 = t55 ^ 2 * t41;
	t82 = t69 * t75;
	t78 = -t75 * t63 + t73 * t64;
	t62 = t77 * t69;
	t60 = 0.1e1 / t61 ^ 2;
	t59 = 0.1e1 / t61;
	t48 = 0.1e1 / t51;
	t46 = 0.1e1 / (t53 ^ 2 * t60 + 0.1e1);
	t43 = 0.1e1 / (t50 ^ 2 * t49 + 0.1e1);
	t40 = 0.1e1 / t42;
	t39 = 0.1e1 / (0.1e1 + t84);
	t38 = (-t53 * t60 * t62 + t59 * t78) * t46;
	t1 = [t55 * t59 * t46, t38, 0, 0, 0, 0; (t53 * t40 + (t44 + (t59 * t86 - t44) * t46) * t84) * t39, (t57 * t40 + (t44 * t78 + t45 * t62 + (-t44 * t61 + t86) * t38) * t55 * t41) * t39, 0, 0, 0, 0; ((t67 * t78 - t70 * t82) * t48 - (t67 * t82 + t70 * t78) * t85) * t43, (t67 * t48 - t70 * t85) * t55 * t43, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (422->27), mult. (1051->69), div. (55->9), fcn. (1490->13), ass. (0->46)
	t82 = sin(qJ(1));
	t84 = cos(qJ(1));
	t78 = sin(pkin(11));
	t80 = cos(pkin(11));
	t83 = cos(qJ(2));
	t91 = cos(pkin(6));
	t88 = t83 * t91;
	t81 = sin(qJ(2));
	t89 = t81 * t91;
	t85 = -t78 * t89 + t80 * t88;
	t86 = t78 * t83 + t80 * t81;
	t61 = -t82 * t86 + t84 * t85;
	t72 = t78 * t81 - t83 * t80;
	t79 = sin(pkin(6));
	t69 = t72 * t79;
	t55 = atan2(t61, t69);
	t52 = sin(t55);
	t53 = cos(t55);
	t50 = t52 * t61 + t53 * t69;
	t49 = 0.1e1 / t50 ^ 2;
	t63 = -t82 * t85 - t84 * t86;
	t96 = t49 * t63 ^ 2;
	t95 = t53 * t61;
	t71 = t78 * t88 + t80 * t89;
	t65 = -t82 * t71 - t72 * t84;
	t77 = pkin(12) + qJ(5);
	t75 = sin(t77);
	t76 = cos(t77);
	t93 = t79 * t82;
	t59 = t65 * t76 + t75 * t93;
	t57 = 0.1e1 / t59 ^ 2;
	t58 = t65 * t75 - t76 * t93;
	t94 = t57 * t58;
	t92 = t79 * t84;
	t90 = t57 * t58 ^ 2 + 0.1e1;
	t87 = -t71 * t84 + t82 * t72;
	t70 = t86 * t79;
	t68 = 0.1e1 / t69 ^ 2;
	t67 = 0.1e1 / t69;
	t56 = 0.1e1 / t59;
	t54 = 0.1e1 / (t61 ^ 2 * t68 + 0.1e1);
	t51 = 0.1e1 / t90;
	t48 = 0.1e1 / t50;
	t47 = 0.1e1 / (0.1e1 + t96);
	t46 = (-t61 * t68 * t70 + t67 * t87) * t54;
	t1 = [t63 * t67 * t54, t46, 0, 0, 0, 0; (t61 * t48 + (t52 + (t67 * t95 - t52) * t54) * t96) * t47, (t65 * t48 + (t52 * t87 + t53 * t70 + (-t52 * t69 + t95) * t46) * t63 * t49) * t47, 0, 0, 0, 0; ((t75 * t87 - t76 * t92) * t56 - (t75 * t92 + t76 * t87) * t94) * t51, (t75 * t56 - t76 * t94) * t63 * t51, 0, 0, t90 * t51, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (1286->40), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->58)
	t106 = pkin(12) + qJ(5);
	t104 = sin(t106);
	t105 = cos(t106);
	t108 = sin(pkin(6));
	t116 = cos(qJ(1));
	t122 = t108 * t116;
	t113 = sin(qJ(1));
	t110 = cos(pkin(6));
	t107 = sin(pkin(11));
	t109 = cos(pkin(11));
	t112 = sin(qJ(2));
	t115 = cos(qJ(2));
	t118 = t115 * t107 + t112 * t109;
	t98 = t118 * t110;
	t99 = t112 * t107 - t115 * t109;
	t87 = -t113 * t99 + t116 * t98;
	t80 = t87 * t104 + t105 * t122;
	t97 = t118 * t108;
	t93 = t97 * t104 - t110 * t105;
	t79 = atan2(-t80, t93);
	t76 = sin(t79);
	t77 = cos(t79);
	t70 = -t76 * t80 + t77 * t93;
	t69 = 0.1e1 / t70 ^ 2;
	t119 = -t113 * t98 - t116 * t99;
	t123 = t108 * t113;
	t84 = t104 * t119 - t105 * t123;
	t130 = t69 * t84;
	t114 = cos(qJ(6));
	t111 = sin(qJ(6));
	t117 = t99 * t110;
	t89 = t113 * t117 - t116 * t118;
	t125 = t89 * t111;
	t85 = t104 * t123 + t105 * t119;
	t75 = t85 * t114 - t125;
	t73 = 0.1e1 / t75 ^ 2;
	t124 = t89 * t114;
	t74 = t85 * t111 + t124;
	t129 = t73 * t74;
	t128 = t77 * t80;
	t92 = 0.1e1 / t93 ^ 2;
	t127 = t80 * t92;
	t126 = t84 ^ 2 * t69;
	t121 = t74 ^ 2 * t73 + 0.1e1;
	t82 = -t104 * t122 + t87 * t105;
	t120 = -t76 * t93 - t128;
	t96 = t99 * t108;
	t94 = t110 * t104 + t97 * t105;
	t91 = 0.1e1 / t93;
	t86 = -t113 * t118 - t116 * t117;
	t78 = 0.1e1 / (t80 ^ 2 * t92 + 0.1e1);
	t72 = 0.1e1 / t75;
	t71 = 0.1e1 / t121;
	t68 = 0.1e1 / t70;
	t67 = 0.1e1 / (0.1e1 + t126);
	t66 = (-t96 * t127 - t86 * t91) * t78 * t104;
	t65 = (t94 * t127 - t82 * t91) * t78;
	t1 = [-t84 * t91 * t78, t66, 0, 0, t65, 0; (-t80 * t68 - (-t76 + (t91 * t128 + t76) * t78) * t126) * t67, (t89 * t104 * t68 - (t120 * t66 + (-t76 * t86 - t77 * t96) * t104) * t130) * t67, 0, 0, (t85 * t68 - (t120 * t65 - t76 * t82 + t77 * t94) * t130) * t67, 0; ((-t111 * t82 - t86 * t114) * t72 - (t86 * t111 - t114 * t82) * t129) * t71, ((t105 * t125 - t114 * t119) * t72 - (t105 * t124 + t111 * t119) * t129) * t71, 0, 0, (-t111 * t72 + t114 * t129) * t84 * t71, t121 * t71;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end