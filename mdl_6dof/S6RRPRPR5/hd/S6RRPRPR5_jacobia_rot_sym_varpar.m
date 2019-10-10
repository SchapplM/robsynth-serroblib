% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR5
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
%   Wie in S6RRPRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (808->39), mult. (2210->96), div. (80->9), fcn. (3119->15), ass. (0->55)
	t100 = cos(qJ(4));
	t102 = cos(qJ(1));
	t93 = sin(pkin(6));
	t108 = t102 * t93;
	t101 = cos(qJ(2));
	t92 = sin(pkin(11));
	t95 = cos(pkin(11));
	t98 = sin(qJ(2));
	t106 = t101 * t95 - t98 * t92;
	t105 = t101 * t92 + t98 * t95;
	t96 = cos(pkin(6));
	t85 = t105 * t96;
	t99 = sin(qJ(1));
	t74 = t102 * t85 + t106 * t99;
	t97 = sin(qJ(4));
	t67 = t100 * t108 + t74 * t97;
	t84 = t105 * t93;
	t80 = -t96 * t100 + t84 * t97;
	t66 = atan2(-t67, t80);
	t63 = sin(t66);
	t64 = cos(t66);
	t57 = -t63 * t67 + t64 * t80;
	t56 = 0.1e1 / t57 ^ 2;
	t104 = t102 * t106 - t99 * t85;
	t110 = t93 * t99;
	t71 = -t100 * t110 + t104 * t97;
	t115 = t56 * t71;
	t72 = t100 * t104 + t97 * t110;
	t103 = t106 * t96;
	t76 = -t102 * t105 - t99 * t103;
	t91 = sin(pkin(12));
	t94 = cos(pkin(12));
	t62 = t72 * t94 - t76 * t91;
	t60 = 0.1e1 / t62 ^ 2;
	t61 = t72 * t91 + t76 * t94;
	t114 = t60 * t61;
	t113 = t64 * t67;
	t79 = 0.1e1 / t80 ^ 2;
	t112 = t67 * t79;
	t111 = t71 ^ 2 * t56;
	t109 = t100 * t76;
	t69 = t74 * t100 - t97 * t108;
	t107 = -t63 * t80 - t113;
	t83 = t106 * t93;
	t81 = t84 * t100 + t96 * t97;
	t78 = 0.1e1 / t80;
	t73 = t102 * t103 - t105 * t99;
	t65 = 0.1e1 / (t67 ^ 2 * t79 + 0.1e1);
	t59 = 0.1e1 / t62;
	t58 = 0.1e1 / (t61 ^ 2 * t60 + 0.1e1);
	t55 = 0.1e1 / t57;
	t54 = 0.1e1 / (0.1e1 + t111);
	t53 = (t83 * t112 - t73 * t78) * t97 * t65;
	t52 = (t81 * t112 - t69 * t78) * t65;
	t1 = [-t71 * t78 * t65, t53, 0, t52, 0, 0; (-t67 * t55 - (-t63 + (t78 * t113 + t63) * t65) * t111) * t54, (t76 * t97 * t55 - ((-t63 * t73 + t64 * t83) * t97 + t107 * t53) * t115) * t54, 0, (t72 * t55 - (t107 * t52 - t63 * t69 + t64 * t81) * t115) * t54, 0, 0; ((-t69 * t91 - t73 * t94) * t59 - (-t69 * t94 + t73 * t91) * t114) * t58, ((-t104 * t94 + t91 * t109) * t59 - (t104 * t91 + t94 * t109) * t114) * t58, 0, (t94 * t114 - t91 * t59) * t71 * t58, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (933->40), mult. (2363->96), div. (85->9), fcn. (3329->15), ass. (0->58)
	t107 = sin(qJ(4));
	t110 = cos(qJ(4));
	t104 = sin(pkin(6));
	t112 = cos(qJ(1));
	t118 = t104 * t112;
	t109 = sin(qJ(1));
	t106 = cos(pkin(6));
	t103 = sin(pkin(11));
	t105 = cos(pkin(11));
	t108 = sin(qJ(2));
	t111 = cos(qJ(2));
	t114 = t111 * t103 + t108 * t105;
	t94 = t114 * t106;
	t95 = t108 * t103 - t111 * t105;
	t83 = -t109 * t95 + t112 * t94;
	t76 = t83 * t107 + t110 * t118;
	t93 = t114 * t104;
	t89 = -t106 * t110 + t93 * t107;
	t75 = atan2(-t76, t89);
	t72 = sin(t75);
	t73 = cos(t75);
	t66 = -t72 * t76 + t73 * t89;
	t65 = 0.1e1 / t66 ^ 2;
	t115 = -t109 * t94 - t112 * t95;
	t119 = t104 * t109;
	t80 = t107 * t115 - t110 * t119;
	t126 = t65 * t80;
	t102 = pkin(12) + qJ(6);
	t101 = cos(t102);
	t100 = sin(t102);
	t113 = t95 * t106;
	t85 = t109 * t113 - t112 * t114;
	t121 = t85 * t100;
	t81 = t107 * t119 + t110 * t115;
	t71 = t81 * t101 - t121;
	t69 = 0.1e1 / t71 ^ 2;
	t120 = t85 * t101;
	t70 = t81 * t100 + t120;
	t125 = t69 * t70;
	t124 = t73 * t76;
	t88 = 0.1e1 / t89 ^ 2;
	t123 = t76 * t88;
	t122 = t80 ^ 2 * t65;
	t117 = t70 ^ 2 * t69 + 0.1e1;
	t78 = -t107 * t118 + t83 * t110;
	t116 = -t72 * t89 - t124;
	t92 = t95 * t104;
	t90 = t106 * t107 + t93 * t110;
	t87 = 0.1e1 / t89;
	t82 = -t109 * t114 - t112 * t113;
	t74 = 0.1e1 / (t76 ^ 2 * t88 + 0.1e1);
	t68 = 0.1e1 / t71;
	t67 = 0.1e1 / t117;
	t64 = 0.1e1 / t66;
	t63 = 0.1e1 / (0.1e1 + t122);
	t62 = (-t92 * t123 - t82 * t87) * t74 * t107;
	t61 = (t90 * t123 - t78 * t87) * t74;
	t1 = [-t80 * t87 * t74, t62, 0, t61, 0, 0; (-t76 * t64 - (-t72 + (t87 * t124 + t72) * t74) * t122) * t63, (t85 * t107 * t64 - (t116 * t62 + (-t72 * t82 - t73 * t92) * t107) * t126) * t63, 0, (t81 * t64 - (t116 * t61 - t72 * t78 + t73 * t90) * t126) * t63, 0, 0; ((-t100 * t78 - t82 * t101) * t68 - (t82 * t100 - t101 * t78) * t125) * t67, ((-t101 * t115 + t110 * t121) * t68 - (t100 * t115 + t110 * t120) * t125) * t67, 0, (-t100 * t68 + t101 * t125) * t80 * t67, 0, t117 * t67;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end