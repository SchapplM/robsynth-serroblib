% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR11
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
%   Wie in S6RRRPRR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR11_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR11_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
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
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
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
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (359->30), mult. (1033->74), div. (77->9), fcn. (1491->11), ass. (0->49)
	t67 = cos(pkin(6));
	t69 = sin(qJ(2));
	t73 = cos(qJ(1));
	t76 = t73 * t69;
	t70 = sin(qJ(1));
	t72 = cos(qJ(2));
	t77 = t70 * t72;
	t60 = t67 * t76 + t77;
	t68 = sin(qJ(3));
	t71 = cos(qJ(3));
	t66 = sin(pkin(6));
	t79 = t66 * t73;
	t49 = t60 * t68 + t71 * t79;
	t82 = t66 * t68;
	t57 = -t67 * t71 + t69 * t82;
	t46 = atan2(-t49, t57);
	t42 = sin(t46);
	t43 = cos(t46);
	t41 = -t42 * t49 + t43 * t57;
	t40 = 0.1e1 / t41 ^ 2;
	t75 = t73 * t72;
	t78 = t70 * t69;
	t62 = -t67 * t78 + t75;
	t81 = t66 * t71;
	t52 = t62 * t68 - t70 * t81;
	t88 = t40 * t52;
	t87 = t43 * t49;
	t53 = t62 * t71 + t70 * t82;
	t48 = 0.1e1 / t53 ^ 2;
	t61 = -t67 * t77 - t76;
	t86 = t48 * t61;
	t55 = 0.1e1 / t57 ^ 2;
	t85 = t49 * t55;
	t84 = t52 ^ 2 * t40;
	t83 = t61 ^ 2 * t48;
	t80 = t66 * t72;
	t51 = t60 * t71 - t68 * t79;
	t74 = -t42 * t57 - t87;
	t59 = -t67 * t75 + t78;
	t58 = t67 * t68 + t69 * t81;
	t54 = 0.1e1 / t57;
	t47 = 0.1e1 / t53;
	t45 = 0.1e1 / (0.1e1 + t83);
	t44 = 0.1e1 / (t49 ^ 2 * t55 + 0.1e1);
	t39 = 0.1e1 / t41;
	t38 = 0.1e1 / (0.1e1 + t84);
	t37 = (t54 * t59 + t80 * t85) * t68 * t44;
	t36 = (-t51 * t54 + t58 * t85) * t44;
	t1 = [-t52 * t54 * t44, t37, t36, 0, 0, 0; (-t49 * t39 - (-t42 + (t54 * t87 + t42) * t44) * t84) * t38, (t61 * t68 * t39 - ((t42 * t59 + t43 * t80) * t68 + t74 * t37) * t88) * t38, (t53 * t39 - (t74 * t36 - t42 * t51 + t43 * t58) * t88) * t38, 0, 0, 0; (t59 * t47 + t51 * t86) * t45, (-t47 * t62 - t71 * t83) * t45, t52 * t45 * t86, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (270->29), mult. (779->75), div. (78->11), fcn. (1122->13), ass. (0->49)
	t70 = cos(pkin(6));
	t77 = cos(qJ(2));
	t78 = cos(qJ(1));
	t80 = t78 * t77;
	t73 = sin(qJ(2));
	t74 = sin(qJ(1));
	t83 = t74 * t73;
	t63 = -t70 * t83 + t80;
	t72 = sin(qJ(3));
	t76 = cos(qJ(3));
	t69 = sin(pkin(6));
	t86 = t69 * t74;
	t53 = t63 * t72 - t76 * t86;
	t54 = t63 * t76 + t72 * t86;
	t71 = sin(qJ(5));
	t75 = cos(qJ(5));
	t50 = t53 * t71 + t54 * t75;
	t48 = 0.1e1 / t50 ^ 2;
	t49 = -t53 * t75 + t54 * t71;
	t91 = t48 * t49;
	t90 = t49 ^ 2 * t48;
	t59 = -t70 * t80 + t83;
	t85 = t69 * t77;
	t58 = atan2(t59, t85);
	t56 = cos(t58);
	t89 = t56 * t59;
	t55 = sin(t58);
	t46 = t55 * t59 + t56 * t85;
	t45 = 0.1e1 / t46 ^ 2;
	t81 = t78 * t73;
	t82 = t74 * t77;
	t61 = t70 * t82 + t81;
	t88 = t61 ^ 2 * t45;
	t66 = 0.1e1 / t69;
	t67 = 0.1e1 / t77;
	t87 = t66 * t67;
	t84 = t69 * t78;
	t79 = 0.1e1 + t90;
	t68 = 0.1e1 / t77 ^ 2;
	t60 = t70 * t81 + t82;
	t57 = 0.1e1 / (0.1e1 + t59 ^ 2 / t69 ^ 2 * t68);
	t52 = -t60 * t76 + t72 * t84;
	t51 = -t60 * t72 - t76 * t84;
	t47 = 0.1e1 / t50;
	t44 = 0.1e1 / t46;
	t43 = 0.1e1 / (0.1e1 + t88);
	t42 = (t59 * t68 * t73 + t60 * t67) * t66 * t57;
	t41 = 0.1e1 / t79;
	t1 = [t61 * t57 * t87, t42, 0, 0, 0, 0; (t59 * t44 + (t55 + (t87 * t89 - t55) * t57) * t88) * t43, (-t63 * t44 + (-t56 * t69 * t73 + t55 * t60 + (-t55 * t85 + t89) * t42) * t61 * t45) * t43, 0, 0, 0, 0; ((-t51 * t75 + t52 * t71) * t47 - (t51 * t71 + t52 * t75) * t91) * t41, ((-t71 * t76 + t72 * t75) * t47 - (-t71 * t72 - t75 * t76) * t91) * t41 * t61, (-t47 * t50 - t90) * t41, 0, t79 * t41, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (1357->47), mult. (3691->106), div. (115->9), fcn. (5181->15), ass. (0->66)
	t130 = cos(pkin(6));
	t134 = sin(qJ(2));
	t140 = cos(qJ(1));
	t148 = t140 * t134;
	t135 = sin(qJ(1));
	t139 = cos(qJ(2));
	t149 = t135 * t139;
	t123 = t130 * t148 + t149;
	t133 = sin(qJ(3));
	t138 = cos(qJ(3));
	t129 = sin(pkin(6));
	t151 = t129 * t140;
	t117 = -t123 * t138 + t133 * t151;
	t132 = sin(qJ(5));
	t137 = cos(qJ(5));
	t141 = t123 * t133 + t138 * t151;
	t102 = t117 * t137 - t141 * t132;
	t147 = t140 * t139;
	t150 = t135 * t134;
	t125 = -t130 * t150 + t147;
	t153 = t129 * t133;
	t118 = t125 * t138 + t135 * t153;
	t152 = t129 * t138;
	t145 = t125 * t133 - t135 * t152;
	t105 = t118 * t137 + t145 * t132;
	t120 = -t130 * t138 + t134 * t153;
	t121 = t130 * t133 + t134 * t152;
	t112 = t120 * t132 + t121 * t137;
	t111 = -t120 * t137 + t121 * t132;
	t164 = t117 * t132 + t141 * t137;
	t93 = atan2(t164, t111);
	t91 = cos(t93);
	t161 = t91 * t164;
	t90 = sin(t93);
	t144 = -t111 * t90 + t161;
	t104 = t118 * t132 - t145 * t137;
	t88 = t91 * t111 + t164 * t90;
	t87 = 0.1e1 / t88 ^ 2;
	t159 = t104 * t87;
	t108 = 0.1e1 / t111;
	t109 = 0.1e1 / t111 ^ 2;
	t157 = t109 * t164;
	t92 = 0.1e1 / (t109 * t164 ^ 2 + 0.1e1);
	t169 = (t102 * t108 - t112 * t157) * t92;
	t158 = t104 ^ 2 * t87;
	t85 = 0.1e1 / (0.1e1 + t158);
	t86 = 0.1e1 / t88;
	t172 = (-t105 * t86 + (t102 * t90 + t91 * t112 + t144 * t169) * t159) * t85;
	t131 = sin(qJ(6));
	t136 = cos(qJ(6));
	t124 = -t130 * t149 - t148;
	t97 = t105 * t136 + t124 * t131;
	t95 = 0.1e1 / t97 ^ 2;
	t96 = t105 * t131 - t124 * t136;
	t160 = t95 * t96;
	t146 = t96 ^ 2 * t95 + 0.1e1;
	t89 = 0.1e1 / t146;
	t94 = 0.1e1 / t97;
	t162 = (-t131 * t94 + t136 * t160) * t89 * t104;
	t143 = t132 * t138 - t133 * t137;
	t122 = -t130 * t147 + t150;
	t119 = t143 * t139 * t129;
	t107 = (t132 * t133 + t137 * t138) * t124;
	t106 = t143 * t122;
	t84 = (t106 * t108 - t119 * t157) * t92;
	t1 = [-t104 * t108 * t92, t84, -t169, 0, t169, 0; (t164 * t86 - (-t90 + (-t108 * t161 + t90) * t92) * t158) * t85, (-(t90 * t106 + t91 * t119 + t144 * t84) * t159 + t143 * t86 * t124) * t85, t172, 0, -t172, 0; ((t102 * t131 - t122 * t136) * t94 - (t102 * t136 + t122 * t131) * t160) * t89, ((t107 * t131 + t125 * t136) * t94 - (t107 * t136 - t125 * t131) * t160) * t89, -t162, 0, t162, t146 * t89;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end