% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRP8
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
%   Wie in S6RRRRRP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
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
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
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
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (252->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
	t69 = sin(qJ(2));
	t71 = cos(qJ(2));
	t72 = cos(qJ(1));
	t70 = sin(qJ(1));
	t76 = cos(pkin(6));
	t74 = t70 * t76;
	t60 = -t69 * t74 + t72 * t71;
	t67 = qJ(3) + qJ(4);
	t62 = sin(t67);
	t63 = cos(t67);
	t68 = sin(pkin(6));
	t79 = t68 * t70;
	t51 = t60 * t63 + t62 * t79;
	t49 = 0.1e1 / t51 ^ 2;
	t50 = t60 * t62 - t63 * t79;
	t83 = t49 * t50;
	t73 = t72 * t76;
	t56 = t70 * t69 - t71 * t73;
	t78 = t68 * t71;
	t54 = atan2(-t56, -t78);
	t53 = cos(t54);
	t82 = t53 * t56;
	t52 = sin(t54);
	t46 = -t52 * t56 - t53 * t78;
	t45 = 0.1e1 / t46 ^ 2;
	t59 = t72 * t69 + t71 * t74;
	t81 = t59 ^ 2 * t45;
	t64 = 0.1e1 / t68;
	t65 = 0.1e1 / t71;
	t80 = t64 * t65;
	t77 = t68 * t72;
	t75 = t50 ^ 2 * t49 + 0.1e1;
	t66 = 0.1e1 / t71 ^ 2;
	t58 = t69 * t73 + t70 * t71;
	t55 = 0.1e1 / (0.1e1 + t56 ^ 2 / t68 ^ 2 * t66);
	t48 = 0.1e1 / t51;
	t47 = 0.1e1 / t75;
	t44 = 0.1e1 / t46;
	t43 = 0.1e1 / (0.1e1 + t81);
	t42 = (t56 * t66 * t69 + t58 * t65) * t64 * t55;
	t41 = t75 * t47;
	t1 = [t59 * t55 * t80, t42, 0, 0, 0, 0; (-t56 * t44 - (-t52 + (-t80 * t82 + t52) * t55) * t81) * t43, (t60 * t44 - (t53 * t68 * t69 - t52 * t58 + (t52 * t78 - t82) * t42) * t59 * t45) * t43, 0, 0, 0, 0; ((-t58 * t62 - t63 * t77) * t48 - (-t58 * t63 + t62 * t77) * t83) * t47, (-t62 * t48 + t63 * t83) * t59 * t47, t41, t41, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (1185->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
	t100 = sin(pkin(6));
	t107 = cos(qJ(1));
	t114 = t100 * t107;
	t101 = cos(pkin(6));
	t103 = sin(qJ(2));
	t111 = t107 * t103;
	t104 = sin(qJ(1));
	t106 = cos(qJ(2));
	t112 = t104 * t106;
	t92 = t101 * t111 + t112;
	t99 = qJ(3) + qJ(4);
	t97 = sin(t99);
	t98 = cos(t99);
	t81 = t98 * t114 + t92 * t97;
	t117 = t100 * t103;
	t89 = -t101 * t98 + t97 * t117;
	t80 = atan2(-t81, t89);
	t75 = sin(t80);
	t76 = cos(t80);
	t71 = -t75 * t81 + t76 * t89;
	t70 = 0.1e1 / t71 ^ 2;
	t116 = t100 * t104;
	t110 = t107 * t106;
	t113 = t104 * t103;
	t94 = -t101 * t113 + t110;
	t85 = -t98 * t116 + t94 * t97;
	t124 = t70 * t85;
	t105 = cos(qJ(5));
	t102 = sin(qJ(5));
	t93 = t101 * t112 + t111;
	t119 = t93 * t102;
	t86 = t97 * t116 + t94 * t98;
	t78 = t86 * t105 + t119;
	t74 = 0.1e1 / t78 ^ 2;
	t118 = t93 * t105;
	t77 = t86 * t102 - t118;
	t123 = t74 * t77;
	t122 = t76 * t81;
	t88 = 0.1e1 / t89 ^ 2;
	t121 = t81 * t88;
	t120 = t85 ^ 2 * t70;
	t115 = t100 * t106;
	t109 = t77 ^ 2 * t74 + 0.1e1;
	t83 = -t97 * t114 + t92 * t98;
	t108 = -t75 * t89 - t122;
	t91 = t101 * t110 - t113;
	t90 = t101 * t97 + t98 * t117;
	t87 = 0.1e1 / t89;
	t79 = 0.1e1 / (t81 ^ 2 * t88 + 0.1e1);
	t73 = 0.1e1 / t78;
	t72 = 0.1e1 / t109;
	t69 = 0.1e1 / t71;
	t68 = 0.1e1 / (0.1e1 + t120);
	t67 = (t115 * t121 - t87 * t91) * t97 * t79;
	t66 = (t90 * t121 - t83 * t87) * t79;
	t65 = (-t102 * t73 + t105 * t123) * t85 * t72;
	t64 = (t86 * t69 - (t108 * t66 - t75 * t83 + t76 * t90) * t124) * t68;
	t1 = [-t85 * t87 * t79, t67, t66, t66, 0, 0; (-t81 * t69 - (-t75 + (t87 * t122 + t75) * t79) * t120) * t68, (-t93 * t97 * t69 - ((t76 * t115 - t75 * t91) * t97 + t108 * t67) * t124) * t68, t64, t64, 0, 0; ((-t102 * t83 - t91 * t105) * t73 - (t91 * t102 - t105 * t83) * t123) * t72, ((-t94 * t105 - t98 * t119) * t73 - (t94 * t102 - t98 * t118) * t123) * t72, t65, t65, t109 * t72, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:42
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (1836->47), mult. (3176->112), div. (137->9), fcn. (4481->13), ass. (0->64)
	t140 = cos(pkin(6));
	t142 = sin(qJ(2));
	t146 = cos(qJ(1));
	t151 = t146 * t142;
	t143 = sin(qJ(1));
	t145 = cos(qJ(2));
	t152 = t143 * t145;
	t132 = t140 * t151 + t152;
	t138 = qJ(3) + qJ(4);
	t136 = sin(t138);
	t137 = cos(t138);
	t139 = sin(pkin(6));
	t155 = t139 * t146;
	t124 = -t132 * t137 + t136 * t155;
	t141 = sin(qJ(5));
	t144 = cos(qJ(5));
	t150 = t146 * t145;
	t153 = t143 * t142;
	t149 = -t140 * t150 + t153;
	t167 = t124 * t141 + t149 * t144;
	t148 = t149 * t141;
	t166 = t124 * t144 - t148;
	t157 = t139 * t142;
	t129 = t140 * t136 + t137 * t157;
	t120 = t139 * t145 * t144 + t129 * t141;
	t108 = atan2(t167, t120);
	t104 = sin(t108);
	t105 = cos(t108);
	t103 = t104 * t167 + t105 * t120;
	t102 = 0.1e1 / t103 ^ 2;
	t134 = -t140 * t153 + t150;
	t156 = t139 * t143;
	t126 = t134 * t137 + t136 * t156;
	t133 = t140 * t152 + t151;
	t158 = t133 * t144;
	t114 = t126 * t141 - t158;
	t164 = t102 * t114;
	t159 = t133 * t141;
	t115 = t126 * t144 + t159;
	t110 = 0.1e1 / t115 ^ 2;
	t125 = -t134 * t136 + t137 * t156;
	t163 = t110 * t125;
	t118 = 0.1e1 / t120 ^ 2;
	t162 = t167 * t118;
	t161 = t114 ^ 2 * t102;
	t160 = t125 ^ 2 * t110;
	t154 = t141 * t145;
	t147 = t132 * t136 + t137 * t155;
	t128 = -t136 * t157 + t140 * t137;
	t127 = (t137 * t154 - t142 * t144) * t139;
	t121 = t129 * t144 - t139 * t154;
	t117 = 0.1e1 / t120;
	t116 = -t132 * t144 - t137 * t148;
	t109 = 0.1e1 / t115;
	t107 = 0.1e1 / (0.1e1 + t160);
	t106 = 0.1e1 / (t118 * t167 ^ 2 + 0.1e1);
	t101 = 0.1e1 / t103;
	t100 = 0.1e1 / (0.1e1 + t161);
	t99 = (-t109 * t126 - t144 * t160) * t107;
	t98 = (t117 * t147 - t128 * t162) * t141 * t106;
	t97 = (-t116 * t117 - t127 * t162) * t106;
	t96 = (t117 * t166 - t121 * t162) * t106;
	t95 = (t125 * t141 * t101 - ((t128 * t141 + t167 * t98) * t105 + (-t120 * t98 + t141 * t147) * t104) * t164) * t100;
	t1 = [-t114 * t117 * t106, t97, t98, t98, t96, 0; (t167 * t101 - (-t104 + (-t105 * t117 * t167 + t104) * t106) * t161) * t100, ((-t134 * t144 - t137 * t159) * t101 - ((t167 * t97 + t127) * t105 + (-t120 * t97 - t116) * t104) * t164) * t100, t95, t95, (t115 * t101 - ((t167 * t96 + t121) * t105 + (-t120 * t96 + t166) * t104) * t164) * t100, 0; (t147 * t109 - t166 * t163) * t107, (t133 * t136 * t109 - (t134 * t141 - t137 * t158) * t163) * t107, t99, t99, t114 * t107 * t163, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end