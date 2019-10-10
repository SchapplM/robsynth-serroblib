% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S7RRRRRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% Ja_rot [3x7]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 17:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S7RRRRRRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_rot_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_rot_sym_varpar: pkin has to be [4x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (63->18), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
	t38 = cos(qJ(2));
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t44 = t36 * t35;
	t30 = atan2(t44, t38);
	t27 = sin(t30);
	t28 = cos(t30);
	t20 = t27 * t44 + t28 * t38;
	t19 = 0.1e1 / t20 ^ 2;
	t39 = cos(qJ(1));
	t49 = t19 * t39 ^ 2;
	t34 = sin(qJ(3));
	t37 = cos(qJ(3));
	t41 = t39 * t37;
	t26 = -t36 * t34 + t38 * t41;
	t24 = 0.1e1 / t26 ^ 2;
	t42 = t39 * t34;
	t25 = t36 * t37 + t38 * t42;
	t48 = t24 * t25;
	t31 = t35 ^ 2;
	t47 = t31 / t38 ^ 2;
	t46 = t35 * t39;
	t29 = 0.1e1 / (t36 ^ 2 * t47 + 0.1e1);
	t45 = t36 * t29;
	t43 = t36 * t38;
	t40 = t25 ^ 2 * t24 + 0.1e1;
	t32 = 0.1e1 / t38;
	t23 = 0.1e1 / t26;
	t22 = (0.1e1 + t47) * t45;
	t21 = 0.1e1 / t40;
	t18 = 0.1e1 / t20;
	t17 = 0.1e1 / (t31 * t49 + 0.1e1);
	t1 = [t32 * t29 * t46, t22, 0, 0, 0, 0, 0; (t18 * t44 + (t28 * t31 * t32 * t45 + (-t29 + 0.1e1) * t35 * t27) * t35 * t49) * t17, (-t38 * t18 + (t27 * t43 - t28 * t35 + (-t27 * t38 + t28 * t44) * t22) * t35 * t19) * t39 * t17, 0, 0, 0, 0, 0; ((-t34 * t43 + t41) * t23 - (-t37 * t43 - t42) * t48) * t21, (-t23 * t34 + t37 * t48) * t21 * t46, t40 * t21, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (252->32), mult. (662->85), div. (108->11), fcn. (985->11), ass. (0->45)
	t56 = sin(qJ(3));
	t60 = cos(qJ(3));
	t62 = cos(qJ(1));
	t64 = t62 * t60;
	t58 = sin(qJ(1));
	t61 = cos(qJ(2));
	t66 = t58 * t61;
	t45 = t56 * t66 - t64;
	t57 = sin(qJ(2));
	t70 = t57 * t56;
	t43 = atan2(t45, -t70);
	t41 = sin(t43);
	t42 = cos(t43);
	t35 = t41 * t45 - t42 * t70;
	t33 = 0.1e1 / t35 ^ 2;
	t65 = t62 * t56;
	t47 = t58 * t60 + t61 * t65;
	t75 = t33 * t47;
	t49 = -t58 * t56 + t61 * t64;
	t55 = sin(qJ(4));
	t59 = cos(qJ(4));
	t67 = t57 * t62;
	t40 = t49 * t59 + t55 * t67;
	t38 = 0.1e1 / t40 ^ 2;
	t39 = t49 * t55 - t59 * t67;
	t74 = t38 * t39;
	t73 = t42 * t45;
	t72 = t47 ^ 2 * t33;
	t51 = 0.1e1 / t56;
	t53 = 0.1e1 / t57;
	t71 = t51 * t53;
	t69 = t57 * t58;
	t68 = t57 * t60;
	t63 = t39 ^ 2 * t38 + 0.1e1;
	t54 = 0.1e1 / t57 ^ 2;
	t52 = 0.1e1 / t56 ^ 2;
	t46 = t60 * t66 + t65;
	t44 = 0.1e1 / (t45 ^ 2 * t54 * t52 + 0.1e1);
	t37 = 0.1e1 / t40;
	t36 = 0.1e1 / t63;
	t34 = (t45 * t51 * t54 * t61 + t58) * t44;
	t32 = 0.1e1 / t35;
	t31 = 0.1e1 / (0.1e1 + t72);
	t30 = (t45 * t52 * t60 - t46 * t51) * t53 * t44;
	t1 = [-t47 * t44 * t71, t34, t30, 0, 0, 0, 0; (t45 * t32 + (t41 + (-t71 * t73 - t41) * t44) * t72) * t31, (t34 * t73 * t75 + (t32 * t67 + (-t42 * t61 + (t34 * t57 - t69) * t41) * t75) * t56) * t31, (-t49 * t32 + (-t42 * t68 + t41 * t46 + (t41 * t70 + t73) * t30) * t75) * t31, 0, 0, 0, 0; ((-t46 * t55 + t59 * t69) * t37 - (-t46 * t59 - t55 * t69) * t74) * t36, ((-t55 * t68 - t59 * t61) * t37 - (t55 * t61 - t59 * t68) * t74) * t36 * t62, (-t55 * t37 + t59 * t74) * t47 * t36, t63 * t36, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:05
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (609->45), mult. (1737->113), div. (115->9), fcn. (2479->13), ass. (0->60)
	t102 = sin(qJ(1));
	t101 = sin(qJ(2));
	t104 = cos(qJ(4));
	t117 = t101 * t104;
	t100 = sin(qJ(3));
	t107 = cos(qJ(1));
	t112 = t107 * t100;
	t105 = cos(qJ(3));
	t106 = cos(qJ(2));
	t113 = t105 * t106;
	t90 = t102 * t113 + t112;
	t99 = sin(qJ(4));
	t77 = -t102 * t117 + t90 * t99;
	t116 = t101 * t105;
	t87 = t106 * t104 + t99 * t116;
	t76 = atan2(-t77, t87);
	t73 = sin(t76);
	t74 = cos(t76);
	t67 = -t73 * t77 + t74 * t87;
	t66 = 0.1e1 / t67 ^ 2;
	t115 = t101 * t107;
	t111 = t107 * t105;
	t114 = t102 * t100;
	t93 = t106 * t111 - t114;
	t81 = -t104 * t115 + t93 * t99;
	t125 = t66 * t81;
	t103 = cos(qJ(5));
	t92 = -t102 * t105 - t106 * t112;
	t98 = sin(qJ(5));
	t120 = t92 * t98;
	t82 = t93 * t104 + t99 * t115;
	t72 = t82 * t103 + t120;
	t70 = 0.1e1 / t72 ^ 2;
	t119 = t92 * t103;
	t71 = t82 * t98 - t119;
	t124 = t70 * t71;
	t123 = t74 * t77;
	t86 = 0.1e1 / t87 ^ 2;
	t122 = t77 * t86;
	t121 = t81 ^ 2 * t66;
	t118 = t100 * t101;
	t110 = t71 ^ 2 * t70 + 0.1e1;
	t109 = t101 * t112;
	t108 = -t73 * t87 - t123;
	t79 = t102 * t101 * t99 + t90 * t104;
	t88 = t104 * t116 - t106 * t99;
	t91 = t99 * t113 - t117;
	t89 = t106 * t114 - t111;
	t85 = 0.1e1 / t87;
	t84 = t88 * t107;
	t83 = t87 * t102;
	t75 = 0.1e1 / (t77 ^ 2 * t86 + 0.1e1);
	t69 = 0.1e1 / t72;
	t68 = 0.1e1 / t110;
	t65 = 0.1e1 / t67;
	t64 = 0.1e1 / (0.1e1 + t121);
	t63 = (-t118 * t122 + t85 * t89) * t99 * t75;
	t62 = (t91 * t122 + t83 * t85) * t75;
	t61 = (t88 * t122 - t79 * t85) * t75;
	t1 = [-t81 * t85 * t75, t62, t63, t61, 0, 0, 0; (-t77 * t65 - (-t73 + (t85 * t123 + t73) * t75) * t121) * t64, (-(t108 * t62 + t73 * t83 + t74 * t91) * t125 - t87 * t65 * t107) * t64, (t92 * t99 * t65 - ((-t74 * t118 + t73 * t89) * t99 + t108 * t63) * t125) * t64, (t82 * t65 - (t108 * t61 - t73 * t79 + t74 * t88) * t125) * t64, 0, 0, 0; ((-t89 * t103 - t79 * t98) * t69 - (-t103 * t79 + t89 * t98) * t124) * t68, ((-t103 * t109 - t84 * t98) * t69 - (-t84 * t103 + t98 * t109) * t124) * t68, ((t93 * t103 + t104 * t120) * t69 - (t104 * t119 - t93 * t98) * t124) * t68, (t103 * t124 - t98 * t69) * t81 * t68, t110 * t68, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:05
	% EndTime: 2019-10-10 17:10:05
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (1373->65), mult. (3817->158), div. (145->9), fcn. (5357->15), ass. (0->78)
	t141 = cos(qJ(3));
	t135 = sin(qJ(3));
	t143 = cos(qJ(1));
	t152 = t143 * t135;
	t137 = sin(qJ(1));
	t142 = cos(qJ(2));
	t154 = t137 * t142;
	t125 = t141 * t154 + t152;
	t134 = sin(qJ(4));
	t140 = cos(qJ(4));
	t136 = sin(qJ(2));
	t157 = t136 * t137;
	t112 = t125 * t140 + t134 * t157;
	t133 = sin(qJ(5));
	t139 = cos(qJ(5));
	t151 = t143 * t141;
	t146 = -t135 * t154 + t151;
	t97 = t112 * t133 - t146 * t139;
	t144 = t146 * t133;
	t99 = t112 * t139 + t144;
	t127 = -t137 * t135 + t142 * t151;
	t155 = t136 * t143;
	t117 = t127 * t140 + t134 * t155;
	t126 = -t137 * t141 - t142 * t152;
	t102 = t117 * t139 + t126 * t133;
	t116 = t127 * t134 - t140 * t155;
	t132 = sin(qJ(6));
	t138 = cos(qJ(6));
	t92 = t102 * t138 + t116 * t132;
	t90 = 0.1e1 / t92 ^ 2;
	t91 = t102 * t132 - t116 * t138;
	t167 = t90 * t91;
	t156 = t136 * t141;
	t124 = -t142 * t134 + t140 * t156;
	t158 = t135 * t139;
	t148 = t136 * t158;
	t109 = t124 * t133 + t148;
	t96 = atan2(-t97, t109);
	t94 = cos(t96);
	t166 = t94 * t97;
	t160 = t126 * t139;
	t101 = t117 * t133 - t160;
	t93 = sin(t96);
	t87 = t94 * t109 - t93 * t97;
	t86 = 0.1e1 / t87 ^ 2;
	t165 = t101 * t86;
	t164 = t101 ^ 2 * t86;
	t108 = 0.1e1 / t109 ^ 2;
	t163 = t108 * t97;
	t162 = t116 * t139;
	t161 = t126 * t134;
	t159 = t133 * t140;
	t153 = t142 * t140;
	t150 = t91 ^ 2 * t90 + 0.1e1;
	t149 = t136 * t135 * t133;
	t147 = -t125 * t134 + t140 * t157;
	t145 = -t109 * t93 - t166;
	t123 = -t134 * t156 - t153;
	t120 = t124 * t143;
	t119 = t123 * t143;
	t118 = (-t135 * t159 + t139 * t141) * t136;
	t115 = (t136 * t134 + t141 * t153) * t133 + t142 * t158;
	t110 = t124 * t139 - t149;
	t107 = 0.1e1 / t109;
	t106 = -t120 * t139 + t143 * t149;
	t105 = t109 * t137;
	t104 = -t127 * t133 + t140 * t160;
	t103 = t125 * t139 + t140 * t144;
	t95 = 0.1e1 / (t97 ^ 2 * t108 + 0.1e1);
	t89 = 0.1e1 / t92;
	t88 = 0.1e1 / t150;
	t85 = 0.1e1 / t87;
	t84 = 0.1e1 / (0.1e1 + t164);
	t83 = (-t107 * t147 + t123 * t163) * t95 * t133;
	t82 = (-t103 * t107 + t118 * t163) * t95;
	t81 = (t105 * t107 + t115 * t163) * t95;
	t80 = (-t107 * t99 + t110 * t163) * t95;
	t1 = [-t101 * t107 * t95, t81, t82, t83, t80, 0, 0; (-t97 * t85 - (-t93 + (t107 * t166 + t93) * t95) * t164) * t84, ((-t120 * t133 - t143 * t148) * t85 - (t93 * t105 + t94 * t115 + t145 * t81) * t165) * t84, ((t126 * t159 + t127 * t139) * t85 - (-t93 * t103 + t94 * t118 + t145 * t82) * t165) * t84, (-t116 * t133 * t85 - (t145 * t83 + (t123 * t94 - t147 * t93) * t133) * t165) * t84, (t102 * t85 - (t94 * t110 + t145 * t80 - t93 * t99) * t165) * t84, 0, 0; ((-t132 * t99 - t138 * t147) * t89 - (t132 * t147 - t138 * t99) * t167) * t88, ((t106 * t132 - t119 * t138) * t89 - (t106 * t138 + t119 * t132) * t167) * t88, ((t104 * t132 - t138 * t161) * t89 - (t104 * t138 + t132 * t161) * t167) * t88, ((-t117 * t138 - t132 * t162) * t89 - (t117 * t132 - t138 * t162) * t167) * t88, (-t132 * t89 + t138 * t167) * t88 * t101, t150 * t88, 0;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:05
	% EndTime: 2019-10-10 17:10:07
	% DurationCPUTime: 1.50s
	% Computational Cost: add. (2869->90), mult. (7788->218), div. (175->9), fcn. (10858->17), ass. (0->100)
	t196 = sin(qJ(1));
	t200 = cos(qJ(4));
	t201 = cos(qJ(3));
	t194 = sin(qJ(3));
	t237 = cos(qJ(1));
	t214 = t237 * t194;
	t202 = cos(qJ(2));
	t221 = t196 * t202;
	t207 = t201 * t221 + t214;
	t193 = sin(qJ(4));
	t195 = sin(qJ(2));
	t224 = t195 * t193;
	t175 = t196 * t224 + t207 * t200;
	t213 = t237 * t201;
	t184 = t194 * t221 - t213;
	t192 = sin(qJ(5));
	t199 = cos(qJ(5));
	t159 = t175 * t199 - t184 * t192;
	t191 = sin(qJ(6));
	t198 = cos(qJ(6));
	t223 = t195 * t200;
	t203 = t207 * t193 - t196 * t223;
	t142 = t159 * t191 - t203 * t198;
	t143 = t159 * t198 + t203 * t191;
	t238 = -t175 * t192 - t184 * t199;
	t219 = t202 * t193;
	t222 = t195 * t201;
	t183 = t200 * t222 - t219;
	t225 = t194 * t195;
	t174 = t183 * t199 - t192 * t225;
	t218 = t202 * t200;
	t182 = t193 * t222 + t218;
	t156 = -t174 * t191 + t182 * t198;
	t141 = atan2(t142, t156);
	t138 = sin(t141);
	t139 = cos(t141);
	t132 = t138 * t142 + t139 * t156;
	t131 = 0.1e1 / t132 ^ 2;
	t208 = -t196 * t194 + t202 * t213;
	t216 = t195 * t237;
	t210 = t193 * t216;
	t178 = t208 * t200 + t210;
	t185 = -t196 * t201 - t202 * t214;
	t228 = t185 * t192;
	t164 = t178 * t199 + t228;
	t187 = t200 * t216;
	t206 = t208 * t193 - t187;
	t145 = t164 * t191 - t206 * t198;
	t236 = t131 * t145;
	t147 = t164 * t198 + t206 * t191;
	t163 = t178 * t192 - t185 * t199;
	t190 = sin(qJ(7));
	t197 = cos(qJ(7));
	t137 = t147 * t197 - t163 * t190;
	t135 = 0.1e1 / t137 ^ 2;
	t136 = t147 * t190 + t163 * t197;
	t235 = t135 * t136;
	t234 = t139 * t142;
	t155 = 0.1e1 / t156 ^ 2;
	t233 = t142 * t155;
	t232 = t145 ^ 2 * t131;
	t231 = t163 * t198;
	t227 = t191 * t199;
	t226 = t193 * t198;
	t220 = t200 * t199;
	t215 = t202 * t237;
	t212 = t136 ^ 2 * t135 + 0.1e1;
	t211 = t195 * t214;
	t209 = -t138 * t156 + t234;
	t205 = t192 * t206;
	t204 = t199 * t206;
	t180 = -t201 * t187 + t193 * t215;
	t179 = -t200 * t215 - t201 * t210;
	t173 = -t183 * t192 - t199 * t225;
	t170 = t180 * t199 + t192 * t211;
	t169 = t180 * t192 - t199 * t211;
	t168 = (-(-t192 * t201 - t194 * t220) * t191 - t194 * t226) * t195;
	t167 = t185 * t220 - t208 * t192;
	t166 = t208 * t199 + t200 * t228;
	t165 = t182 * t227 + t183 * t198;
	t162 = -((t201 * t218 + t224) * t199 - t202 * t194 * t192) * t191 + (t201 * t219 - t223) * t198;
	t157 = -t174 * t198 - t182 * t191;
	t154 = 0.1e1 / t156;
	t153 = t170 * t198 + t179 * t191;
	t152 = t156 * t196;
	t151 = t185 * t193 * t191 + t167 * t198;
	t150 = (-t184 * t220 - t207 * t192) * t191 + t184 * t226;
	t149 = t178 * t191 - t198 * t204;
	t148 = -t175 * t198 - t203 * t227;
	t140 = 0.1e1 / (t142 ^ 2 * t155 + 0.1e1);
	t134 = 0.1e1 / t137;
	t133 = 0.1e1 / t212;
	t130 = 0.1e1 / t132;
	t129 = 0.1e1 / (0.1e1 + t232);
	t128 = (t154 * t238 + t173 * t233) * t191 * t140;
	t127 = (t150 * t154 - t168 * t233) * t140;
	t126 = (t148 * t154 - t165 * t233) * t140;
	t125 = (t152 * t154 - t162 * t233) * t140;
	t124 = (t143 * t154 - t157 * t233) * t140;
	t1 = [t145 * t154 * t140, t125, t127, t126, t128, t124, 0; (t142 * t130 + (t138 + (t154 * t234 - t138) * t140) * t232) * t129, ((-t170 * t191 + t179 * t198) * t130 + (t209 * t125 + t138 * t152 + t139 * t162) * t236) * t129, ((-t167 * t191 + t185 * t226) * t130 + (t209 * t127 + t138 * t150 + t139 * t168) * t236) * t129, ((t178 * t198 + t191 * t204) * t130 + (t209 * t126 + t138 * t148 + t139 * t165) * t236) * t129, (t163 * t191 * t130 + ((t138 * t238 - t139 * t173) * t191 + t209 * t128) * t236) * t129, (-t147 * t130 + (t209 * t124 + t138 * t143 + t139 * t157) * t236) * t129, 0; ((-t143 * t190 + t197 * t238) * t134 - (-t143 * t197 - t190 * t238) * t235) * t133, ((t153 * t190 + t169 * t197) * t134 - (t153 * t197 - t169 * t190) * t235) * t133, ((t151 * t190 + t166 * t197) * t134 - (t151 * t197 - t166 * t190) * t235) * t133, ((t149 * t190 - t197 * t205) * t134 - (t149 * t197 + t190 * t205) * t235) * t133, ((t164 * t197 - t190 * t231) * t134 - (-t164 * t190 - t197 * t231) * t235) * t133, (-t190 * t134 + t197 * t235) * t145 * t133, t212 * t133;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,7);
end