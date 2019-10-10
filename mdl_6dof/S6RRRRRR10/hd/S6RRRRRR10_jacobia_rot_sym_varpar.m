% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR10
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
%   Wie in S6RRRRRR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
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
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:53
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (304->29), mult. (885->75), div. (55->9), fcn. (1254->13), ass. (0->49)
	t69 = cos(pkin(6));
	t74 = cos(qJ(2));
	t75 = cos(qJ(1));
	t79 = t74 * t75;
	t71 = sin(qJ(2));
	t72 = sin(qJ(1));
	t81 = t72 * t71;
	t60 = -t69 * t79 + t81;
	t66 = sin(pkin(7));
	t68 = cos(pkin(7));
	t67 = sin(pkin(6));
	t83 = t67 * t75;
	t54 = -t60 * t66 + t68 * t83;
	t59 = -t66 * t67 * t74 + t68 * t69;
	t53 = atan2(t54, t59);
	t50 = sin(t53);
	t51 = cos(t53);
	t44 = t50 * t54 + t51 * t59;
	t43 = 0.1e1 / t44 ^ 2;
	t80 = t72 * t74;
	t82 = t71 * t75;
	t62 = -t69 * t80 - t82;
	t84 = t67 * t72;
	t55 = t62 * t66 - t68 * t84;
	t90 = t43 * t55 ^ 2;
	t70 = sin(qJ(3));
	t76 = t62 * t68 + t66 * t84;
	t63 = -t69 * t81 + t79;
	t73 = cos(qJ(3));
	t86 = t63 * t73;
	t49 = t76 * t70 + t86;
	t47 = 0.1e1 / t49 ^ 2;
	t87 = t63 * t70;
	t48 = -t76 * t73 + t87;
	t89 = t47 * t48;
	t88 = t51 * t54;
	t85 = t67 * t71;
	t78 = t47 * t48 ^ 2 + 0.1e1;
	t77 = t60 * t68 + t66 * t83;
	t61 = -t69 * t82 - t80;
	t58 = 0.1e1 / t59 ^ 2;
	t57 = 0.1e1 / t59;
	t52 = 0.1e1 / (t54 ^ 2 * t58 + 0.1e1);
	t46 = 0.1e1 / t49;
	t45 = 0.1e1 / t78;
	t42 = 0.1e1 / t44;
	t41 = 0.1e1 / (0.1e1 + t90);
	t40 = (-t54 * t58 * t85 + t57 * t61) * t66 * t52;
	t1 = [t55 * t57 * t52, t40, 0, 0, 0, 0; (t54 * t42 + (t50 + (t57 * t88 - t50) * t52) * t90) * t41, (t63 * t66 * t42 + ((t50 * t61 + t51 * t85) * t66 + (-t50 * t59 + t88) * t40) * t55 * t43) * t41, 0, 0, 0, 0; ((t61 * t70 - t77 * t73) * t46 - (t61 * t73 + t77 * t70) * t89) * t45, ((t62 * t70 + t68 * t86) * t46 - (t62 * t73 - t68 * t87) * t89) * t45, t78 * t45, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:53
	% EndTime: 2019-10-10 13:34:54
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (1426->58), mult. (4126->139), div. (85->9), fcn. (5651->17), ass. (0->76)
	t128 = cos(pkin(6));
	t131 = sin(qJ(2));
	t167 = sin(qJ(1));
	t145 = t167 * t131;
	t134 = cos(qJ(2));
	t135 = cos(qJ(1));
	t149 = t135 * t134;
	t117 = -t128 * t149 + t145;
	t144 = t167 * t134;
	t150 = t135 * t131;
	t118 = t128 * t150 + t144;
	t127 = cos(pkin(7));
	t130 = sin(qJ(3));
	t133 = cos(qJ(3));
	t124 = sin(pkin(7));
	t125 = sin(pkin(6));
	t155 = t125 * t135;
	t147 = t124 * t155;
	t106 = (t117 * t127 + t147) * t133 + t118 * t130;
	t116 = -t117 * t124 + t127 * t155;
	t123 = sin(pkin(8));
	t126 = cos(pkin(8));
	t98 = -t106 * t123 + t116 * t126;
	t153 = t127 * t133;
	t156 = t124 * t128;
	t104 = -(t133 * t156 + (-t130 * t131 + t134 * t153) * t125) * t123 + (-t125 * t134 * t124 + t128 * t127) * t126;
	t97 = atan2(t98, t104);
	t94 = sin(t97);
	t95 = cos(t97);
	t88 = t95 * t104 + t94 * t98;
	t87 = 0.1e1 / t88 ^ 2;
	t120 = -t128 * t145 + t149;
	t119 = -t128 * t144 - t150;
	t146 = t125 * t167;
	t137 = t119 * t127 + t124 * t146;
	t108 = -t120 * t130 + t137 * t133;
	t138 = t119 * t124 - t127 * t146;
	t99 = t108 * t123 + t138 * t126;
	t166 = t87 * t99;
	t129 = sin(qJ(4));
	t136 = t108 * t126 - t138 * t123;
	t109 = t120 * t133 + t137 * t130;
	t132 = cos(qJ(4));
	t160 = t109 * t132;
	t93 = t136 * t129 + t160;
	t91 = 0.1e1 / t93 ^ 2;
	t161 = t109 * t129;
	t92 = -t136 * t132 + t161;
	t165 = t91 * t92;
	t164 = t95 * t98;
	t163 = t99 ^ 2 * t87;
	t103 = 0.1e1 / t104 ^ 2;
	t162 = t103 * t98;
	t157 = t124 * t126;
	t154 = t127 * t130;
	t152 = t130 * t134;
	t151 = t131 * t133;
	t148 = t92 ^ 2 * t91 + 0.1e1;
	t143 = -t104 * t94 + t164;
	t142 = t106 * t126 + t116 * t123;
	t110 = -t119 * t130 - t120 * t153;
	t141 = t120 * t123 * t124 + t110 * t126;
	t139 = t117 * t154 - t118 * t133 + t130 * t147;
	t115 = -t130 * t156 + (-t127 * t152 - t151) * t125;
	t112 = (-(-t127 * t151 - t152) * t123 + t131 * t157) * t125;
	t111 = t119 * t133 - t120 * t154;
	t102 = 0.1e1 / t104;
	t101 = (t117 * t130 - t118 * t153) * t123 - t118 * t157;
	t96 = 0.1e1 / (t98 ^ 2 * t103 + 0.1e1);
	t90 = 0.1e1 / t93;
	t89 = 0.1e1 / t148;
	t86 = 0.1e1 / t88;
	t85 = 0.1e1 / (0.1e1 + t163);
	t84 = (t102 * t139 + t115 * t162) * t96 * t123;
	t83 = (t101 * t102 - t112 * t162) * t96;
	t1 = [t99 * t102 * t96, t83, t84, 0, 0, 0; (t98 * t86 + (t94 + (t102 * t164 - t94) * t96) * t163) * t85, ((-t110 * t123 + t120 * t157) * t86 + (t94 * t101 + t95 * t112 + t143 * t83) * t166) * t85, (t109 * t123 * t86 + (t143 * t84 + (-t115 * t95 + t139 * t94) * t123) * t166) * t85, 0, 0, 0; ((t129 * t139 - t142 * t132) * t90 - (t142 * t129 + t132 * t139) * t165) * t89, ((t111 * t129 - t141 * t132) * t90 - (t111 * t132 + t141 * t129) * t165) * t89, ((t108 * t129 + t126 * t160) * t90 - (t108 * t132 - t126 * t161) * t165) * t89, t148 * t89, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:54
	% EndTime: 2019-10-10 13:34:55
	% DurationCPUTime: 1.34s
	% Computational Cost: add. (3237->82), mult. (9534->184), div. (115->9), fcn. (12904->19), ass. (0->95)
	t177 = cos(pkin(6));
	t181 = sin(qJ(2));
	t226 = sin(qJ(1));
	t202 = t226 * t181;
	t185 = cos(qJ(2));
	t186 = cos(qJ(1));
	t205 = t186 * t185;
	t166 = -t177 * t205 + t202;
	t201 = t226 * t185;
	t206 = t186 * t181;
	t167 = t177 * t206 + t201;
	t180 = sin(qJ(3));
	t184 = cos(qJ(3));
	t173 = sin(pkin(7));
	t174 = sin(pkin(6));
	t214 = t174 * t186;
	t204 = t173 * t214;
	t176 = cos(pkin(7));
	t212 = t176 * t180;
	t156 = t166 * t212 - t167 * t184 + t180 * t204;
	t179 = sin(qJ(4));
	t183 = cos(qJ(4));
	t155 = (t166 * t176 + t204) * t184 + t167 * t180;
	t164 = -t166 * t173 + t176 * t214;
	t172 = sin(pkin(8));
	t175 = cos(pkin(8));
	t198 = t155 * t175 + t164 * t172;
	t136 = t156 * t183 + t198 * t179;
	t133 = -t156 * t179 + t198 * t183;
	t208 = t181 * t184;
	t209 = t180 * t185;
	t215 = t173 * t177;
	t163 = t180 * t215 + (t176 * t209 + t208) * t174;
	t207 = t184 * t185;
	t210 = t180 * t181;
	t162 = t184 * t215 + (t176 * t207 - t210) * t174;
	t197 = t162 * t175 + (-t174 * t185 * t173 + t177 * t176) * t172;
	t145 = t163 * t179 - t197 * t183;
	t132 = atan2(-t133, t145);
	t127 = sin(t132);
	t128 = cos(t132);
	t123 = -t127 * t133 + t128 * t145;
	t122 = 0.1e1 / t123 ^ 2;
	t169 = -t177 * t202 + t205;
	t168 = -t177 * t201 - t206;
	t203 = t174 * t226;
	t193 = t168 * t176 + t173 * t203;
	t189 = t169 * t180 - t193 * t184;
	t187 = t189 * t183;
	t194 = -t168 * t173 + t176 * t203;
	t191 = t194 * t172;
	t157 = t169 * t184 + t193 * t180;
	t219 = t157 * t179;
	t137 = t175 * t187 - t183 * t191 + t219;
	t225 = t122 * t137;
	t138 = t157 * t183 + (-t189 * t175 + t191) * t179;
	t148 = t189 * t172 + t194 * t175;
	t178 = sin(qJ(5));
	t182 = cos(qJ(5));
	t130 = t138 * t182 + t148 * t178;
	t126 = 0.1e1 / t130 ^ 2;
	t129 = t138 * t178 - t148 * t182;
	t224 = t126 * t129;
	t223 = t128 * t133;
	t144 = 0.1e1 / t145 ^ 2;
	t222 = t133 * t144;
	t221 = t137 ^ 2 * t122;
	t220 = t157 * t172;
	t216 = t173 * t172;
	t213 = t175 * t183;
	t211 = t176 * t184;
	t200 = t129 ^ 2 * t126 + 0.1e1;
	t199 = -t127 * t145 - t223;
	t158 = -t168 * t180 - t169 * t211;
	t196 = t158 * t175 + t169 * t216;
	t159 = t168 * t184 - t169 * t212;
	t151 = t169 * t173 * t175 - t158 * t172;
	t150 = ((-t176 * t210 + t207) * t179 + (-(-t176 * t208 - t209) * t175 - t181 * t216) * t183) * t174;
	t149 = t162 * t179 + t163 * t213;
	t147 = -t155 * t172 + t164 * t175;
	t146 = t163 * t183 + t197 * t179;
	t143 = 0.1e1 / t145;
	t142 = -t175 * t219 - t187;
	t141 = -t155 * t179 - t156 * t213;
	t140 = t159 * t183 + t196 * t179;
	t139 = (-t166 * t184 - t167 * t212) * t179 + (-(t166 * t180 - t167 * t211) * t175 - t167 * t216) * t183;
	t131 = 0.1e1 / (t133 ^ 2 * t144 + 0.1e1);
	t125 = 0.1e1 / t130;
	t124 = 0.1e1 / t200;
	t121 = 0.1e1 / t123;
	t120 = 0.1e1 / (0.1e1 + t221);
	t119 = (-t139 * t143 + t150 * t222) * t131;
	t118 = (-t141 * t143 + t149 * t222) * t131;
	t117 = (t136 * t143 + t146 * t222) * t131;
	t1 = [-t137 * t143 * t131, t119, t118, t117, 0, 0; (-t133 * t121 - (-t127 + (t143 * t223 + t127) * t131) * t221) * t120, ((t159 * t179 - t196 * t183) * t121 - (t199 * t119 - t127 * t139 + t128 * t150) * t225) * t120, ((t157 * t213 - t189 * t179) * t121 - (t199 * t118 - t127 * t141 + t128 * t149) * t225) * t120, (t138 * t121 - (t199 * t117 + t127 * t136 + t128 * t146) * t225) * t120, 0, 0; ((t136 * t178 - t147 * t182) * t125 - (t136 * t182 + t147 * t178) * t224) * t124, ((t140 * t178 - t151 * t182) * t125 - (t140 * t182 + t151 * t178) * t224) * t124, ((t142 * t178 - t182 * t220) * t125 - (t142 * t182 + t178 * t220) * t224) * t124, (-t178 * t125 + t182 * t224) * t137 * t124, t200 * t124, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:34:55
	% EndTime: 2019-10-10 13:34:57
	% DurationCPUTime: 2.66s
	% Computational Cost: add. (6679->105), mult. (19296->240), div. (145->9), fcn. (26142->21), ass. (0->114)
	t227 = cos(pkin(6));
	t232 = sin(qJ(2));
	t284 = sin(qJ(1));
	t256 = t284 * t232;
	t237 = cos(qJ(2));
	t238 = cos(qJ(1));
	t259 = t238 * t237;
	t216 = -t227 * t259 + t256;
	t255 = t284 * t237;
	t260 = t238 * t232;
	t217 = t227 * t260 + t255;
	t231 = sin(qJ(3));
	t236 = cos(qJ(3));
	t223 = sin(pkin(7));
	t224 = sin(pkin(6));
	t269 = t224 * t238;
	t258 = t223 * t269;
	t226 = cos(pkin(7));
	t267 = t226 * t231;
	t205 = t216 * t267 - t217 * t236 + t231 * t258;
	t230 = sin(qJ(4));
	t235 = cos(qJ(4));
	t204 = (t216 * t226 + t258) * t236 + t217 * t231;
	t213 = -t216 * t223 + t226 * t269;
	t222 = sin(pkin(8));
	t225 = cos(pkin(8));
	t252 = t204 * t225 + t213 * t222;
	t186 = t205 * t235 + t230 * t252;
	t229 = sin(qJ(5));
	t234 = cos(qJ(5));
	t241 = t204 * t222 - t213 * t225;
	t290 = t186 * t229 + t234 * t241;
	t172 = t186 * t234 - t229 * t241;
	t287 = t205 * t230 - t235 * t252;
	t262 = t232 * t236;
	t264 = t231 * t237;
	t270 = t223 * t227;
	t212 = t231 * t270 + (t226 * t264 + t262) * t224;
	t261 = t236 * t237;
	t265 = t231 * t232;
	t211 = t236 * t270 + (t226 * t261 - t265) * t224;
	t215 = -t224 * t237 * t223 + t227 * t226;
	t251 = t211 * t225 + t215 * t222;
	t196 = t212 * t235 + t230 * t251;
	t202 = -t211 * t222 + t215 * t225;
	t181 = t196 * t229 - t202 * t234;
	t168 = atan2(t290, t181);
	t161 = sin(t168);
	t162 = cos(t168);
	t159 = t161 * t290 + t162 * t181;
	t158 = 0.1e1 / t159 ^ 2;
	t218 = -t227 * t256 + t259;
	t248 = t227 * t255 + t260;
	t257 = t224 * t284;
	t245 = t223 * t257 - t248 * t226;
	t243 = t218 * t231 - t245 * t236;
	t246 = t248 * t223 + t226 * t257;
	t244 = t246 * t222;
	t206 = t218 * t236 + t245 * t231;
	t277 = t206 * t235;
	t188 = t277 + (-t225 * t243 + t244) * t230;
	t239 = t243 * t222 + t246 * t225;
	t173 = t188 * t229 - t234 * t239;
	t283 = t158 * t173;
	t282 = t162 * t290;
	t174 = t188 * t234 + t229 * t239;
	t242 = t243 * t235;
	t187 = t206 * t230 + t225 * t242 - t235 * t244;
	t228 = sin(qJ(6));
	t233 = cos(qJ(6));
	t167 = t174 * t233 + t187 * t228;
	t164 = 0.1e1 / t167 ^ 2;
	t166 = t174 * t228 - t187 * t233;
	t281 = t164 * t166;
	t180 = 0.1e1 / t181 ^ 2;
	t280 = t290 * t180;
	t279 = t173 ^ 2 * t158;
	t278 = t187 * t234;
	t273 = t222 * t234;
	t272 = t223 * t222;
	t271 = t223 * t225;
	t268 = t225 * t230;
	t266 = t226 * t236;
	t263 = t232 * t223;
	t254 = t166 ^ 2 * t164 + 0.1e1;
	t253 = -t161 * t181 + t282;
	t208 = -t218 * t266 + t248 * t231;
	t250 = t208 * t225 + t218 * t272;
	t209 = -t218 * t267 - t248 * t236;
	t207 = t216 * t231 - t217 * t266;
	t200 = -t208 * t222 + t218 * t271;
	t195 = -t212 * t230 + t235 * t251;
	t194 = (t211 * t235 - t212 * t268) * t229 - t212 * t273;
	t193 = ((t229 * t268 + t273) * (-t226 * t262 - t264) + ((-t226 * t265 + t261) * t235 + t222 * t230 * t263) * t229 - t225 * t234 * t263) * t224;
	t192 = -t206 * t268 - t242;
	t191 = t225 * t277 - t230 * t243;
	t190 = t209 * t235 + t230 * t250;
	t189 = t209 * t230 - t235 * t250;
	t182 = t196 * t234 + t202 * t229;
	t179 = 0.1e1 / t181;
	t178 = t206 * t222 * t229 + t192 * t234;
	t177 = (-t204 * t235 + t205 * t268) * t229 + t205 * t273;
	t176 = t190 * t234 + t200 * t229;
	t175 = ((-t216 * t236 - t217 * t267) * t235 + (t207 * t225 + t217 * t272) * t230) * t229 - (-t207 * t222 + t217 * t271) * t234;
	t165 = 0.1e1 / (t180 * t290 ^ 2 + 0.1e1);
	t163 = 0.1e1 / t167;
	t160 = 0.1e1 / t254;
	t157 = 0.1e1 / t159;
	t156 = 0.1e1 / (0.1e1 + t279);
	t155 = (-t179 * t287 - t195 * t280) * t229 * t165;
	t154 = (-t177 * t179 - t194 * t280) * t165;
	t153 = (-t175 * t179 - t193 * t280) * t165;
	t152 = (t172 * t179 - t182 * t280) * t165;
	t1 = [-t173 * t179 * t165, t153, t154, t155, t152, 0; (t290 * t157 - (-t161 + (-t179 * t282 + t161) * t165) * t279) * t156, ((t190 * t229 - t200 * t234) * t157 - (t153 * t253 - t161 * t175 + t162 * t193) * t283) * t156, ((t192 * t229 - t206 * t273) * t157 - (t154 * t253 - t161 * t177 + t162 * t194) * t283) * t156, (-t187 * t229 * t157 - ((-t161 * t287 + t162 * t195) * t229 + t253 * t155) * t283) * t156, (t174 * t157 - (t152 * t253 + t161 * t172 + t162 * t182) * t283) * t156, 0; ((t172 * t228 - t233 * t287) * t163 - (t172 * t233 + t228 * t287) * t281) * t160, ((t176 * t228 - t189 * t233) * t163 - (t176 * t233 + t189 * t228) * t281) * t160, ((t178 * t228 - t191 * t233) * t163 - (t178 * t233 + t191 * t228) * t281) * t160, ((-t188 * t233 - t228 * t278) * t163 - (t188 * t228 - t233 * t278) * t281) * t160, (-t228 * t163 + t233 * t281) * t173 * t160, t254 * t160;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end