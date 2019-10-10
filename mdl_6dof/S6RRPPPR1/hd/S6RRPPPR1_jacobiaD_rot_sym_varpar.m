% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR1
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
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRPPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:09
	% EndTime: 2019-10-10 09:18:10
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (2086->84), mult. (2191->187), div. (456->12), fcn. (2616->9), ass. (0->87)
	t107 = qJ(2) + pkin(9);
	t105 = sin(t107);
	t101 = t105 ^ 2;
	t106 = cos(t107);
	t150 = t101 / t106 ^ 2;
	t112 = sin(qJ(1));
	t126 = 0.1e1 + t150;
	t108 = t112 ^ 2;
	t99 = t108 * t150 + 0.1e1;
	t97 = 0.1e1 / t99;
	t123 = t126 * t97;
	t80 = t112 * t123;
	t171 = t112 * t80 - 0.1e1;
	t111 = cos(pkin(10));
	t113 = cos(qJ(1));
	t139 = qJD(2) * t113;
	t128 = t105 * t139;
	t145 = t112 * t111;
	t110 = sin(pkin(10));
	t149 = t110 * t113;
	t91 = -t106 * t145 + t149;
	t85 = t91 * qJD(1) - t111 * t128;
	t146 = t112 * t110;
	t148 = t111 * t113;
	t93 = t106 * t148 + t146;
	t88 = 0.1e1 / t93 ^ 2;
	t170 = t85 * t88;
	t92 = t106 * t149 - t145;
	t157 = t88 * t92;
	t86 = t92 ^ 2;
	t83 = t86 * t88 + 0.1e1;
	t81 = 0.1e1 / t83;
	t87 = 0.1e1 / t93;
	t169 = (-t110 * t87 + t111 * t157) * t81;
	t168 = t105 * t150;
	t147 = t112 * t105;
	t96 = atan2(-t147, -t106);
	t94 = sin(t96);
	t133 = t94 * t147;
	t95 = cos(t96);
	t78 = -t106 * t95 - t133;
	t75 = 0.1e1 / t78;
	t102 = 0.1e1 / t106;
	t76 = 0.1e1 / t78 ^ 2;
	t167 = 0.2e1 * t105;
	t166 = t97 - 0.1e1;
	t109 = t113 ^ 2;
	t142 = qJD(1) * t113;
	t130 = t112 * t142;
	t141 = qJD(2) * t106;
	t131 = t76 * t141;
	t140 = qJD(2) * t112;
	t154 = t106 * t94;
	t71 = (-(-t105 * t142 - t106 * t140) * t102 + t140 * t150) * t97;
	t66 = (t71 - t140) * t154 + (-t94 * t142 + (-t112 * t71 + qJD(2)) * t95) * t105;
	t164 = t66 * t75 * t76;
	t156 = t101 * t76;
	t74 = t109 * t156 + 0.1e1;
	t165 = (t109 * t105 * t131 + (-t109 * t164 - t76 * t130) * t101) / t74 ^ 2;
	t158 = t87 * t170;
	t90 = -t106 * t146 - t148;
	t84 = t90 * qJD(1) - t110 * t128;
	t163 = (t84 * t157 - t86 * t158) / t83 ^ 2;
	t72 = 0.1e1 / t74;
	t161 = t72 * t76;
	t160 = t75 * t72;
	t121 = qJD(2) * (t105 + t168) * t102;
	t159 = (t108 * t121 + t130 * t150) / t99 ^ 2;
	t155 = t102 * t97;
	t152 = t113 * t76;
	t151 = qJD(2) * t80;
	t144 = qJD(1) * t105;
	t143 = qJD(1) * t112;
	t138 = 0.2e1 * t164;
	t137 = 0.2e1 * t163;
	t136 = t75 * t165;
	t135 = t92 * t158;
	t134 = t112 * t155;
	t132 = t166 * t105;
	t129 = t105 * t140;
	t127 = 0.2e1 * t76 * t165;
	t125 = -0.2e1 * t102 * t159;
	t124 = t101 * t134;
	t70 = (-t95 * t124 + t94 * t132) * t113;
	t68 = (-t112 + t80) * t154 - t171 * t95 * t105;
	t67 = t123 * t142 + 0.2e1 * (t121 * t97 - t126 * t159) * t112;
	t1 = [-t134 * t144 + (qJD(2) * t123 + t105 * t125) * t113, t67, 0, 0, 0, 0; (-t141 * t160 + (0.2e1 * t136 + (qJD(1) * t70 + t66) * t161) * t105) * t112 + (t70 * t127 * t105 + (-t70 * t131 + (t70 * t138 + ((-t71 * t124 - t166 * t141 + t159 * t167) * t94 + (-t71 * t132 + (t101 * t125 + (t167 + t168) * t97 * qJD(2)) * t112) * t95) * t152) * t105 + (-t75 + t166 * t76 * t133 - (t108 - t109) * t95 * t155 * t156) * t144) * t72) * t113, (-t143 * t160 + (-0.2e1 * t136 + (-qJD(2) * t68 - t66) * t161) * t113) * t106 + (t68 * t113 * t127 + (-t75 * t139 + (t113 * t138 + t76 * t143) * t68 + (-((-t112 * t67 - t142 * t80) * t95 + (t171 * t71 + t140 - t151) * t94) * t105 - ((t67 - t142) * t94 + (t71 * t80 + qJD(2) + (-t71 - t151) * t112) * t95) * t106) * t152) * t72) * t105, 0, 0, 0, 0; (t91 * t157 - t87 * t90) * t137 + ((-t92 * qJD(1) + t110 * t129) * t87 + 0.2e1 * t91 * t135 + (-t90 * t85 - (-t93 * qJD(1) + t111 * t129) * t92 - t91 * t84) * t88) * t81, t106 * t139 * t169 + (-t143 * t169 + ((t87 * t137 + t81 * t170) * t110 + (-0.2e1 * t157 * t163 + (t84 * t88 - 0.2e1 * t135) * t81) * t111) * t113) * t105, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:09
	% EndTime: 2019-10-10 09:18:09
	% DurationCPUTime: 0.86s
	% Computational Cost: add. (2051->79), mult. (3228->188), div. (613->15), fcn. (4191->9), ass. (0->85)
	t133 = cos(pkin(10));
	t181 = 0.2e1 * t133;
	t130 = qJ(2) + pkin(9);
	t127 = cos(t130);
	t132 = sin(pkin(10));
	t134 = sin(qJ(1));
	t162 = t134 * t132;
	t135 = cos(qJ(1));
	t163 = t133 * t135;
	t114 = t127 * t162 + t163;
	t126 = sin(t130);
	t165 = t126 * t132;
	t104 = atan2(-t114, t165);
	t101 = cos(t104);
	t100 = sin(t104);
	t171 = t100 * t114;
	t96 = t101 * t165 - t171;
	t93 = 0.1e1 / t96;
	t109 = t114 ^ 2;
	t124 = 0.1e1 / t126 ^ 2;
	t129 = 0.1e1 / t132 ^ 2;
	t107 = t109 * t124 * t129 + 0.1e1;
	t102 = 0.1e1 / t107;
	t164 = t132 * t135;
	t150 = t127 * t164;
	t161 = t134 * t133;
	t117 = t150 - t161;
	t123 = 0.1e1 / t126;
	t128 = 0.1e1 / t132;
	t169 = t114 * t128;
	t151 = t123 * t169;
	t85 = (-t100 + (t101 * t151 + t100) * t102) * t117;
	t180 = 0.2e1 * t85;
	t118 = t127 * t163 + t162;
	t111 = 0.1e1 / t118;
	t112 = 0.1e1 / t118 ^ 2;
	t94 = 0.1e1 / t96 ^ 2;
	t110 = t117 ^ 2;
	t156 = qJD(2) * t135;
	t147 = t126 * t156;
	t97 = t114 * qJD(1) + t132 * t147;
	t176 = t97 * t94;
	t158 = qJD(2) * t127;
	t149 = t124 * t158;
	t157 = qJD(2) * t134;
	t160 = qJD(1) * t134;
	t99 = qJD(1) * t150 - t133 * t160 - t157 * t165;
	t141 = t114 * t149 - t123 * t99;
	t86 = t141 * t128 * t102;
	t83 = (-t114 * t86 + t132 * t158) * t101 + (-t86 * t165 - t99) * t100;
	t178 = t83 * t93 * t94;
	t90 = t110 * t94 + 0.1e1;
	t179 = (-t110 * t178 - t117 * t176) / t90 ^ 2;
	t88 = 0.1e1 / t90;
	t177 = t88 * t94;
	t122 = t126 ^ 2;
	t125 = t123 / t122;
	t175 = 0.1e1 / t107 ^ 2 * (-t109 * t125 * t158 + t114 * t124 * t99) * t129;
	t116 = -t127 * t161 + t164;
	t143 = t133 * t147;
	t98 = t116 * qJD(1) - t143;
	t174 = t111 * t112 * t98;
	t173 = t123 * t97;
	t166 = t124 * t127;
	t142 = t166 * t169 + t134;
	t92 = t142 * t102;
	t172 = t134 - t92;
	t170 = t100 * t117;
	t168 = t116 * t135;
	t131 = t135 ^ 2;
	t167 = t122 * t131;
	t159 = qJD(1) * t135;
	t155 = -0.2e1 * t175;
	t152 = t112 * t167;
	t108 = 0.1e1 + t152;
	t144 = t122 * t134 * t159;
	t145 = t167 * t174;
	t148 = qJD(2) * t126 * t131;
	t154 = 0.2e1 / t108 ^ 2 * (-t145 + (t127 * t148 - t144) * t112);
	t153 = t94 * t179;
	t146 = 0.2e1 * t93 * t179;
	t105 = 0.1e1 / t108;
	t84 = -t101 * t114 * t92 + (t172 * t126 * t100 + t101 * t127) * t132;
	t82 = t142 * t155 + (t159 + (t99 * t166 + (-0.2e1 * t125 * t127 ^ 2 - t123) * t114 * qJD(2)) * t128) * t102;
	t1 = [(0.2e1 * t117 * t123 * t175 + (t117 * t149 + t173) * t102) * t128, t82, 0, 0, 0, 0; t114 * t146 + (-t99 * t93 + (t114 * t83 + t85 * t97) * t94) * t88 + (t153 * t180 + (t178 * t180 + ((t102 * t97 - t97 - (-t102 * t86 * t151 + t155) * t117) * t100 + (-(t151 * t155 - t86) * t117 + (-t117 * t86 + (t114 * t173 + t141 * t117) * t128) * t102) * t101) * t94) * t88) * t117, t84 * t88 * t176 + (-(t92 * t86 * t171 + (-t114 * t82 - t92 * t99) * t101) * t177 + 0.2e1 * (t88 * t178 + t153) * t84) * t117 + ((-t93 * t88 * t156 - (t172 * qJD(2) - t86) * t170 * t177) * t127 + (t135 * t146 + (t93 * t160 + (t135 * t83 - (-t82 + t159) * t170 - (t172 * t86 - qJD(2)) * t117 * t101) * t94) * t88) * t126) * t132, 0, 0, 0, 0; (-t111 * t134 - t112 * t168) * t126 * t154 + (-0.2e1 * t126 * t168 * t174 + (t126 * t159 + t127 * t157) * t111 + (((-t98 + t143) * t134 - t118 * t159) * t126 + (-t126 * t160 + t127 * t156) * t116) * t112) * t105, (t111 * t127 * t135 + t133 * t152) * t154 + (t145 * t181 + (t127 * t160 + t147) * t111 + (t144 * t181 + (-0.2e1 * t133 * t148 + t135 * t98) * t127) * t112) * t105, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:09
	% EndTime: 2019-10-10 09:18:10
	% DurationCPUTime: 1.21s
	% Computational Cost: add. (2357->105), mult. (3345->232), div. (468->12), fcn. (4023->11), ass. (0->107)
	t168 = qJ(2) + pkin(9);
	t166 = sin(t168);
	t162 = t166 ^ 2;
	t167 = cos(t168);
	t164 = 0.1e1 / t167 ^ 2;
	t223 = t162 * t164;
	t174 = sin(qJ(1));
	t169 = t174 ^ 2;
	t160 = t169 * t223 + 0.1e1;
	t163 = 0.1e1 / t167;
	t220 = t163 * t166;
	t242 = t166 * t223;
	t184 = qJD(2) * (t163 * t242 + t220);
	t176 = cos(qJ(1));
	t212 = qJD(1) * t176;
	t221 = t162 * t174;
	t190 = t212 * t221;
	t228 = (t164 * t190 + t169 * t184) / t160 ^ 2;
	t243 = -0.2e1 * t228;
	t194 = 0.1e1 + t223;
	t241 = t174 * t194;
	t157 = 0.1e1 / t160;
	t199 = t166 * t212;
	t210 = qJD(2) * t174;
	t126 = ((t167 * t210 + t199) * t163 + t210 * t223) * t157;
	t240 = t126 - t210;
	t173 = sin(qJ(6));
	t175 = cos(qJ(6));
	t208 = qJD(6) * t176;
	t213 = qJD(1) * t174;
	t239 = t173 * t213 - t175 * t208;
	t238 = t173 * t208 + t175 * t213;
	t216 = t174 * t166;
	t159 = atan2(t216, t167);
	t156 = cos(t159);
	t155 = sin(t159);
	t203 = t155 * t216;
	t133 = t156 * t167 + t203;
	t130 = 0.1e1 / t133;
	t172 = cos(pkin(10));
	t214 = t174 * t172;
	t171 = sin(pkin(10));
	t218 = t171 * t176;
	t153 = t167 * t218 - t214;
	t215 = t174 * t171;
	t217 = t172 * t176;
	t154 = t167 * t217 + t215;
	t143 = t153 * t173 + t154 * t175;
	t137 = 0.1e1 / t143;
	t131 = 0.1e1 / t133 ^ 2;
	t138 = 0.1e1 / t143 ^ 2;
	t237 = 0.2e1 * t166;
	t236 = t157 - 0.1e1;
	t170 = t176 ^ 2;
	t222 = t162 * t170;
	t129 = t131 * t222 + 0.1e1;
	t211 = qJD(2) * t167;
	t224 = t156 * t166;
	t117 = (t126 * t174 - qJD(2)) * t224 + (-t240 * t167 + t199) * t155;
	t234 = t117 * t130 * t131;
	t235 = (-t222 * t234 + (t166 * t170 * t211 - t190) * t131) / t129 ^ 2;
	t151 = -t167 * t215 - t217;
	t209 = qJD(2) * t176;
	t197 = t166 * t209;
	t144 = t151 * qJD(1) - t171 * t197;
	t152 = -t167 * t214 + t218;
	t145 = t152 * qJD(1) - t172 * t197;
	t187 = t153 * t175 - t154 * t173;
	t121 = t187 * qJD(6) + t144 * t173 + t145 * t175;
	t139 = t137 * t138;
	t233 = t121 * t139;
	t120 = t143 * qJD(6) - t144 * t175 + t145 * t173;
	t136 = t187 ^ 2;
	t125 = t136 * t138 + 0.1e1;
	t227 = t138 * t187;
	t232 = 0.1e1 / t125 ^ 2 * (-t120 * t227 - t136 * t233);
	t231 = t126 * t166;
	t230 = t131 * t166;
	t229 = t131 * t176;
	t185 = -t171 * t173 - t172 * t175;
	t219 = t166 * t176;
	t149 = t185 * t219;
	t226 = t138 * t149;
	t225 = t155 * t174;
	t207 = -0.2e1 * t234;
	t206 = 0.2e1 * t232;
	t205 = -0.2e1 * t139 * t187;
	t204 = t131 * t219;
	t202 = t157 * t162 * t163;
	t198 = t166 * t210;
	t193 = -0.2e1 * t166 * t235;
	t192 = t163 * t243;
	t191 = t174 * t202;
	t189 = t194 * t176;
	t188 = t151 * t175 - t152 * t173;
	t141 = t151 * t173 + t152 * t175;
	t186 = t171 * t175 - t172 * t173;
	t148 = t186 * t219;
	t147 = -t154 * qJD(1) + t172 * t198;
	t146 = -t153 * qJD(1) + t171 * t198;
	t135 = t157 * t241;
	t127 = 0.1e1 / t129;
	t123 = 0.1e1 / t125;
	t122 = (-t236 * t166 * t155 + t156 * t191) * t176;
	t119 = t167 * t225 - t224 + (-t155 * t167 + t156 * t216) * t135;
	t118 = t241 * t243 + (qJD(1) * t189 + 0.2e1 * t174 * t184) * t157;
	t1 = [t192 * t219 + (qJD(2) * t189 - t213 * t220) * t157, t118, 0, 0, 0, 0; (t130 * t193 + (t130 * t211 + (-qJD(1) * t122 - t117) * t230) * t127) * t174 + (t131 * t193 * t122 + (((-t126 * t191 - t236 * t211 + t228 * t237) * t155 + (t192 * t221 + t231 + (-t231 + (t237 + t242) * t210) * t157) * t156) * t204 + (t131 * t211 + t166 * t207) * t122 + (t130 + ((-t169 + t170) * t156 * t202 + t236 * t203) * t131) * t166 * qJD(1)) * t127) * t176, 0.2e1 * (-t119 * t230 + t130 * t167) * t176 * t235 + ((t130 * t213 + (qJD(2) * t119 + t117) * t229) * t167 + (t130 * t209 + (t118 * t156 * t174 + t240 * t155 + (qJD(2) * t155 - t126 * t225 + t156 * t212) * t135) * t204 + (-t131 * t213 + t176 * t207) * t119 + ((-t118 + t212) * t155 + ((t135 * t174 - 0.1e1) * qJD(2) + (-t135 + t174) * t126) * t156) * t167 * t229) * t166) * t127, 0, 0, 0, 0; (t137 * t188 - t141 * t227) * t206 + ((t141 * qJD(6) - t146 * t175 + t147 * t173) * t137 + t141 * t121 * t205 + (t188 * t121 + (t188 * qJD(6) + t146 * t173 + t147 * t175) * t187 - t141 * t120) * t138) * t123, (-t137 * t148 - t187 * t226) * t206 + (-t120 * t226 + (-t138 * t148 + t149 * t205) * t121 + (t186 * t137 + t185 * t227) * t167 * t209 + ((t239 * t137 + t238 * t227) * t172 + (-t238 * t137 + t239 * t227) * t171) * t166) * t123, 0, 0, 0, -0.2e1 * t232 - 0.2e1 * (t120 * t138 * t123 - (-t123 * t233 - t138 * t232) * t187) * t187;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end