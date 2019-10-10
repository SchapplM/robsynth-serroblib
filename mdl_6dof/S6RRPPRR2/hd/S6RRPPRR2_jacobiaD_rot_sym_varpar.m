% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR2
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
%   Wie in S6RRPPRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:37
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:32
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (2086->84), mult. (2191->187), div. (456->12), fcn. (2616->9), ass. (0->87)
	t107 = qJ(2) + pkin(10);
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
	t111 = cos(pkin(11));
	t113 = cos(qJ(1));
	t139 = qJD(2) * t113;
	t128 = t105 * t139;
	t145 = t112 * t111;
	t110 = sin(pkin(11));
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
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:32
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (2618->95), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->93)
	t145 = qJ(2) + pkin(10);
	t141 = sin(t145);
	t136 = t141 ^ 2;
	t143 = cos(t145);
	t138 = 0.1e1 / t143 ^ 2;
	t193 = t136 * t138;
	t148 = sin(qJ(1));
	t210 = 0.2e1 * t148;
	t209 = t141 * t193;
	t146 = t148 ^ 2;
	t131 = t146 * t193 + 0.1e1;
	t129 = 0.1e1 / t131;
	t137 = 0.1e1 / t143;
	t149 = cos(qJ(1));
	t183 = qJD(1) * t149;
	t171 = t141 * t183;
	t181 = qJD(2) * t148;
	t103 = (-(-t143 * t181 - t171) * t137 + t181 * t193) * t129;
	t208 = t103 - t181;
	t144 = pkin(11) + qJ(5);
	t142 = cos(t144);
	t140 = sin(t144);
	t187 = t148 * t140;
	t188 = t143 * t149;
	t125 = t142 * t188 + t187;
	t186 = t148 * t141;
	t128 = atan2(-t186, -t143);
	t127 = cos(t128);
	t126 = sin(t128);
	t174 = t126 * t186;
	t112 = -t127 * t143 - t174;
	t109 = 0.1e1 / t112;
	t119 = 0.1e1 / t125;
	t110 = 0.1e1 / t112 ^ 2;
	t120 = 0.1e1 / t125 ^ 2;
	t207 = -0.2e1 * t141;
	t206 = t129 - 0.1e1;
	t195 = t127 * t141;
	t98 = (-t103 * t148 + qJD(2)) * t195 + (t143 * t208 - t171) * t126;
	t205 = t109 * t110 * t98;
	t159 = t142 * t149 + t143 * t187;
	t180 = qJD(2) * t149;
	t170 = t141 * t180;
	t104 = t159 * qJD(1) - qJD(5) * t125 + t140 * t170;
	t185 = t148 * t142;
	t124 = t140 * t188 - t185;
	t118 = t124 ^ 2;
	t117 = t118 * t120 + 0.1e1;
	t198 = t120 * t124;
	t164 = -qJD(1) * t143 + qJD(5);
	t165 = qJD(5) * t143 - qJD(1);
	t190 = t140 * t149;
	t105 = -t165 * t190 + (t164 * t148 - t170) * t142;
	t202 = t105 * t119 * t120;
	t204 = (-t104 * t198 - t118 * t202) / t117 ^ 2;
	t203 = t103 * t141;
	t201 = t110 * t141;
	t191 = t137 * t141;
	t158 = qJD(2) * (t137 * t209 + t191);
	t162 = t136 * t148 * t183;
	t200 = (t138 * t162 + t146 * t158) / t131 ^ 2;
	t199 = t119 * t140;
	t197 = t124 * t142;
	t196 = t126 * t148;
	t194 = t136 * t137;
	t147 = t149 ^ 2;
	t192 = t136 * t147;
	t189 = t141 * t149;
	t184 = qJD(1) * t148;
	t182 = qJD(2) * t143;
	t108 = t110 * t192 + 0.1e1;
	t179 = 0.2e1 / t108 ^ 2 * (-t192 * t205 + (t141 * t147 * t182 - t162) * t110);
	t178 = 0.2e1 * t205;
	t177 = -0.2e1 * t204;
	t176 = t124 * t202;
	t175 = t110 * t189;
	t173 = t129 * t194;
	t169 = 0.1e1 + t193;
	t168 = t141 * t179;
	t167 = t200 * t207;
	t166 = t200 * t210;
	t163 = t148 * t173;
	t161 = t169 * t149;
	t160 = t120 * t197 - t199;
	t157 = t141 * t181 + t164 * t149;
	t123 = -t143 * t185 + t190;
	t116 = t169 * t148 * t129;
	t114 = 0.1e1 / t117;
	t106 = 0.1e1 / t108;
	t102 = (t206 * t141 * t126 - t127 * t163) * t149;
	t101 = -t143 * t196 + t195 + (t126 * t143 - t127 * t186) * t116;
	t99 = -t169 * t166 + (qJD(1) * t161 + t158 * t210) * t129;
	t1 = [t137 * t149 * t167 + (qJD(2) * t161 - t184 * t191) * t129, t99, 0, 0, 0, 0; (t109 * t168 + (-t109 * t182 + (qJD(1) * t102 + t98) * t201) * t106) * t148 + (t110 * t168 * t102 + (-((t103 * t163 + t206 * t182 + t167) * t126 + (t166 * t194 - t203 + (t203 + (t207 - t209) * t181) * t129) * t127) * t175 + (-t110 * t182 + t141 * t178) * t102 + (-t109 + ((-t146 + t147) * t127 * t173 + t206 * t174) * t110) * t141 * qJD(1)) * t106) * t149, (t101 * t201 - t109 * t143) * t149 * t179 + ((-t109 * t184 + (-qJD(2) * t101 - t98) * t149 * t110) * t143 + (-t109 * t180 - (-t127 * t148 * t99 - t208 * t126 + (-qJD(2) * t126 + t103 * t196 - t127 * t183) * t116) * t175 + (t110 * t184 + t149 * t178) * t101 - ((t99 - t183) * t126 + ((-t116 * t148 + 0.1e1) * qJD(2) + (t116 - t148) * t103) * t127) * t110 * t188) * t141) * t106, 0, 0, 0, 0; 0.2e1 * (t119 * t159 + t123 * t198) * t204 + (0.2e1 * t123 * t176 - t165 * t119 * t185 + t157 * t199 + (-t165 * t124 * t187 + t123 * t104 + t105 * t159 - t157 * t197) * t120) * t114, t160 * t177 * t189 + (t160 * t143 * t180 + (-t160 * t184 + ((-qJD(5) * t119 - 0.2e1 * t176) * t142 + (-t104 * t142 + (-qJD(5) * t124 + t105) * t140) * t120) * t149) * t141) * t114, 0, 0, t177 + 0.2e1 * (-t104 * t114 * t120 + (-t114 * t202 - t120 * t204) * t124) * t124, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:32
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (3326->96), mult. (2734->202), div. (498->12), fcn. (3199->9), ass. (0->96)
	t166 = qJ(2) + pkin(10);
	t162 = sin(t166);
	t158 = t162 ^ 2;
	t163 = cos(t166);
	t160 = 0.1e1 / t163 ^ 2;
	t214 = t158 * t160;
	t169 = sin(qJ(1));
	t232 = 0.2e1 * t169;
	t231 = t162 * t214;
	t167 = t169 ^ 2;
	t152 = t167 * t214 + 0.1e1;
	t150 = 0.1e1 / t152;
	t159 = 0.1e1 / t163;
	t170 = cos(qJ(1));
	t204 = qJD(1) * t170;
	t192 = t162 * t204;
	t202 = qJD(2) * t169;
	t125 = (-(-t202 * t163 - t192) * t159 + t202 * t214) * t150;
	t230 = t125 - t202;
	t164 = pkin(11) + qJ(5) + qJ(6);
	t156 = cos(t164);
	t206 = t170 * t156;
	t155 = sin(t164);
	t210 = t169 * t155;
	t145 = t163 * t206 + t210;
	t208 = t169 * t162;
	t149 = atan2(-t208, -t163);
	t148 = cos(t149);
	t147 = sin(t149);
	t195 = t147 * t208;
	t135 = -t148 * t163 - t195;
	t132 = 0.1e1 / t135;
	t139 = 0.1e1 / t145;
	t133 = 0.1e1 / t135 ^ 2;
	t140 = 0.1e1 / t145 ^ 2;
	t229 = -0.2e1 * t162;
	t228 = t150 - 0.1e1;
	t216 = t148 * t162;
	t118 = (-t125 * t169 + qJD(2)) * t216 + (t230 * t163 - t192) * t147;
	t227 = t118 * t132 * t133;
	t165 = qJD(5) + qJD(6);
	t180 = t163 * t210 + t206;
	t201 = qJD(2) * t170;
	t191 = t162 * t201;
	t123 = qJD(1) * t180 - t145 * t165 + t155 * t191;
	t207 = t170 * t155;
	t209 = t169 * t156;
	t144 = t163 * t207 - t209;
	t138 = t144 ^ 2;
	t131 = t138 * t140 + 0.1e1;
	t219 = t140 * t144;
	t185 = -qJD(1) * t163 + t165;
	t186 = t163 * t165 - qJD(1);
	t124 = -t186 * t207 + (t169 * t185 - t191) * t156;
	t225 = t124 * t139 * t140;
	t226 = (-t123 * t219 - t138 * t225) / t131 ^ 2;
	t224 = t125 * t162;
	t223 = t133 * t162;
	t222 = t133 * t170;
	t212 = t159 * t162;
	t179 = qJD(2) * (t159 * t231 + t212);
	t183 = t158 * t169 * t204;
	t221 = (t160 * t183 + t167 * t179) / t152 ^ 2;
	t220 = t139 * t155;
	t218 = t144 * t156;
	t217 = t147 * t169;
	t215 = t158 * t159;
	t168 = t170 ^ 2;
	t213 = t158 * t168;
	t211 = t162 * t170;
	t205 = qJD(1) * t169;
	t203 = qJD(2) * t163;
	t128 = t133 * t213 + 0.1e1;
	t200 = 0.2e1 * (-t213 * t227 + (t162 * t168 * t203 - t183) * t133) / t128 ^ 2;
	t199 = 0.2e1 * t227;
	t198 = -0.2e1 * t226;
	t197 = t133 * t211;
	t196 = t144 * t225;
	t194 = t150 * t215;
	t190 = 0.1e1 + t214;
	t189 = t162 * t200;
	t188 = t221 * t229;
	t187 = t221 * t232;
	t184 = t169 * t194;
	t182 = t190 * t170;
	t181 = t218 * t140 - t220;
	t178 = t162 * t202 + t170 * t185;
	t143 = -t163 * t209 + t207;
	t137 = t190 * t169 * t150;
	t129 = 0.1e1 / t131;
	t126 = 0.1e1 / t128;
	t122 = (t147 * t162 * t228 - t148 * t184) * t170;
	t121 = -t163 * t217 + t216 + (t147 * t163 - t148 * t208) * t137;
	t119 = -t190 * t187 + (qJD(1) * t182 + t179 * t232) * t150;
	t116 = t198 + 0.2e1 * (-t123 * t140 * t129 + (-t129 * t225 - t140 * t226) * t144) * t144;
	t1 = [t170 * t159 * t188 + (qJD(2) * t182 - t205 * t212) * t150, t119, 0, 0, 0, 0; (t132 * t189 + (-t132 * t203 + (qJD(1) * t122 + t118) * t223) * t126) * t169 + (t133 * t189 * t122 + (-((t125 * t184 + t203 * t228 + t188) * t147 + (t187 * t215 - t224 + (t224 + (t229 - t231) * t202) * t150) * t148) * t197 + (-t133 * t203 + t162 * t199) * t122 + (-t132 + ((-t167 + t168) * t148 * t194 + t228 * t195) * t133) * t162 * qJD(1)) * t126) * t170, (t121 * t223 - t132 * t163) * t170 * t200 + ((-t132 * t205 + (-qJD(2) * t121 - t118) * t222) * t163 + (-t132 * t201 - (-t119 * t148 * t169 - t230 * t147 + (-qJD(2) * t147 + t125 * t217 - t148 * t204) * t137) * t197 + (t133 * t205 + t170 * t199) * t121 - ((t119 - t204) * t147 + ((-t137 * t169 + 0.1e1) * qJD(2) + (t137 - t169) * t125) * t148) * t163 * t222) * t162) * t126, 0, 0, 0, 0; 0.2e1 * (t139 * t180 + t143 * t219) * t226 + (0.2e1 * t143 * t196 - t186 * t139 * t209 + t178 * t220 + (-t144 * t186 * t210 + t143 * t123 + t124 * t180 - t178 * t218) * t140) * t129, t181 * t198 * t211 + (t181 * t163 * t201 + (-t181 * t205 + ((-t139 * t165 - 0.2e1 * t196) * t156 + (-t123 * t156 + (-t144 * t165 + t124) * t155) * t140) * t170) * t162) * t129, 0, 0, t116, t116;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end