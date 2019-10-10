% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR3
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
%   Wie in S6RRPRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:54
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
	t139 = sin(qJ(1));
	t136 = t139 ^ 2;
	t135 = qJ(2) + pkin(10);
	t133 = sin(t135);
	t129 = t133 ^ 2;
	t134 = cos(t135);
	t131 = 0.1e1 / t134 ^ 2;
	t188 = t129 * t131;
	t124 = t136 * t188 + 0.1e1;
	t128 = t133 * t129;
	t130 = 0.1e1 / t134;
	t185 = t130 * t133;
	t149 = qJD(2) * (t128 * t130 * t131 + t185);
	t141 = cos(qJ(1));
	t177 = qJD(1) * t141;
	t186 = t129 * t139;
	t154 = t177 * t186;
	t194 = (t131 * t154 + t136 * t149) / t124 ^ 2;
	t204 = -0.2e1 * t194;
	t161 = 0.1e1 + t188;
	t203 = t139 * t161;
	t140 = cos(qJ(4));
	t179 = t141 * t140;
	t138 = sin(qJ(4));
	t182 = t139 * t138;
	t120 = t134 * t179 + t182;
	t183 = t139 * t133;
	t121 = atan2(-t183, -t134);
	t116 = cos(t121);
	t115 = sin(t121);
	t168 = t115 * t183;
	t105 = -t116 * t134 - t168;
	t102 = 0.1e1 / t105;
	t112 = 0.1e1 / t120;
	t103 = 0.1e1 / t105 ^ 2;
	t113 = 0.1e1 / t120 ^ 2;
	t122 = 0.1e1 / t124;
	t202 = t122 - 0.1e1;
	t137 = t141 ^ 2;
	t176 = qJD(2) * t134;
	t187 = t129 * t137;
	t175 = qJD(2) * t139;
	t163 = t131 * t175;
	t164 = t133 * t177;
	t96 = (-(-t134 * t175 - t164) * t130 + t129 * t163) * t122;
	t158 = t96 - t175;
	t159 = -t139 * t96 + qJD(2);
	t190 = t116 * t133;
	t91 = t159 * t190 + (t158 * t134 - t164) * t115;
	t199 = t102 * t103 * t91;
	t99 = t103 * t187 + 0.1e1;
	t201 = (-t187 * t199 + (t133 * t137 * t176 - t154) * t103) / t99 ^ 2;
	t97 = 0.1e1 / t99;
	t200 = t103 * t97;
	t150 = t134 * t182 + t179;
	t174 = qJD(2) * t141;
	t162 = t133 * t174;
	t100 = t150 * qJD(1) - t120 * qJD(4) + t138 * t162;
	t180 = t141 * t138;
	t181 = t139 * t140;
	t119 = t134 * t180 - t181;
	t111 = t119 ^ 2;
	t110 = t111 * t113 + 0.1e1;
	t192 = t113 * t119;
	t156 = -qJD(1) * t134 + qJD(4);
	t157 = qJD(4) * t134 - qJD(1);
	t101 = -t157 * t180 + (t156 * t139 - t162) * t140;
	t196 = t101 * t112 * t113;
	t198 = 0.1e1 / t110 ^ 2 * (-t100 * t192 - t111 * t196);
	t193 = t112 * t138;
	t191 = t115 * t134;
	t189 = t119 * t140;
	t184 = t133 * t141;
	t178 = qJD(1) * t139;
	t173 = 0.2e1 * t201;
	t172 = 0.2e1 * t199;
	t171 = -0.2e1 * t198;
	t170 = t102 * t201;
	t169 = t97 * t176;
	t167 = t119 * t196;
	t166 = t122 * t129 * t130;
	t160 = t130 * t204;
	t155 = t139 * t166;
	t153 = t161 * t141;
	t152 = t156 * t141;
	t151 = t113 * t189 - t193;
	t118 = -t134 * t181 + t180;
	t108 = 0.1e1 / t110;
	t107 = t122 * t203;
	t95 = (t202 * t133 * t115 - t116 * t155) * t141;
	t93 = -t139 * t191 + t190 + (-t116 * t183 + t191) * t107;
	t92 = t203 * t204 + (qJD(1) * t153 + 0.2e1 * t139 * t149) * t122;
	t1 = [t160 * t184 + (qJD(2) * t153 - t178 * t185) * t122, t92, 0, 0, 0, 0; (-t102 * t169 + (0.2e1 * t170 + (qJD(1) * t95 + t91) * t200) * t133) * t139 + ((-t95 * t169 + (t95 * t173 + ((0.2e1 * t133 * t194 - t96 * t155 - t202 * t176) * t115 + (t160 * t186 + t133 * t96 + (t128 * t163 - (t96 - 0.2e1 * t175) * t133) * t122) * t116) * t97 * t141) * t133) * t103 + (t95 * t172 + (-t102 + ((-t136 + t137) * t116 * t166 + t202 * t168) * t103) * qJD(1)) * t133 * t97) * t141, (-t102 * t97 * t178 + (-0.2e1 * t170 + (-qJD(2) * t93 - t91) * t200) * t141) * t134 + (((-qJD(2) * t102 + t93 * t172) * t141 + (t93 * t178 + (-(-t107 * t177 - t139 * t92) * t116 - ((t107 * t139 - 0.1e1) * t96 + (-t107 + t139) * qJD(2)) * t115) * t184) * t103) * t97 + (t93 * t173 - ((t92 - t177) * t115 + (t158 * t107 + t159) * t116) * t97 * t134) * t103 * t141) * t133, 0, 0, 0, 0; 0.2e1 * (t112 * t150 + t118 * t192) * t198 + (0.2e1 * t118 * t167 - t157 * t112 * t181 + (t133 * t175 + t152) * t193 + (t118 * t100 + t150 * t101 - t152 * t189 - (qJD(2) * t133 * t140 + t157 * t138) * t119 * t139) * t113) * t108, t151 * t171 * t184 + (t151 * t134 * t174 + (-t151 * t178 + ((-qJD(4) * t112 - 0.2e1 * t167) * t140 + (-t100 * t140 + (-qJD(4) * t119 + t101) * t138) * t113) * t141) * t133) * t108, 0, t171 + 0.2e1 * (-t100 * t113 * t108 + (-t108 * t196 - t113 * t198) * t119) * t119, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:54
	% DurationCPUTime: 1.19s
	% Computational Cost: add. (2618->95), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->93)
	t147 = qJ(2) + pkin(10);
	t143 = sin(t147);
	t138 = t143 ^ 2;
	t145 = cos(t147);
	t140 = 0.1e1 / t145 ^ 2;
	t195 = t138 * t140;
	t150 = sin(qJ(1));
	t212 = 0.2e1 * t150;
	t211 = t143 * t195;
	t148 = t150 ^ 2;
	t133 = t148 * t195 + 0.1e1;
	t131 = 0.1e1 / t133;
	t139 = 0.1e1 / t145;
	t151 = cos(qJ(1));
	t185 = qJD(1) * t151;
	t173 = t143 * t185;
	t183 = qJD(2) * t150;
	t105 = (-(-t145 * t183 - t173) * t139 + t183 * t195) * t131;
	t210 = t105 - t183;
	t146 = qJ(4) + pkin(11);
	t144 = cos(t146);
	t142 = sin(t146);
	t189 = t150 * t142;
	t190 = t145 * t151;
	t127 = t144 * t190 + t189;
	t188 = t150 * t143;
	t130 = atan2(-t188, -t145);
	t129 = cos(t130);
	t128 = sin(t130);
	t176 = t128 * t188;
	t114 = -t129 * t145 - t176;
	t111 = 0.1e1 / t114;
	t121 = 0.1e1 / t127;
	t112 = 0.1e1 / t114 ^ 2;
	t122 = 0.1e1 / t127 ^ 2;
	t209 = -0.2e1 * t143;
	t208 = t131 - 0.1e1;
	t197 = t129 * t143;
	t100 = (-t105 * t150 + qJD(2)) * t197 + (t210 * t145 - t173) * t128;
	t207 = t100 * t111 * t112;
	t161 = t144 * t151 + t145 * t189;
	t182 = qJD(2) * t151;
	t172 = t143 * t182;
	t106 = t161 * qJD(1) - t127 * qJD(4) + t142 * t172;
	t187 = t150 * t144;
	t126 = t142 * t190 - t187;
	t120 = t126 ^ 2;
	t119 = t120 * t122 + 0.1e1;
	t200 = t122 * t126;
	t166 = -qJD(1) * t145 + qJD(4);
	t167 = qJD(4) * t145 - qJD(1);
	t192 = t142 * t151;
	t107 = -t167 * t192 + (t166 * t150 - t172) * t144;
	t204 = t107 * t121 * t122;
	t206 = (-t106 * t200 - t120 * t204) / t119 ^ 2;
	t205 = t105 * t143;
	t203 = t112 * t143;
	t193 = t139 * t143;
	t160 = qJD(2) * (t139 * t211 + t193);
	t164 = t138 * t150 * t185;
	t202 = (t140 * t164 + t148 * t160) / t133 ^ 2;
	t201 = t121 * t142;
	t199 = t126 * t144;
	t198 = t128 * t150;
	t196 = t138 * t139;
	t149 = t151 ^ 2;
	t194 = t138 * t149;
	t191 = t143 * t151;
	t186 = qJD(1) * t150;
	t184 = qJD(2) * t145;
	t110 = t112 * t194 + 0.1e1;
	t181 = 0.2e1 / t110 ^ 2 * (-t194 * t207 + (t143 * t149 * t184 - t164) * t112);
	t180 = 0.2e1 * t207;
	t179 = -0.2e1 * t206;
	t178 = t126 * t204;
	t177 = t112 * t191;
	t175 = t131 * t196;
	t171 = 0.1e1 + t195;
	t170 = t143 * t181;
	t169 = t202 * t209;
	t168 = t202 * t212;
	t165 = t150 * t175;
	t163 = t171 * t151;
	t162 = t122 * t199 - t201;
	t159 = t143 * t183 + t166 * t151;
	t125 = -t145 * t187 + t192;
	t118 = t171 * t150 * t131;
	t116 = 0.1e1 / t119;
	t108 = 0.1e1 / t110;
	t104 = (t208 * t143 * t128 - t129 * t165) * t151;
	t103 = -t145 * t198 + t197 + (t128 * t145 - t129 * t188) * t118;
	t101 = -t171 * t168 + (qJD(1) * t163 + t160 * t212) * t131;
	t1 = [t139 * t151 * t169 + (qJD(2) * t163 - t186 * t193) * t131, t101, 0, 0, 0, 0; (t111 * t170 + (-t111 * t184 + (qJD(1) * t104 + t100) * t203) * t108) * t150 + (t112 * t170 * t104 + (-((t105 * t165 + t208 * t184 + t169) * t128 + (t168 * t196 - t205 + (t205 + (t209 - t211) * t183) * t131) * t129) * t177 + (-t112 * t184 + t143 * t180) * t104 + (-t111 + ((-t148 + t149) * t129 * t175 + t208 * t176) * t112) * t143 * qJD(1)) * t108) * t151, (t103 * t203 - t111 * t145) * t151 * t181 + ((-t111 * t186 + (-qJD(2) * t103 - t100) * t151 * t112) * t145 + (-t111 * t182 - (-t101 * t129 * t150 - t210 * t128 + (-qJD(2) * t128 + t105 * t198 - t129 * t185) * t118) * t177 + (t112 * t186 + t151 * t180) * t103 - ((t101 - t185) * t128 + ((-t118 * t150 + 0.1e1) * qJD(2) + (t118 - t150) * t105) * t129) * t112 * t190) * t143) * t108, 0, 0, 0, 0; 0.2e1 * (t121 * t161 + t125 * t200) * t206 + (0.2e1 * t125 * t178 - t167 * t121 * t187 + t159 * t201 + (-t167 * t126 * t189 + t125 * t106 + t107 * t161 - t159 * t199) * t122) * t116, t162 * t179 * t191 + (t162 * t145 * t182 + (-t162 * t186 + ((-qJD(4) * t121 - 0.2e1 * t178) * t144 + (-t106 * t144 + (-qJD(4) * t126 + t107) * t142) * t122) * t151) * t143) * t116, 0, t179 + 0.2e1 * (-t106 * t116 * t122 + (-t116 * t204 - t122 * t206) * t126) * t126, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:54
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
	t164 = qJ(4) + pkin(11) + qJ(6);
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
	t165 = qJD(4) + qJD(6);
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
	t1 = [t170 * t159 * t188 + (qJD(2) * t182 - t205 * t212) * t150, t119, 0, 0, 0, 0; (t132 * t189 + (-t132 * t203 + (qJD(1) * t122 + t118) * t223) * t126) * t169 + (t133 * t189 * t122 + (-((t125 * t184 + t228 * t203 + t188) * t147 + (t187 * t215 - t224 + (t224 + (t229 - t231) * t202) * t150) * t148) * t197 + (-t133 * t203 + t162 * t199) * t122 + (-t132 + ((-t167 + t168) * t148 * t194 + t228 * t195) * t133) * t162 * qJD(1)) * t126) * t170, (t121 * t223 - t132 * t163) * t170 * t200 + ((-t132 * t205 + (-qJD(2) * t121 - t118) * t222) * t163 + (-t132 * t201 - (-t119 * t148 * t169 - t230 * t147 + (-qJD(2) * t147 + t125 * t217 - t148 * t204) * t137) * t197 + (t133 * t205 + t170 * t199) * t121 - ((t119 - t204) * t147 + ((-t137 * t169 + 0.1e1) * qJD(2) + (t137 - t169) * t125) * t148) * t163 * t222) * t162) * t126, 0, 0, 0, 0; 0.2e1 * (t139 * t180 + t143 * t219) * t226 + (0.2e1 * t143 * t196 - t186 * t139 * t209 + t178 * t220 + (-t144 * t186 * t210 + t143 * t123 + t124 * t180 - t178 * t218) * t140) * t129, t181 * t198 * t211 + (t181 * t163 * t201 + (-t181 * t205 + ((-t139 * t165 - 0.2e1 * t196) * t156 + (-t123 * t156 + (-t144 * t165 + t124) * t155) * t140) * t170) * t162) * t129, 0, t116, 0, t116;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end