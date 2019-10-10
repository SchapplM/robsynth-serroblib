% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP4
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
%   Wie in S6RPRRPP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:56
	% EndTime: 2019-10-10 01:14:57
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
	t136 = sin(qJ(1));
	t133 = t136 ^ 2;
	t132 = pkin(9) + qJ(3);
	t130 = sin(t132);
	t126 = t130 ^ 2;
	t131 = cos(t132);
	t128 = 0.1e1 / t131 ^ 2;
	t185 = t126 * t128;
	t121 = t133 * t185 + 0.1e1;
	t125 = t130 * t126;
	t127 = 0.1e1 / t131;
	t182 = t127 * t130;
	t146 = qJD(3) * (t125 * t127 * t128 + t182);
	t138 = cos(qJ(1));
	t174 = qJD(1) * t138;
	t183 = t126 * t136;
	t151 = t174 * t183;
	t191 = (t128 * t151 + t133 * t146) / t121 ^ 2;
	t201 = -0.2e1 * t191;
	t158 = 0.1e1 + t185;
	t200 = t136 * t158;
	t137 = cos(qJ(4));
	t176 = t137 * t138;
	t135 = sin(qJ(4));
	t178 = t136 * t135;
	t117 = t131 * t176 + t178;
	t179 = t136 * t130;
	t118 = atan2(-t179, -t131);
	t113 = cos(t118);
	t112 = sin(t118);
	t164 = t112 * t179;
	t102 = -t113 * t131 - t164;
	t99 = 0.1e1 / t102;
	t109 = 0.1e1 / t117;
	t100 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t117 ^ 2;
	t119 = 0.1e1 / t121;
	t199 = t119 - 0.1e1;
	t134 = t138 ^ 2;
	t173 = qJD(3) * t131;
	t184 = t126 * t134;
	t172 = qJD(3) * t136;
	t160 = t128 * t172;
	t161 = t130 * t174;
	t93 = (-(-t131 * t172 - t161) * t127 + t126 * t160) * t119;
	t155 = t93 - t172;
	t156 = -t136 * t93 + qJD(3);
	t187 = t113 * t130;
	t88 = t156 * t187 + (t155 * t131 - t161) * t112;
	t196 = t99 * t100 * t88;
	t96 = t100 * t184 + 0.1e1;
	t198 = (-t184 * t196 + (t130 * t134 * t173 - t151) * t100) / t96 ^ 2;
	t94 = 0.1e1 / t96;
	t197 = t100 * t94;
	t177 = t136 * t137;
	t180 = t135 * t138;
	t116 = t131 * t180 - t177;
	t108 = t116 ^ 2;
	t107 = t108 * t110 + 0.1e1;
	t189 = t110 * t116;
	t153 = -qJD(1) * t131 + qJD(4);
	t154 = qJD(4) * t131 - qJD(1);
	t171 = qJD(3) * t138;
	t159 = t130 * t171;
	t98 = -t154 * t180 + (t153 * t136 - t159) * t137;
	t194 = t109 * t110 * t98;
	t147 = t131 * t178 + t176;
	t97 = t147 * qJD(1) - t117 * qJD(4) + t135 * t159;
	t195 = 0.1e1 / t107 ^ 2 * (-t108 * t194 - t97 * t189);
	t190 = t109 * t135;
	t188 = t112 * t131;
	t186 = t116 * t137;
	t181 = t130 * t138;
	t175 = qJD(1) * t136;
	t170 = 0.2e1 * t198;
	t169 = 0.2e1 * t196;
	t168 = -0.2e1 * t195;
	t167 = t99 * t198;
	t166 = t116 * t194;
	t165 = t94 * t173;
	t163 = t119 * t126 * t127;
	t157 = t127 * t201;
	t152 = t136 * t163;
	t150 = t158 * t138;
	t149 = t153 * t138;
	t148 = t110 * t186 - t190;
	t115 = -t131 * t177 + t180;
	t105 = 0.1e1 / t107;
	t104 = t119 * t200;
	t92 = (t199 * t130 * t112 - t113 * t152) * t138;
	t90 = -t136 * t188 + t187 + (-t113 * t179 + t188) * t104;
	t89 = t200 * t201 + (qJD(1) * t150 + 0.2e1 * t136 * t146) * t119;
	t1 = [t157 * t181 + (qJD(3) * t150 - t175 * t182) * t119, 0, t89, 0, 0, 0; (-t99 * t165 + (0.2e1 * t167 + (qJD(1) * t92 + t88) * t197) * t130) * t136 + ((-t92 * t165 + (t92 * t170 + ((0.2e1 * t130 * t191 - t93 * t152 - t199 * t173) * t112 + (t157 * t183 + t130 * t93 + (t125 * t160 - (t93 - 0.2e1 * t172) * t130) * t119) * t113) * t94 * t138) * t130) * t100 + (t92 * t169 + (-t99 + ((-t133 + t134) * t113 * t163 + t199 * t164) * t100) * qJD(1)) * t130 * t94) * t138, 0, (-t99 * t94 * t175 + (-0.2e1 * t167 + (-qJD(3) * t90 - t88) * t197) * t138) * t131 + (((-qJD(3) * t99 + t90 * t169) * t138 + (t90 * t175 + (-(-t104 * t174 - t136 * t89) * t113 - ((t104 * t136 - 0.1e1) * t93 + (-t104 + t136) * qJD(3)) * t112) * t181) * t100) * t94 + (t90 * t170 - ((t89 - t174) * t112 + (t155 * t104 + t156) * t113) * t94 * t131) * t100 * t138) * t130, 0, 0, 0; 0.2e1 * (t109 * t147 + t115 * t189) * t195 + (0.2e1 * t115 * t166 - t154 * t109 * t177 + (t130 * t172 + t149) * t190 + (t147 * t98 + t115 * t97 - t149 * t186 - (qJD(3) * t130 * t137 + t154 * t135) * t116 * t136) * t110) * t105, 0, t148 * t168 * t181 + (t148 * t131 * t171 + (-t148 * t175 + ((-qJD(4) * t109 - 0.2e1 * t166) * t137 + (-t137 * t97 + (-qJD(4) * t116 + t98) * t135) * t110) * t138) * t130) * t105, t168 + 0.2e1 * (-t105 * t110 * t97 + (-t105 * t194 - t110 * t195) * t116) * t116, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:56
	% EndTime: 2019-10-10 01:14:57
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (2618->95), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->94)
	t145 = pkin(9) + qJ(3);
	t141 = sin(t145);
	t137 = t141 ^ 2;
	t143 = cos(t145);
	t139 = 0.1e1 / t143 ^ 2;
	t195 = t137 * t139;
	t149 = sin(qJ(1));
	t212 = 0.2e1 * t149;
	t211 = t141 * t195;
	t147 = t149 ^ 2;
	t132 = t147 * t195 + 0.1e1;
	t130 = 0.1e1 / t132;
	t138 = 0.1e1 / t143;
	t150 = cos(qJ(1));
	t184 = qJD(1) * t150;
	t172 = t141 * t184;
	t182 = qJD(3) * t149;
	t104 = (-(-t143 * t182 - t172) * t138 + t182 * t195) * t130;
	t210 = t104 - t182;
	t146 = qJ(4) + pkin(10);
	t142 = sin(t146);
	t187 = t149 * t142;
	t144 = cos(t146);
	t189 = t144 * t150;
	t126 = t143 * t189 + t187;
	t188 = t149 * t141;
	t129 = atan2(-t188, -t143);
	t128 = cos(t129);
	t127 = sin(t129);
	t175 = t127 * t188;
	t113 = -t128 * t143 - t175;
	t110 = 0.1e1 / t113;
	t120 = 0.1e1 / t126;
	t111 = 0.1e1 / t113 ^ 2;
	t121 = 0.1e1 / t126 ^ 2;
	t209 = -0.2e1 * t141;
	t208 = t130 - 0.1e1;
	t197 = t128 * t141;
	t99 = (-t104 * t149 + qJD(3)) * t197 + (t143 * t210 - t172) * t127;
	t207 = t110 * t111 * t99;
	t160 = t143 * t187 + t189;
	t181 = qJD(3) * t150;
	t171 = t141 * t181;
	t105 = t160 * qJD(1) - qJD(4) * t126 + t142 * t171;
	t186 = t149 * t144;
	t190 = t143 * t150;
	t125 = t142 * t190 - t186;
	t119 = t125 ^ 2;
	t118 = t119 * t121 + 0.1e1;
	t200 = t121 * t125;
	t165 = -qJD(1) * t143 + qJD(4);
	t166 = qJD(4) * t143 - qJD(1);
	t191 = t142 * t150;
	t106 = -t166 * t191 + (t165 * t149 - t171) * t144;
	t204 = t106 * t120 * t121;
	t206 = (-t105 * t200 - t119 * t204) / t118 ^ 2;
	t205 = t104 * t141;
	t203 = t111 * t141;
	t193 = t138 * t141;
	t159 = qJD(3) * (t138 * t211 + t193);
	t163 = t137 * t149 * t184;
	t202 = (t139 * t163 + t147 * t159) / t132 ^ 2;
	t201 = t120 * t142;
	t199 = t125 * t144;
	t198 = t127 * t149;
	t196 = t137 * t138;
	t148 = t150 ^ 2;
	t194 = t137 * t148;
	t192 = t141 * t150;
	t185 = qJD(1) * t149;
	t183 = qJD(3) * t143;
	t109 = t111 * t194 + 0.1e1;
	t180 = 0.2e1 / t109 ^ 2 * (-t194 * t207 + (t141 * t148 * t183 - t163) * t111);
	t179 = 0.2e1 * t207;
	t178 = -0.2e1 * t206;
	t177 = t125 * t204;
	t176 = t111 * t192;
	t174 = t130 * t196;
	t170 = 0.1e1 + t195;
	t169 = t141 * t180;
	t168 = t202 * t209;
	t167 = t202 * t212;
	t164 = t149 * t174;
	t162 = t170 * t150;
	t161 = t121 * t199 - t201;
	t158 = t141 * t182 + t165 * t150;
	t124 = -t143 * t186 + t191;
	t117 = t170 * t149 * t130;
	t115 = 0.1e1 / t118;
	t107 = 0.1e1 / t109;
	t103 = (t208 * t141 * t127 - t128 * t164) * t150;
	t102 = -t143 * t198 + t197 + (t127 * t143 - t128 * t188) * t117;
	t100 = -t170 * t167 + (qJD(1) * t162 + t159 * t212) * t130;
	t1 = [t138 * t150 * t168 + (qJD(3) * t162 - t185 * t193) * t130, 0, t100, 0, 0, 0; (t110 * t169 + (-t110 * t183 + (qJD(1) * t103 + t99) * t203) * t107) * t149 + (t111 * t169 * t103 + (-((t104 * t164 + t208 * t183 + t168) * t127 + (t167 * t196 - t205 + (t205 + (t209 - t211) * t182) * t130) * t128) * t176 + (-t111 * t183 + t141 * t179) * t103 + (-t110 + ((-t147 + t148) * t128 * t174 + t208 * t175) * t111) * t141 * qJD(1)) * t107) * t150, 0, (t102 * t203 - t110 * t143) * t150 * t180 + ((-t110 * t185 + (-qJD(3) * t102 - t99) * t150 * t111) * t143 + (-t110 * t181 - (-t100 * t128 * t149 - t210 * t127 + (-qJD(3) * t127 + t104 * t198 - t128 * t184) * t117) * t176 + (t111 * t185 + t150 * t179) * t102 - ((t100 - t184) * t127 + ((-t117 * t149 + 0.1e1) * qJD(3) + (t117 - t149) * t104) * t128) * t111 * t190) * t141) * t107, 0, 0, 0; 0.2e1 * (t120 * t160 + t124 * t200) * t206 + (0.2e1 * t124 * t177 - t166 * t120 * t186 + t158 * t201 + (-t166 * t125 * t187 + t124 * t105 + t106 * t160 - t158 * t199) * t121) * t115, 0, t161 * t178 * t192 + (t161 * t143 * t181 + (-t161 * t185 + ((-qJD(4) * t120 - 0.2e1 * t177) * t144 + (-t105 * t144 + (-qJD(4) * t125 + t106) * t142) * t121) * t150) * t141) * t115, t178 + 0.2e1 * (-t105 * t115 * t121 + (-t115 * t204 - t121 * t206) * t125) * t125, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:56
	% EndTime: 2019-10-10 01:14:57
	% DurationCPUTime: 1.53s
	% Computational Cost: add. (6874->125), mult. (6168->272), div. (1114->15), fcn. (7752->9), ass. (0->116)
	t176 = pkin(9) + qJ(3);
	t174 = cos(t176);
	t177 = qJ(4) + pkin(10);
	t173 = sin(t177);
	t251 = sin(qJ(1));
	t211 = t251 * t173;
	t175 = cos(t177);
	t179 = cos(qJ(1));
	t231 = t179 * t175;
	t154 = t174 * t231 + t211;
	t148 = 0.1e1 / t154 ^ 2;
	t172 = sin(t176);
	t165 = t172 ^ 2;
	t178 = t179 ^ 2;
	t239 = t165 * t178;
	t216 = t148 * t239;
	t144 = 0.1e1 + t216;
	t203 = qJD(1) * t251;
	t228 = qJD(3) * t179;
	t207 = t172 * t228;
	t189 = t174 * t203 + t207;
	t202 = t251 * qJD(4);
	t232 = t179 * t173;
	t133 = (-qJD(4) * t174 + qJD(1)) * t232 + (t202 - t189) * t175;
	t147 = 0.1e1 / t154;
	t246 = t133 * t147 * t148;
	t197 = t239 * t246;
	t208 = qJD(3) * t172 * t178;
	t254 = (-t197 + (-t165 * t179 * t203 + t174 * t208) * t148) / t144 ^ 2;
	t233 = t172 * t179;
	t150 = t174 * t211 + t231;
	t194 = t173 * t202;
	t225 = qJD(4) * t179;
	t205 = t175 * t225;
	t132 = t150 * qJD(1) + t173 * t207 - t174 * t205 - t194;
	t210 = t251 * t175;
	t153 = t174 * t232 - t210;
	t166 = 0.1e1 / t172;
	t169 = 0.1e1 / t173;
	t170 = 0.1e1 / t173 ^ 2;
	t226 = qJD(4) * t175;
	t206 = t170 * t226;
	t167 = 0.1e1 / t172 ^ 2;
	t229 = qJD(3) * t174;
	t209 = t167 * t229;
	t238 = t166 * t169;
	t253 = (t166 * t206 + t169 * t209) * t153 + t132 * t238;
	t234 = t172 * t173;
	t140 = atan2(-t150, t234);
	t137 = cos(t140);
	t136 = sin(t140);
	t245 = t136 * t150;
	t131 = t137 * t234 - t245;
	t128 = 0.1e1 / t131;
	t129 = 0.1e1 / t131 ^ 2;
	t252 = 0.2e1 * t153;
	t145 = t150 ^ 2;
	t237 = t167 * t170;
	t141 = t145 * t237 + 0.1e1;
	t138 = 0.1e1 / t141;
	t190 = t172 * t226 + t173 * t229;
	t214 = t150 * t237;
	t212 = t172 * t251;
	t195 = qJD(3) * t212;
	t196 = t175 * t203;
	t230 = qJD(1) * t179;
	t134 = t175 * t202 * t174 - t196 + (t230 * t174 - t195 - t225) * t173;
	t217 = t134 * t238;
	t120 = (t190 * t214 - t217) * t138;
	t187 = -t120 * t150 + t190;
	t116 = (-t120 * t234 - t134) * t136 + t187 * t137;
	t130 = t128 * t129;
	t250 = t116 * t130;
	t168 = t166 / t165;
	t171 = t169 * t170;
	t249 = (t134 * t214 + (-t167 * t171 * t226 - t168 * t170 * t229) * t145) / t141 ^ 2;
	t248 = t129 * t153;
	t247 = t132 * t129;
	t244 = t136 * t153;
	t243 = t136 * t172;
	t242 = t137 * t150;
	t241 = t137 * t153;
	t240 = t137 * t174;
	t236 = t167 * t174;
	t235 = t170 * t175;
	t227 = qJD(4) * t173;
	t146 = t153 ^ 2;
	t126 = t129 * t146 + 0.1e1;
	t224 = 0.2e1 * (-t146 * t250 - t153 * t247) / t126 ^ 2;
	t223 = -0.2e1 * t249;
	t222 = 0.2e1 * t254;
	t221 = t130 * t252;
	t220 = t166 * t249;
	t219 = t129 * t244;
	t215 = t150 * t238;
	t213 = t169 * t236;
	t192 = t150 * t213 + t251;
	t127 = t192 * t138;
	t204 = t251 - t127;
	t201 = t128 * t224;
	t200 = t129 * t224;
	t199 = t233 * t252;
	t198 = t169 * t220;
	t152 = t174 * t210 - t232;
	t193 = t150 * t235 - t152 * t169;
	t191 = t148 * t152 * t179 - t251 * t147;
	t142 = 0.1e1 / t144;
	t135 = t154 * qJD(1) - t174 * t194 - t175 * t195 - t205;
	t124 = 0.1e1 / t126;
	t123 = t193 * t166 * t138;
	t119 = (-t136 + (t137 * t215 + t136) * t138) * t153;
	t118 = -t127 * t242 + (t204 * t243 + t240) * t173;
	t117 = t137 * t172 * t175 - t136 * t152 + (-t136 * t234 - t242) * t123;
	t115 = t192 * t223 + (t134 * t213 + t230 + (-t206 * t236 + (-0.2e1 * t168 * t174 ^ 2 - t166) * t169 * qJD(3)) * t150) * t138;
	t113 = -0.2e1 * t193 * t220 + (-t193 * t209 + (t134 * t235 - t135 * t169 + (t152 * t235 + (-0.2e1 * t171 * t175 ^ 2 - t169) * t150) * qJD(4)) * t166) * t138;
	t1 = [t253 * t138 + t198 * t252, 0, t115, t113, 0, 0; t150 * t201 + (-t134 * t128 + (t116 * t150 + t119 * t132) * t129) * t124 + (t119 * t200 + (0.2e1 * t119 * t250 + (t132 * t138 - t132 - (-t120 * t138 * t215 + t223) * t153) * t129 * t136 + (-(-0.2e1 * t150 * t198 - t120) * t248 + (-(t120 + t217) * t153 + t253 * t150) * t129 * t138) * t137) * t124) * t153, 0, t118 * t153 * t200 + (-(-t115 * t242 + (t120 * t245 - t134 * t137) * t127) * t248 + (t116 * t221 + t247) * t118 + (-t128 * t233 - (-t127 * t243 + t136 * t212 + t240) * t248) * t226) * t124 + (t201 * t233 + ((-t128 * t228 - (t204 * qJD(3) - t120) * t219) * t174 + (t128 * t203 + (t179 * t116 - (-t115 + t230) * t244 - (t204 * t120 - qJD(3)) * t241) * t129) * t172) * t124) * t173, (t117 * t248 - t128 * t154) * t224 + (t117 * t247 + t133 * t128 + (t117 * t221 - t129 * t154) * t116 - (t175 * t229 - t172 * t227 - t113 * t150 - t123 * t134 + (-t123 * t234 - t152) * t120) * t129 * t241 - (-t135 + (-t113 * t173 - t120 * t175) * t172 - t187 * t123) * t219) * t124, 0, 0; t191 * t172 * t222 + (-t191 * t229 + ((qJD(1) * t147 + 0.2e1 * t152 * t246) * t179 + (-t251 * t133 - t135 * t179 + t152 * t203) * t148) * t172) * t142, 0, (t147 * t174 * t179 + t175 * t216) * t222 + (0.2e1 * t175 * t197 + t189 * t147 + ((t133 * t179 - 0.2e1 * t175 * t208) * t174 + (t178 * t227 + 0.2e1 * t179 * t196) * t165) * t148) * t142, t148 * t199 * t254 + (t199 * t246 + (t132 * t233 + (t172 * t203 - t174 * t228) * t153) * t148) * t142, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end