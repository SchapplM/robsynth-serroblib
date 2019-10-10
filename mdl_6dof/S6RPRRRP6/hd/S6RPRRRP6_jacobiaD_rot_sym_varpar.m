% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP6
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
%   Wie in S6RPRRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:09
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
	t136 = sin(qJ(1));
	t133 = t136 ^ 2;
	t132 = pkin(10) + qJ(3);
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
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:09
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (2887->96), mult. (2734->202), div. (498->12), fcn. (3199->9), ass. (0->96)
	t161 = pkin(10) + qJ(3);
	t157 = sin(t161);
	t153 = t157 ^ 2;
	t158 = cos(t161);
	t155 = 0.1e1 / t158 ^ 2;
	t211 = t153 * t155;
	t166 = sin(qJ(1));
	t229 = 0.2e1 * t166;
	t228 = t157 * t211;
	t163 = t166 ^ 2;
	t149 = t163 * t211 + 0.1e1;
	t147 = 0.1e1 / t149;
	t154 = 0.1e1 / t158;
	t167 = cos(qJ(1));
	t201 = qJD(1) * t167;
	t189 = t157 * t201;
	t199 = qJD(3) * t166;
	t120 = (-(-t158 * t199 - t189) * t154 + t199 * t211) * t147;
	t227 = t120 - t199;
	t165 = qJ(4) + qJ(5);
	t159 = sin(t165);
	t204 = t166 * t159;
	t160 = cos(t165);
	t206 = t160 * t167;
	t142 = t158 * t206 + t204;
	t205 = t166 * t157;
	t145 = atan2(-t205, -t158);
	t144 = cos(t145);
	t143 = sin(t145);
	t192 = t143 * t205;
	t129 = -t144 * t158 - t192;
	t126 = 0.1e1 / t129;
	t136 = 0.1e1 / t142;
	t127 = 0.1e1 / t129 ^ 2;
	t137 = 0.1e1 / t142 ^ 2;
	t226 = -0.2e1 * t157;
	t225 = t147 - 0.1e1;
	t213 = t144 * t157;
	t115 = (-t120 * t166 + qJD(3)) * t213 + (t158 * t227 - t189) * t143;
	t224 = t115 * t126 * t127;
	t162 = qJD(4) + qJD(5);
	t177 = t158 * t204 + t206;
	t198 = qJD(3) * t167;
	t188 = t157 * t198;
	t121 = t177 * qJD(1) - t142 * t162 + t159 * t188;
	t203 = t166 * t160;
	t207 = t159 * t167;
	t141 = t158 * t207 - t203;
	t135 = t141 ^ 2;
	t134 = t135 * t137 + 0.1e1;
	t216 = t137 * t141;
	t182 = -qJD(1) * t158 + t162;
	t183 = t158 * t162 - qJD(1);
	t122 = -t183 * t207 + (t182 * t166 - t188) * t160;
	t221 = t122 * t136 * t137;
	t223 = (-t121 * t216 - t135 * t221) / t134 ^ 2;
	t222 = t120 * t157;
	t220 = t127 * t157;
	t219 = t127 * t167;
	t209 = t154 * t157;
	t176 = qJD(3) * (t154 * t228 + t209);
	t180 = t153 * t166 * t201;
	t218 = (t155 * t180 + t163 * t176) / t149 ^ 2;
	t217 = t136 * t159;
	t215 = t141 * t160;
	t214 = t143 * t166;
	t212 = t153 * t154;
	t164 = t167 ^ 2;
	t210 = t153 * t164;
	t208 = t157 * t167;
	t202 = qJD(1) * t166;
	t200 = qJD(3) * t158;
	t125 = t127 * t210 + 0.1e1;
	t197 = 0.2e1 * (-t210 * t224 + (t157 * t164 * t200 - t180) * t127) / t125 ^ 2;
	t196 = 0.2e1 * t224;
	t195 = -0.2e1 * t223;
	t194 = t141 * t221;
	t193 = t127 * t208;
	t191 = t147 * t212;
	t187 = 0.1e1 + t211;
	t186 = t157 * t197;
	t185 = t218 * t226;
	t184 = t218 * t229;
	t181 = t166 * t191;
	t179 = t187 * t167;
	t178 = t137 * t215 - t217;
	t175 = t157 * t199 + t182 * t167;
	t140 = -t158 * t203 + t207;
	t132 = 0.1e1 / t134;
	t131 = t187 * t166 * t147;
	t123 = 0.1e1 / t125;
	t119 = (t225 * t157 * t143 - t144 * t181) * t167;
	t118 = -t158 * t214 + t213 + (t143 * t158 - t144 * t205) * t131;
	t116 = -t187 * t184 + (qJD(1) * t179 + t176 * t229) * t147;
	t113 = t195 + 0.2e1 * (-t121 * t132 * t137 + (-t132 * t221 - t137 * t223) * t141) * t141;
	t1 = [t154 * t167 * t185 + (qJD(3) * t179 - t202 * t209) * t147, 0, t116, 0, 0, 0; (t126 * t186 + (-t126 * t200 + (qJD(1) * t119 + t115) * t220) * t123) * t166 + (t127 * t186 * t119 + (-((t120 * t181 + t225 * t200 + t185) * t143 + (t184 * t212 - t222 + (t222 + (t226 - t228) * t199) * t147) * t144) * t193 + (-t127 * t200 + t157 * t196) * t119 + (-t126 + ((-t163 + t164) * t144 * t191 + t225 * t192) * t127) * t157 * qJD(1)) * t123) * t167, 0, (t118 * t220 - t126 * t158) * t167 * t197 + ((-t126 * t202 + (-qJD(3) * t118 - t115) * t219) * t158 + (-t126 * t198 - (-t116 * t144 * t166 - t227 * t143 + (-qJD(3) * t143 + t120 * t214 - t144 * t201) * t131) * t193 + (t127 * t202 + t167 * t196) * t118 - ((t116 - t201) * t143 + ((-t131 * t166 + 0.1e1) * qJD(3) + (t131 - t166) * t120) * t144) * t158 * t219) * t157) * t123, 0, 0, 0; 0.2e1 * (t136 * t177 + t140 * t216) * t223 + (0.2e1 * t140 * t194 - t183 * t136 * t203 + t175 * t217 + (-t183 * t141 * t204 + t140 * t121 + t122 * t177 - t175 * t215) * t137) * t132, 0, t178 * t195 * t208 + (t178 * t158 * t198 + (-t178 * t202 + ((-t136 * t162 - 0.2e1 * t194) * t160 + (-t121 * t160 + (-t141 * t162 + t122) * t159) * t137) * t167) * t157) * t132, t113, t113, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:54:08
	% EndTime: 2019-10-10 01:54:09
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (2887->96), mult. (2734->202), div. (498->12), fcn. (3199->9), ass. (0->96)
	t162 = pkin(10) + qJ(3);
	t158 = sin(t162);
	t154 = t158 ^ 2;
	t159 = cos(t162);
	t156 = 0.1e1 / t159 ^ 2;
	t212 = t154 * t156;
	t167 = sin(qJ(1));
	t230 = 0.2e1 * t167;
	t229 = t158 * t212;
	t164 = t167 ^ 2;
	t150 = t164 * t212 + 0.1e1;
	t148 = 0.1e1 / t150;
	t155 = 0.1e1 / t159;
	t168 = cos(qJ(1));
	t202 = qJD(1) * t168;
	t190 = t158 * t202;
	t200 = qJD(3) * t167;
	t121 = (-(-t159 * t200 - t190) * t155 + t200 * t212) * t148;
	t228 = t121 - t200;
	t166 = qJ(4) + qJ(5);
	t161 = cos(t166);
	t204 = t168 * t161;
	t160 = sin(t166);
	t207 = t167 * t160;
	t143 = t159 * t204 + t207;
	t208 = t167 * t158;
	t146 = atan2(-t208, -t159);
	t145 = cos(t146);
	t144 = sin(t146);
	t193 = t144 * t208;
	t130 = -t145 * t159 - t193;
	t127 = 0.1e1 / t130;
	t137 = 0.1e1 / t143;
	t128 = 0.1e1 / t130 ^ 2;
	t138 = 0.1e1 / t143 ^ 2;
	t227 = -0.2e1 * t158;
	t226 = t148 - 0.1e1;
	t214 = t145 * t158;
	t116 = (-t121 * t167 + qJD(3)) * t214 + (t228 * t159 - t190) * t144;
	t225 = t116 * t127 * t128;
	t163 = qJD(4) + qJD(5);
	t178 = t159 * t207 + t204;
	t199 = qJD(3) * t168;
	t189 = t158 * t199;
	t122 = qJD(1) * t178 - t143 * t163 + t160 * t189;
	t205 = t168 * t160;
	t206 = t167 * t161;
	t142 = t159 * t205 - t206;
	t136 = t142 ^ 2;
	t135 = t136 * t138 + 0.1e1;
	t217 = t138 * t142;
	t183 = -qJD(1) * t159 + t163;
	t184 = t159 * t163 - qJD(1);
	t123 = -t184 * t205 + (t167 * t183 - t189) * t161;
	t222 = t123 * t137 * t138;
	t224 = (-t122 * t217 - t136 * t222) / t135 ^ 2;
	t223 = t121 * t158;
	t221 = t128 * t158;
	t220 = t128 * t168;
	t210 = t155 * t158;
	t177 = qJD(3) * (t155 * t229 + t210);
	t181 = t154 * t167 * t202;
	t219 = (t156 * t181 + t164 * t177) / t150 ^ 2;
	t218 = t137 * t160;
	t216 = t142 * t161;
	t215 = t144 * t167;
	t213 = t154 * t155;
	t165 = t168 ^ 2;
	t211 = t154 * t165;
	t209 = t158 * t168;
	t203 = qJD(1) * t167;
	t201 = qJD(3) * t159;
	t126 = t128 * t211 + 0.1e1;
	t198 = 0.2e1 * (-t211 * t225 + (t158 * t165 * t201 - t181) * t128) / t126 ^ 2;
	t197 = 0.2e1 * t225;
	t196 = -0.2e1 * t224;
	t195 = t128 * t209;
	t194 = t142 * t222;
	t192 = t148 * t213;
	t188 = 0.1e1 + t212;
	t187 = t158 * t198;
	t186 = t219 * t227;
	t185 = t219 * t230;
	t182 = t167 * t192;
	t180 = t188 * t168;
	t179 = t216 * t138 - t218;
	t176 = t158 * t200 + t168 * t183;
	t141 = -t159 * t206 + t205;
	t133 = 0.1e1 / t135;
	t132 = t188 * t167 * t148;
	t124 = 0.1e1 / t126;
	t120 = (t144 * t158 * t226 - t145 * t182) * t168;
	t119 = -t159 * t215 + t214 + (t144 * t159 - t145 * t208) * t132;
	t117 = -t188 * t185 + (qJD(1) * t180 + t177 * t230) * t148;
	t114 = t196 + 0.2e1 * (-t122 * t138 * t133 + (-t133 * t222 - t138 * t224) * t142) * t142;
	t1 = [t168 * t155 * t186 + (qJD(3) * t180 - t203 * t210) * t148, 0, t117, 0, 0, 0; (t127 * t187 + (-t127 * t201 + (qJD(1) * t120 + t116) * t221) * t124) * t167 + (t128 * t187 * t120 + (-((t121 * t182 + t201 * t226 + t186) * t144 + (t185 * t213 - t223 + (t223 + (t227 - t229) * t200) * t148) * t145) * t195 + (-t128 * t201 + t158 * t197) * t120 + (-t127 + ((-t164 + t165) * t145 * t192 + t226 * t193) * t128) * t158 * qJD(1)) * t124) * t168, 0, (t119 * t221 - t127 * t159) * t168 * t198 + ((-t127 * t203 + (-qJD(3) * t119 - t116) * t220) * t159 + (-t127 * t199 - (-t117 * t145 * t167 - t228 * t144 + (-qJD(3) * t144 + t121 * t215 - t145 * t202) * t132) * t195 + (t128 * t203 + t168 * t197) * t119 - ((t117 - t202) * t144 + ((-t132 * t167 + 0.1e1) * qJD(3) + (t132 - t167) * t121) * t145) * t159 * t220) * t158) * t124, 0, 0, 0; 0.2e1 * (t137 * t178 + t141 * t217) * t224 + (0.2e1 * t141 * t194 - t184 * t137 * t206 + t176 * t218 + (-t142 * t184 * t207 + t141 * t122 + t123 * t178 - t176 * t216) * t138) * t133, 0, t179 * t196 * t209 + (t179 * t159 * t199 + (-t179 * t203 + ((-t137 * t163 - 0.2e1 * t194) * t161 + (-t122 * t161 + (-t142 * t163 + t123) * t160) * t138) * t168) * t158) * t133, t114, t114, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end