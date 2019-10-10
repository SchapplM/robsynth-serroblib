% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR4
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
%   Wie in S6RRRPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:58
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:50
	% EndTime: 2019-10-10 11:58:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:50
	% EndTime: 2019-10-10 11:58:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:50
	% EndTime: 2019-10-10 11:58:51
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (3413->87), mult. (3414->191), div. (723->12), fcn. (4027->9), ass. (0->91)
	t151 = sin(qJ(1));
	t206 = 0.2e1 * t151;
	t148 = qJ(2) + qJ(3);
	t144 = cos(t148);
	t150 = cos(pkin(11));
	t180 = t151 * t150;
	t149 = sin(pkin(11));
	t152 = cos(qJ(1));
	t184 = t149 * t152;
	t130 = t144 * t184 - t180;
	t124 = t130 ^ 2;
	t181 = t151 * t149;
	t183 = t150 * t152;
	t131 = t144 * t183 + t181;
	t126 = 0.1e1 / t131 ^ 2;
	t121 = t124 * t126 + 0.1e1;
	t128 = -t144 * t181 - t183;
	t143 = sin(t148);
	t145 = qJD(2) + qJD(3);
	t185 = t145 * t152;
	t170 = t143 * t185;
	t122 = t128 * qJD(1) - t149 * t170;
	t196 = t126 * t130;
	t129 = -t144 * t180 + t184;
	t123 = t129 * qJD(1) - t150 * t170;
	t125 = 0.1e1 / t131;
	t198 = t123 * t125 * t126;
	t205 = (t122 * t196 - t124 * t198) / t121 ^ 2;
	t146 = t151 ^ 2;
	t139 = t143 ^ 2;
	t141 = 0.1e1 / t144 ^ 2;
	t192 = t139 * t141;
	t137 = t146 * t192 + 0.1e1;
	t135 = 0.1e1 / t137;
	t140 = 0.1e1 / t144;
	t178 = qJD(1) * t152;
	t169 = t143 * t178;
	t186 = t145 * t151;
	t172 = t141 * t186;
	t109 = (-(-t144 * t186 - t169) * t140 + t139 * t172) * t135;
	t204 = t109 - t186;
	t182 = t151 * t143;
	t134 = atan2(-t182, -t144);
	t133 = cos(t134);
	t132 = sin(t134);
	t173 = t132 * t182;
	t117 = -t133 * t144 - t173;
	t114 = 0.1e1 / t117;
	t115 = 0.1e1 / t117 ^ 2;
	t203 = t135 - 0.1e1;
	t194 = t133 * t143;
	t104 = (-t109 * t151 + t145) * t194 + (t204 * t144 - t169) * t132;
	t202 = t104 * t114 * t115;
	t138 = t143 * t139;
	t189 = t140 * t143;
	t160 = t145 * (t138 * t140 * t141 + t189);
	t190 = t139 * t151;
	t163 = t178 * t190;
	t201 = (t141 * t163 + t146 * t160) / t137 ^ 2;
	t200 = t115 * t143;
	t199 = t115 * t152;
	t197 = t125 * t149;
	t195 = t132 * t151;
	t193 = t139 * t140;
	t147 = t152 ^ 2;
	t191 = t139 * t147;
	t188 = t143 * t152;
	t187 = t144 * t145;
	t179 = qJD(1) * t151;
	t112 = t115 * t191 + 0.1e1;
	t177 = 0.2e1 * (-t191 * t202 + (t143 * t147 * t187 - t163) * t115) / t112 ^ 2;
	t176 = 0.2e1 * t202;
	t175 = t115 * t188;
	t174 = t130 * t198;
	t171 = t145 * t182;
	t168 = 0.1e1 + t192;
	t167 = t143 * t177;
	t166 = -0.2e1 * t143 * t201;
	t165 = t201 * t206;
	t164 = t133 * t135 * t193;
	t162 = t168 * t152;
	t161 = t150 * t196 - t197;
	t119 = 0.1e1 / t121;
	t118 = t168 * t151 * t135;
	t110 = 0.1e1 / t112;
	t108 = (t203 * t143 * t132 - t151 * t164) * t152;
	t106 = -t144 * t195 + t194 + (t132 * t144 - t133 * t182) * t118;
	t105 = -t168 * t165 + (qJD(1) * t162 + t160 * t206) * t135;
	t102 = -0.2e1 * t161 * t188 * t205 + (t161 * t144 * t185 + (-0.2e1 * t174 * t183 + t179 * t197 + (t123 * t184 + (t122 * t152 - t130 * t179) * t150) * t126) * t143) * t119;
	t101 = (t106 * t200 - t114 * t144) * t152 * t177 + ((-t114 * t179 + (-t106 * t145 - t104) * t199) * t144 + (-t114 * t185 - (-t105 * t133 * t151 - t204 * t132 + (t109 * t195 - t132 * t145 - t133 * t178) * t118) * t175 + (t115 * t179 + t152 * t176) * t106 - ((t105 - t178) * t132 + ((-t118 * t151 + 0.1e1) * t145 + (t118 - t151) * t109) * t133) * t144 * t199) * t143) * t110;
	t1 = [t140 * t152 * t166 + (t145 * t162 - t179 * t189) * t135, t105, t105, 0, 0, 0; (t114 * t167 + (-t114 * t187 + (qJD(1) * t108 + t104) * t200) * t110) * t151 + (t115 * t167 * t108 + (-((t166 - t187 + (t109 * t140 * t190 + t187) * t135) * t132 + (t165 * t193 - t109 * t143 + (-t138 * t172 + (t109 - 0.2e1 * t186) * t143) * t135) * t133) * t175 + (-t115 * t187 + t143 * t176) * t108 + (-t114 + ((-t146 + t147) * t164 + t203 * t173) * t115) * t143 * qJD(1)) * t110) * t152, t101, t101, 0, 0, 0; 0.2e1 * (-t125 * t128 + t129 * t196) * t205 + ((-t130 * qJD(1) + t149 * t171) * t125 + 0.2e1 * t129 * t174 + (-t128 * t123 - (-t131 * qJD(1) + t150 * t171) * t130 - t129 * t122) * t126) * t119, t102, t102, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:50
	% EndTime: 2019-10-10 11:58:51
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (4131->97), mult. (3810->206), div. (753->12), fcn. (4455->9), ass. (0->95)
	t170 = sin(qJ(1));
	t230 = 0.2e1 * t170;
	t167 = t170 ^ 2;
	t169 = qJ(2) + qJ(3);
	t163 = sin(t169);
	t157 = t163 ^ 2;
	t164 = cos(t169);
	t159 = 0.1e1 / t164 ^ 2;
	t216 = t157 * t159;
	t153 = t167 * t216 + 0.1e1;
	t150 = 0.1e1 / t153;
	t158 = 0.1e1 / t164;
	t171 = cos(qJ(1));
	t202 = qJD(1) * t171;
	t192 = t163 * t202;
	t166 = qJD(2) + qJD(3);
	t208 = t166 * t170;
	t195 = t159 * t208;
	t124 = (-(-t164 * t208 - t192) * t158 + t157 * t195) * t150;
	t229 = t124 - t208;
	t165 = pkin(11) + qJ(5);
	t162 = cos(t165);
	t161 = sin(t165);
	t206 = t170 * t161;
	t209 = t164 * t171;
	t146 = t162 * t209 + t206;
	t204 = t170 * t163;
	t149 = atan2(-t204, -t164);
	t148 = cos(t149);
	t147 = sin(t149);
	t196 = t147 * t204;
	t134 = -t148 * t164 - t196;
	t131 = 0.1e1 / t134;
	t140 = 0.1e1 / t146;
	t132 = 0.1e1 / t134 ^ 2;
	t141 = 0.1e1 / t146 ^ 2;
	t228 = t150 - 0.1e1;
	t218 = t148 * t163;
	t119 = (-t124 * t170 + t166) * t218 + (t229 * t164 - t192) * t147;
	t227 = t119 * t131 * t132;
	t181 = t162 * t171 + t164 * t206;
	t207 = t166 * t171;
	t193 = t163 * t207;
	t125 = t181 * qJD(1) - t146 * qJD(5) + t161 * t193;
	t205 = t170 * t162;
	t145 = t161 * t209 - t205;
	t139 = t145 ^ 2;
	t137 = t139 * t141 + 0.1e1;
	t221 = t141 * t145;
	t186 = -qJD(1) * t164 + qJD(5);
	t187 = qJD(5) * t164 - qJD(1);
	t212 = t161 * t171;
	t126 = -t187 * t212 + (t186 * t170 - t193) * t162;
	t225 = t126 * t140 * t141;
	t226 = (-t125 * t221 - t139 * t225) / t137 ^ 2;
	t156 = t163 * t157;
	t213 = t158 * t163;
	t180 = t166 * (t156 * t158 * t159 + t213);
	t214 = t157 * t170;
	t184 = t202 * t214;
	t224 = (t159 * t184 + t167 * t180) / t153 ^ 2;
	t223 = t132 * t163;
	t222 = t140 * t161;
	t220 = t145 * t162;
	t219 = t147 * t170;
	t217 = t157 * t158;
	t168 = t171 ^ 2;
	t215 = t157 * t168;
	t211 = t163 * t171;
	t210 = t164 * t166;
	t203 = qJD(1) * t170;
	t129 = t132 * t215 + 0.1e1;
	t201 = 0.2e1 * (-t215 * t227 + (t163 * t168 * t210 - t184) * t132) / t129 ^ 2;
	t200 = 0.2e1 * t227;
	t199 = -0.2e1 * t226;
	t198 = t132 * t211;
	t197 = t145 * t225;
	t191 = 0.1e1 + t216;
	t190 = t163 * t201;
	t189 = -0.2e1 * t163 * t224;
	t188 = t224 * t230;
	t185 = t148 * t150 * t217;
	t183 = t191 * t171;
	t182 = t141 * t220 - t222;
	t179 = t166 * t204 + t186 * t171;
	t144 = -t164 * t205 + t212;
	t138 = t191 * t170 * t150;
	t135 = 0.1e1 / t137;
	t127 = 0.1e1 / t129;
	t123 = (t228 * t163 * t147 - t170 * t185) * t171;
	t122 = -t164 * t219 + t218 + (t147 * t164 - t148 * t204) * t138;
	t120 = -t191 * t188 + (qJD(1) * t183 + t180 * t230) * t150;
	t117 = t182 * t199 * t211 + (t182 * t164 * t207 + (-t182 * t203 + ((-qJD(5) * t140 - 0.2e1 * t197) * t162 + (-t125 * t162 + (-qJD(5) * t145 + t126) * t161) * t141) * t171) * t163) * t135;
	t116 = (t122 * t223 - t131 * t164) * t171 * t201 + ((-t131 * t203 + (-t122 * t166 - t119) * t171 * t132) * t164 + (-t131 * t207 - (-t120 * t148 * t170 - t229 * t147 + (t124 * t219 - t147 * t166 - t148 * t202) * t138) * t198 + (t132 * t203 + t171 * t200) * t122 - ((t120 - t202) * t147 + ((-t138 * t170 + 0.1e1) * t166 + (t138 - t170) * t124) * t148) * t132 * t209) * t163) * t127;
	t1 = [t171 * t158 * t189 + (t166 * t183 - t203 * t213) * t150, t120, t120, 0, 0, 0; (t131 * t190 + (-t131 * t210 + (qJD(1) * t123 + t119) * t223) * t127) * t170 + (t132 * t190 * t123 + (-((t189 - t210 + (t124 * t158 * t214 + t210) * t150) * t147 + (t188 * t217 - t124 * t163 + (-t156 * t195 + (t124 - 0.2e1 * t208) * t163) * t150) * t148) * t198 + (-t132 * t210 + t163 * t200) * t123 + (-t131 + ((-t167 + t168) * t185 + t228 * t196) * t132) * t163 * qJD(1)) * t127) * t171, t116, t116, 0, 0, 0; 0.2e1 * (t140 * t181 + t144 * t221) * t226 + (0.2e1 * t144 * t197 - t187 * t140 * t205 + t179 * t222 + (-t187 * t145 * t206 + t144 * t125 + t126 * t181 - t179 * t220) * t141) * t135, t117, t117, 0, t199 + 0.2e1 * (-t125 * t141 * t135 + (-t135 * t225 - t141 * t226) * t145) * t145, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:50
	% EndTime: 2019-10-10 11:58:51
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (5000->98), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->98)
	t196 = sin(qJ(1));
	t257 = 0.2e1 * t196;
	t193 = t196 ^ 2;
	t195 = qJ(2) + qJ(3);
	t189 = sin(t195);
	t184 = t189 ^ 2;
	t190 = cos(t195);
	t186 = 0.1e1 / t190 ^ 2;
	t242 = t184 * t186;
	t178 = t193 * t242 + 0.1e1;
	t176 = 0.1e1 / t178;
	t185 = 0.1e1 / t190;
	t197 = cos(qJ(1));
	t228 = qJD(1) * t197;
	t218 = t189 * t228;
	t192 = qJD(2) + qJD(3);
	t236 = t192 * t196;
	t221 = t186 * t236;
	t151 = (-(-t190 * t236 - t218) * t185 + t184 * t221) * t176;
	t256 = t151 - t236;
	t188 = pkin(11) + qJ(5) + qJ(6);
	t182 = cos(t188);
	t230 = t197 * t182;
	t181 = sin(t188);
	t234 = t196 * t181;
	t171 = t190 * t230 + t234;
	t232 = t196 * t189;
	t175 = atan2(-t232, -t190);
	t174 = cos(t175);
	t173 = sin(t175);
	t222 = t173 * t232;
	t162 = -t174 * t190 - t222;
	t159 = 0.1e1 / t162;
	t165 = 0.1e1 / t171;
	t160 = 0.1e1 / t162 ^ 2;
	t166 = 0.1e1 / t171 ^ 2;
	t255 = t176 - 0.1e1;
	t244 = t174 * t189;
	t144 = (-t151 * t196 + t192) * t244 + (t256 * t190 - t218) * t173;
	t254 = t144 * t159 * t160;
	t191 = qJD(5) + qJD(6);
	t207 = t190 * t234 + t230;
	t235 = t192 * t197;
	t219 = t189 * t235;
	t149 = t207 * qJD(1) - t171 * t191 + t181 * t219;
	t231 = t197 * t181;
	t233 = t196 * t182;
	t170 = t190 * t231 - t233;
	t164 = t170 ^ 2;
	t157 = t164 * t166 + 0.1e1;
	t247 = t166 * t170;
	t212 = -qJD(1) * t190 + t191;
	t213 = t190 * t191 - qJD(1);
	t150 = -t213 * t231 + (t212 * t196 - t219) * t182;
	t252 = t150 * t165 * t166;
	t253 = (-t149 * t247 - t164 * t252) / t157 ^ 2;
	t183 = t189 * t184;
	t239 = t185 * t189;
	t206 = t192 * (t183 * t185 * t186 + t239);
	t240 = t184 * t196;
	t210 = t228 * t240;
	t251 = (t186 * t210 + t193 * t206) / t178 ^ 2;
	t250 = t160 * t189;
	t249 = t160 * t197;
	t248 = t165 * t181;
	t246 = t170 * t182;
	t245 = t173 * t196;
	t243 = t184 * t185;
	t194 = t197 ^ 2;
	t241 = t184 * t194;
	t238 = t189 * t197;
	t237 = t190 * t192;
	t229 = qJD(1) * t196;
	t154 = t160 * t241 + 0.1e1;
	t227 = 0.2e1 * (-t241 * t254 + (t189 * t194 * t237 - t210) * t160) / t154 ^ 2;
	t226 = 0.2e1 * t254;
	t225 = -0.2e1 * t253;
	t224 = t160 * t238;
	t223 = t170 * t252;
	t217 = 0.1e1 + t242;
	t216 = t189 * t227;
	t215 = -0.2e1 * t189 * t251;
	t214 = t251 * t257;
	t211 = t174 * t176 * t243;
	t209 = t217 * t197;
	t208 = t166 * t246 - t248;
	t205 = t192 * t232 + t212 * t197;
	t169 = -t190 * t233 + t231;
	t163 = t217 * t196 * t176;
	t155 = 0.1e1 / t157;
	t152 = 0.1e1 / t154;
	t148 = (t255 * t189 * t173 - t196 * t211) * t197;
	t147 = -t190 * t245 + t244 + (t173 * t190 - t174 * t232) * t163;
	t145 = -t217 * t214 + (qJD(1) * t209 + t206 * t257) * t176;
	t142 = t225 + 0.2e1 * (-t149 * t166 * t155 + (-t155 * t252 - t166 * t253) * t170) * t170;
	t141 = t208 * t225 * t238 + (t208 * t190 * t235 + (-t208 * t229 + ((-t165 * t191 - 0.2e1 * t223) * t182 + (-t149 * t182 + (-t170 * t191 + t150) * t181) * t166) * t197) * t189) * t155;
	t140 = (t147 * t250 - t159 * t190) * t197 * t227 + ((-t159 * t229 + (-t147 * t192 - t144) * t249) * t190 + (-t159 * t235 - (-t145 * t174 * t196 - t256 * t173 + (t151 * t245 - t173 * t192 - t174 * t228) * t163) * t224 + (t160 * t229 + t197 * t226) * t147 - ((t145 - t228) * t173 + ((-t163 * t196 + 0.1e1) * t192 + (t163 - t196) * t151) * t174) * t190 * t249) * t189) * t152;
	t1 = [t197 * t185 * t215 + (t192 * t209 - t229 * t239) * t176, t145, t145, 0, 0, 0; (t159 * t216 + (-t159 * t237 + (qJD(1) * t148 + t144) * t250) * t152) * t196 + (t160 * t216 * t148 + (-((t215 - t237 + (t151 * t185 * t240 + t237) * t176) * t173 + (t214 * t243 - t151 * t189 + (-t183 * t221 + (t151 - 0.2e1 * t236) * t189) * t176) * t174) * t224 + (-t160 * t237 + t189 * t226) * t148 + (-t159 + ((-t193 + t194) * t211 + t255 * t222) * t160) * t189 * qJD(1)) * t152) * t197, t140, t140, 0, 0, 0; 0.2e1 * (t165 * t207 + t169 * t247) * t253 + (0.2e1 * t169 * t223 - t213 * t165 * t233 + t205 * t248 + (-t213 * t170 * t234 + t169 * t149 + t150 * t207 - t205 * t246) * t166) * t155, t141, t141, 0, t142, t142;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end