% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR8
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
%   Wie in S6RPRRRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:55
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (3360->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->97)
	t168 = cos(qJ(1));
	t228 = 0.2e1 * t168;
	t163 = t168 ^ 2;
	t164 = qJ(3) + qJ(4);
	t159 = sin(t164);
	t155 = 0.1e1 / t159 ^ 2;
	t160 = cos(t164);
	t158 = t160 ^ 2;
	t211 = t155 * t158;
	t152 = t163 * t211 + 0.1e1;
	t150 = 0.1e1 / t152;
	t154 = 0.1e1 / t159;
	t166 = sin(qJ(1));
	t200 = qJD(1) * t166;
	t192 = t160 * t200;
	t161 = qJD(3) + qJD(4);
	t206 = t161 * t168;
	t193 = t155 * t206;
	t124 = ((t159 * t206 + t192) * t154 + t158 * t193) * t150;
	t227 = -t124 + t206;
	t202 = t168 * t160;
	t149 = atan2(-t202, t159);
	t147 = sin(t149);
	t148 = cos(t149);
	t134 = -t147 * t202 + t148 * t159;
	t131 = 0.1e1 / t134;
	t167 = cos(qJ(5));
	t203 = t166 * t167;
	t165 = sin(qJ(5));
	t205 = t165 * t168;
	t144 = t159 * t203 + t205;
	t140 = 0.1e1 / t144;
	t132 = 0.1e1 / t134 ^ 2;
	t141 = 0.1e1 / t144 ^ 2;
	t162 = t166 ^ 2;
	t210 = t158 * t162;
	t127 = t132 * t210 + 0.1e1;
	t199 = qJD(1) * t168;
	t184 = t158 * t166 * t199;
	t208 = t160 * t161;
	t214 = t148 * t160;
	t223 = t124 * t168;
	t119 = (t161 - t223) * t214 + (t227 * t159 + t192) * t147;
	t225 = t119 * t131 * t132;
	t226 = (-t210 * t225 + (-t159 * t162 * t208 + t184) * t132) / t127 ^ 2;
	t187 = qJD(1) * t159 + qJD(5);
	t180 = t187 * t168;
	t188 = qJD(5) * t159 + qJD(1);
	t181 = t188 * t167;
	t207 = t161 * t166;
	t128 = t166 * t181 + (t160 * t207 + t180) * t165;
	t201 = t168 * t167;
	t204 = t166 * t165;
	t143 = t159 * t204 - t201;
	t139 = t143 ^ 2;
	t138 = t139 * t141 + 0.1e1;
	t217 = t141 * t143;
	t182 = t188 * t165;
	t129 = t167 * t180 + (t167 * t208 - t182) * t166;
	t222 = t129 * t140 * t141;
	t224 = (t128 * t217 - t139 * t222) / t138 ^ 2;
	t157 = t160 * t158;
	t212 = t154 * t160;
	t178 = t161 * (-t154 * t155 * t157 - t212);
	t221 = (-t155 * t184 + t163 * t178) / t152 ^ 2;
	t220 = t132 * t160;
	t219 = t132 * t166;
	t218 = t140 * t165;
	t216 = t143 * t167;
	t215 = t147 * t168;
	t213 = t154 * t158;
	t209 = t159 * t161;
	t198 = -0.2e1 * t225;
	t197 = 0.2e1 * t224;
	t196 = t160 * t226;
	t195 = t160 * t221;
	t194 = t160 * t219;
	t191 = 0.1e1 + t211;
	t190 = 0.2e1 * t143 * t222;
	t189 = t221 * t228;
	t186 = t148 * t150 * t213;
	t185 = (-t150 + 0.1e1) * t160 * t147;
	t183 = t191 * t166;
	t179 = t141 * t216 - t218;
	t177 = t179 * t166;
	t176 = t161 * t202 - t187 * t166;
	t146 = t159 * t201 - t204;
	t145 = t159 * t205 + t203;
	t136 = 0.1e1 / t138;
	t135 = t191 * t168 * t150;
	t125 = 0.1e1 / t127;
	t123 = (-t168 * t186 + t185) * t166;
	t121 = t159 * t215 + t214 + (-t147 * t159 - t148 * t202) * t135;
	t120 = -t191 * t189 + (-qJD(1) * t183 + t178 * t228) * t150;
	t117 = t160 * t177 * t197 + (t177 * t209 + (-t179 * t199 + ((qJD(5) * t140 + t190) * t167 + (-t128 * t167 + (qJD(5) * t143 - t129) * t165) * t141) * t166) * t160) * t136;
	t116 = 0.2e1 * (-t121 * t220 - t131 * t159) * t166 * t226 + ((t131 * t199 + (-t121 * t161 - t119) * t219) * t159 + (t131 * t207 + (-t120 * t148 * t168 + t227 * t147 + (t124 * t215 - t147 * t161 + t148 * t200) * t135) * t194 + (t132 * t199 + t166 * t198) * t121 + ((-t120 - t200) * t147 + ((t135 * t168 - 0.1e1) * t161 + (-t135 + t168) * t124) * t148) * t159 * t219) * t160) * t125;
	t1 = [-0.2e1 * t166 * t154 * t195 + (-t161 * t183 + t199 * t212) * t150, 0, t120, t120, 0, 0; (0.2e1 * t131 * t196 + (t131 * t209 + (qJD(1) * t123 + t119) * t220) * t125) * t168 + (-0.2e1 * t132 * t196 * t123 + (((0.2e1 * t195 - t209 + (t213 * t223 + t209) * t150) * t147 + (t189 * t213 + t124 * t160 + (t157 * t193 + (-t124 + 0.2e1 * t206) * t160) * t150) * t148) * t194 + (-t132 * t209 + t160 * t198) * t123 + (t131 + ((t162 - t163) * t186 + t168 * t185) * t132) * t160 * qJD(1)) * t125) * t166, 0, t116, t116, 0, 0; (-t140 * t145 + t146 * t217) * t197 + (t146 * t190 + t168 * t140 * t181 + t176 * t218 + (t168 * t143 * t182 - t146 * t128 - t145 * t129 - t176 * t216) * t141) * t136, 0, t117, t117, -0.2e1 * t224 + 0.2e1 * (t128 * t136 * t141 + (-t136 * t222 - t141 * t224) * t143) * t143, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:55
	% DurationCPUTime: 1.05s
	% Computational Cost: add. (4138->96), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->99)
	t196 = cos(qJ(1));
	t255 = 0.2e1 * t196;
	t192 = t196 ^ 2;
	t194 = qJ(3) + qJ(4);
	t186 = sin(t194);
	t181 = 0.1e1 / t186 ^ 2;
	t188 = cos(t194);
	t184 = t188 ^ 2;
	t238 = t181 * t184;
	t178 = t192 * t238 + 0.1e1;
	t176 = 0.1e1 / t178;
	t180 = 0.1e1 / t186;
	t195 = sin(qJ(1));
	t228 = qJD(1) * t195;
	t220 = t188 * t228;
	t190 = qJD(3) + qJD(4);
	t234 = t190 * t196;
	t221 = t181 * t234;
	t150 = ((t186 * t234 + t220) * t180 + t184 * t221) * t176;
	t254 = -t150 + t234;
	t229 = t196 * t188;
	t175 = atan2(-t229, t186);
	t173 = sin(t175);
	t174 = cos(t175);
	t160 = -t173 * t229 + t174 * t186;
	t157 = 0.1e1 / t160;
	t193 = qJ(5) + qJ(6);
	t185 = sin(t193);
	t231 = t196 * t185;
	t187 = cos(t193);
	t232 = t195 * t187;
	t170 = t186 * t232 + t231;
	t166 = 0.1e1 / t170;
	t158 = 0.1e1 / t160 ^ 2;
	t167 = 0.1e1 / t170 ^ 2;
	t191 = t195 ^ 2;
	t237 = t184 * t191;
	t155 = t158 * t237 + 0.1e1;
	t227 = qJD(1) * t196;
	t212 = t184 * t195 * t227;
	t236 = t186 * t190;
	t241 = t174 * t188;
	t250 = t150 * t196;
	t145 = (t190 - t250) * t241 + (t254 * t186 + t220) * t173;
	t252 = t145 * t157 * t158;
	t253 = (-t237 * t252 + (-t188 * t191 * t236 + t212) * t158) / t155 ^ 2;
	t189 = qJD(5) + qJD(6);
	t215 = qJD(1) * t186 + t189;
	t235 = t190 * t195;
	t204 = t235 * t188 + t196 * t215;
	t216 = t186 * t189 + qJD(1);
	t209 = t187 * t216;
	t151 = t185 * t204 + t195 * t209;
	t230 = t196 * t187;
	t233 = t195 * t185;
	t169 = t186 * t233 - t230;
	t165 = t169 ^ 2;
	t164 = t165 * t167 + 0.1e1;
	t244 = t167 * t169;
	t210 = t185 * t216;
	t152 = t187 * t204 - t195 * t210;
	t249 = t152 * t166 * t167;
	t251 = (t151 * t244 - t165 * t249) / t164 ^ 2;
	t183 = t188 * t184;
	t239 = t180 * t188;
	t207 = t190 * (-t180 * t181 * t183 - t239);
	t248 = (-t181 * t212 + t192 * t207) / t178 ^ 2;
	t247 = t158 * t188;
	t246 = t158 * t195;
	t245 = t166 * t185;
	t243 = t169 * t187;
	t242 = t173 * t196;
	t240 = t180 * t184;
	t226 = -0.2e1 * t252;
	t225 = 0.2e1 * t251;
	t224 = t188 * t253;
	t223 = t188 * t248;
	t222 = t188 * t246;
	t219 = 0.1e1 + t238;
	t218 = t248 * t255;
	t217 = 0.2e1 * t169 * t249;
	t214 = t174 * t176 * t240;
	t213 = (-t176 + 0.1e1) * t188 * t173;
	t211 = t219 * t195;
	t208 = t243 * t167 - t245;
	t206 = t208 * t195;
	t205 = t190 * t229 - t195 * t215;
	t172 = t186 * t230 - t233;
	t171 = t186 * t231 + t232;
	t163 = t219 * t196 * t176;
	t161 = 0.1e1 / t164;
	t153 = 0.1e1 / t155;
	t149 = (-t196 * t214 + t213) * t195;
	t148 = t186 * t242 + t241 + (-t173 * t186 - t174 * t229) * t163;
	t146 = -t219 * t218 + (-qJD(1) * t211 + t207 * t255) * t176;
	t143 = -0.2e1 * t251 + 0.2e1 * (t151 * t167 * t161 + (-t161 * t249 - t167 * t251) * t169) * t169;
	t142 = t188 * t206 * t225 + (t206 * t236 + (-t208 * t227 + ((t166 * t189 + t217) * t187 + (-t151 * t187 + (t169 * t189 - t152) * t185) * t167) * t195) * t188) * t161;
	t141 = 0.2e1 * (-t148 * t247 - t157 * t186) * t195 * t253 + ((t157 * t227 + (-t148 * t190 - t145) * t246) * t186 + (t157 * t235 + (-t146 * t174 * t196 + t254 * t173 + (t150 * t242 - t173 * t190 + t174 * t228) * t163) * t222 + (t158 * t227 + t195 * t226) * t148 + ((-t146 - t228) * t173 + ((t163 * t196 - 0.1e1) * t190 + (-t163 + t196) * t150) * t174) * t186 * t246) * t188) * t153;
	t1 = [-0.2e1 * t195 * t180 * t223 + (-t190 * t211 + t227 * t239) * t176, 0, t146, t146, 0, 0; (0.2e1 * t157 * t224 + (t157 * t236 + (qJD(1) * t149 + t145) * t247) * t153) * t196 + (-0.2e1 * t158 * t224 * t149 + (((0.2e1 * t223 - t236 + (t240 * t250 + t236) * t176) * t173 + (t218 * t240 + t150 * t188 + (t183 * t221 + (-t150 + 0.2e1 * t234) * t188) * t176) * t174) * t222 + (-t158 * t236 + t188 * t226) * t149 + (t157 + ((t191 - t192) * t214 + t196 * t213) * t158) * t188 * qJD(1)) * t153) * t195, 0, t141, t141, 0, 0; (-t166 * t171 + t172 * t244) * t225 + (t172 * t217 + t196 * t166 * t209 + t205 * t245 + (t169 * t196 * t210 - t172 * t151 - t171 * t152 - t205 * t243) * t167) * t161, 0, t142, t142, t143, t143;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end