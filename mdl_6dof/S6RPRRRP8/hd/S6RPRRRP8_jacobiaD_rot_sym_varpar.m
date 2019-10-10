% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP8
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
%   Wie in S6RPRRRP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:55
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:52
	% DurationCPUTime: 1.03s
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
	% StartTime: 2019-10-10 01:55:51
	% EndTime: 2019-10-10 01:55:53
	% DurationCPUTime: 1.54s
	% Computational Cost: add. (5529->123), mult. (8382->269), div. (1515->15), fcn. (10508->9), ass. (0->119)
	t188 = qJ(3) + qJ(4);
	t182 = cos(t188);
	t190 = sin(qJ(1));
	t246 = t182 * t190;
	t181 = sin(t188);
	t189 = sin(qJ(5));
	t265 = cos(qJ(1));
	t221 = t265 * t189;
	t191 = cos(qJ(5));
	t239 = t190 * t191;
	t169 = t181 * t221 + t239;
	t247 = t182 * t189;
	t159 = atan2(t169, t247);
	t154 = cos(t159);
	t153 = sin(t159);
	t256 = t153 * t169;
	t148 = t154 * t247 + t256;
	t145 = 0.1e1 / t148;
	t168 = t181 * t239 + t221;
	t163 = 0.1e1 / t168;
	t178 = 0.1e1 / t182;
	t184 = 0.1e1 / t189;
	t146 = 0.1e1 / t148 ^ 2;
	t164 = 0.1e1 / t168 ^ 2;
	t185 = 0.1e1 / t189 ^ 2;
	t220 = t265 * t191;
	t240 = t190 * t189;
	t167 = t181 * t240 - t220;
	t162 = t167 ^ 2;
	t143 = t162 * t146 + 0.1e1;
	t217 = qJD(1) * t265;
	t207 = t181 * t217;
	t201 = t265 * qJD(5) + t207;
	t212 = qJD(5) * t181 + qJD(1);
	t183 = qJD(3) + qJD(4);
	t244 = t183 * t190;
	t225 = t182 * t244;
	t151 = t212 * t239 + (t201 + t225) * t189;
	t259 = t151 * t146;
	t166 = t169 ^ 2;
	t179 = 0.1e1 / t182 ^ 2;
	t249 = t179 * t185;
	t161 = t166 * t249 + 0.1e1;
	t157 = 0.1e1 / t161;
	t236 = qJD(5) * t191;
	t248 = t181 * t183;
	t202 = t182 * t236 - t189 * t248;
	t228 = t169 * t249;
	t208 = t191 * t217;
	t209 = t181 * t220;
	t222 = t265 * t182;
	t210 = t183 * t222;
	t149 = -t189 * t210 - qJD(5) * t209 - t208 + (qJD(1) * t181 + qJD(5)) * t240;
	t250 = t178 * t184;
	t230 = t149 * t250;
	t137 = (-t202 * t228 - t230) * t157;
	t200 = t137 * t169 + t202;
	t133 = (-t137 * t247 - t149) * t153 + t200 * t154;
	t147 = t145 * t146;
	t263 = t133 * t147;
	t264 = (-t162 * t263 + t167 * t259) / t143 ^ 2;
	t177 = t182 ^ 2;
	t187 = t190 ^ 2;
	t251 = t177 * t187;
	t223 = t164 * t251;
	t160 = 0.1e1 + t223;
	t243 = t183 * t191;
	t224 = t182 * t243;
	t152 = t201 * t191 + (-t212 * t189 + t224) * t190;
	t258 = t152 * t163 * t164;
	t211 = t251 * t258;
	t262 = (-t211 + (t177 * t190 * t217 - t182 * t187 * t248) * t164) / t160 ^ 2;
	t180 = t178 / t177;
	t186 = t184 * t185;
	t261 = (-t149 * t228 + (-t179 * t186 * t236 + t180 * t185 * t248) * t166) / t161 ^ 2;
	t260 = t146 * t167;
	t257 = t153 * t167;
	t255 = t153 * t182;
	t254 = t154 * t167;
	t253 = t154 * t169;
	t252 = t154 * t181;
	t245 = t183 * t184;
	t242 = t185 * t191;
	t241 = t190 * t145;
	t238 = qJD(1) * t190;
	t237 = qJD(5) * t189;
	t235 = 0.2e1 * t264;
	t234 = 0.2e1 * t262;
	t233 = -0.2e1 * t261;
	t232 = 0.2e1 * t147 * t167;
	t231 = t146 * t257;
	t229 = t169 * t250;
	t227 = t179 * t181 * t184;
	t226 = t181 * t245;
	t219 = t185 * t236;
	t204 = t169 * t227 + t265;
	t144 = t204 * t157;
	t218 = t265 - t144;
	t216 = -0.2e1 * t145 * t264;
	t215 = t146 * t235;
	t214 = 0.2e1 * t178 * t261;
	t213 = -0.2e1 * t167 * t246;
	t206 = t184 * t214;
	t170 = t209 - t240;
	t205 = t169 * t242 - t170 * t184;
	t203 = t164 * t170 * t190 - t265 * t163;
	t199 = t151 * t250 - (t178 * t219 - t179 * t226) * t167;
	t155 = 0.1e1 / t160;
	t150 = t168 * qJD(1) + t169 * qJD(5) - t191 * t210;
	t141 = 0.1e1 / t143;
	t140 = t205 * t178 * t157;
	t136 = (-t153 + (-t154 * t229 + t153) * t157) * t167;
	t135 = t144 * t253 + (t218 * t255 - t252) * t189;
	t134 = t154 * t182 * t191 + t153 * t170 - (-t153 * t247 + t253) * t140;
	t132 = t204 * t233 + (-t149 * t227 - t238 + (t178 * t245 + (-t179 * t219 + 0.2e1 * t180 * t226) * t181) * t169) * t157;
	t130 = t205 * t214 + (-t205 * t179 * t248 + (t149 * t242 - t150 * t184 + (-t170 * t242 + (0.2e1 * t186 * t191 ^ 2 + t184) * t169) * qJD(5)) * t178) * t157;
	t129 = (t163 * t181 * t190 + t191 * t223) * t234 + (0.2e1 * t191 * t211 + (-t207 - t225) * t163 + ((t152 * t190 + 0.2e1 * t187 * t224) * t181 + (t187 * t237 - 0.2e1 * t190 * t208) * t177) * t164) * t155;
	t128 = t135 * t167 * t215 + (-(t132 * t253 + (-t137 * t256 - t149 * t154) * t144) * t260 + (t133 * t232 - t259) * t135 + (t182 * t241 - (-t144 * t255 + t153 * t222 - t252) * t260) * t236) * t141 + (t216 * t246 + ((-t183 * t241 - (-t218 * t183 + t137) * t231) * t181 + (t145 * t217 + (-t190 * t133 - (-t132 - t238) * t257 - (t218 * t137 - t183) * t254) * t146) * t182) * t141) * t189;
	t1 = [-t199 * t157 + t167 * t206, 0, t132, t132, t130, 0; t169 * t216 + (-t149 * t145 + (-t133 * t169 - t136 * t151) * t146) * t141 + (t136 * t215 + (0.2e1 * t136 * t263 + (-t151 * t157 + t151 - (t137 * t157 * t229 + t233) * t167) * t146 * t153 + (-(t169 * t206 - t137) * t260 + (-(t137 + t230) * t167 + t199 * t169) * t146 * t157) * t154) * t141) * t167, 0, t128, t128, (t134 * t260 - t145 * t168) * t235 + (-t134 * t259 + t152 * t145 + (t134 * t232 - t168 * t146) * t133 - (-t182 * t237 - t181 * t243 + t130 * t169 + t140 * t149 + (t140 * t247 + t170) * t137) * t146 * t254 - (-t150 + (-t130 * t189 - t137 * t191) * t182 + t200 * t140) * t231) * t141, 0; t203 * t182 * t234 + (t203 * t248 + ((-qJD(1) * t163 + 0.2e1 * t170 * t258) * t190 + (t150 * t190 - t265 * t152 - t170 * t217) * t164) * t182) * t155, 0, t129, t129, t164 * t213 * t262 + (t213 * t258 + (t151 * t246 + (-t181 * t244 + t182 * t217) * t167) * t164) * t155, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end