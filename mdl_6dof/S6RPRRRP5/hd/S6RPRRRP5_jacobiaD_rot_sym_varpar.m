% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP5
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
%   Wie in S6RPRRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:24
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (5587->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t174 = sin(qJ(1));
	t236 = 0.2e1 * t174;
	t171 = t174 ^ 2;
	t169 = pkin(10) + qJ(3) + qJ(4);
	t166 = sin(t169);
	t162 = t166 ^ 2;
	t167 = cos(t169);
	t164 = 0.1e1 / t167 ^ 2;
	t221 = t162 * t164;
	t158 = t171 * t221 + 0.1e1;
	t156 = 0.1e1 / t158;
	t163 = 0.1e1 / t167;
	t176 = cos(qJ(1));
	t207 = qJD(1) * t176;
	t197 = t166 * t207;
	t170 = qJD(3) + qJD(4);
	t215 = t170 * t174;
	t200 = t164 * t215;
	t130 = (-(-t167 * t215 - t197) * t163 + t162 * t200) * t156;
	t235 = t130 - t215;
	t175 = cos(qJ(5));
	t209 = t175 * t176;
	t173 = sin(qJ(5));
	t211 = t174 * t173;
	t155 = t167 * t209 + t211;
	t212 = t174 * t166;
	t151 = atan2(-t212, -t167);
	t146 = cos(t151);
	t145 = sin(t151);
	t202 = t145 * t212;
	t140 = -t146 * t167 - t202;
	t137 = 0.1e1 / t140;
	t148 = 0.1e1 / t155;
	t138 = 0.1e1 / t140 ^ 2;
	t149 = 0.1e1 / t155 ^ 2;
	t234 = t156 - 0.1e1;
	t226 = t146 * t166;
	t125 = (-t130 * t174 + t170) * t226 + (t167 * t235 - t197) * t145;
	t233 = t125 * t137 * t138;
	t185 = t167 * t211 + t209;
	t214 = t170 * t176;
	t198 = t166 * t214;
	t135 = t185 * qJD(1) - qJD(5) * t155 + t173 * t198;
	t210 = t174 * t175;
	t213 = t173 * t176;
	t154 = t167 * t213 - t210;
	t147 = t154 ^ 2;
	t144 = t147 * t149 + 0.1e1;
	t224 = t149 * t154;
	t191 = -qJD(1) * t167 + qJD(5);
	t192 = qJD(5) * t167 - qJD(1);
	t136 = -t192 * t213 + (t191 * t174 - t198) * t175;
	t230 = t136 * t148 * t149;
	t232 = (-t135 * t224 - t147 * t230) / t144 ^ 2;
	t161 = t166 * t162;
	t218 = t163 * t166;
	t184 = t170 * (t161 * t163 * t164 + t218);
	t219 = t162 * t174;
	t189 = t207 * t219;
	t231 = (t164 * t189 + t171 * t184) / t158 ^ 2;
	t229 = t138 * t166;
	t228 = t138 * t176;
	t227 = t145 * t174;
	t225 = t148 * t173;
	t223 = t154 * t175;
	t222 = t162 * t163;
	t172 = t176 ^ 2;
	t220 = t162 * t172;
	t217 = t166 * t176;
	t216 = t167 * t170;
	t208 = qJD(1) * t174;
	t133 = t138 * t220 + 0.1e1;
	t206 = 0.2e1 * (-t220 * t233 + (t166 * t172 * t216 - t189) * t138) / t133 ^ 2;
	t205 = 0.2e1 * t233;
	t204 = -0.2e1 * t232;
	t203 = t138 * t217;
	t201 = t154 * t230;
	t196 = 0.1e1 + t221;
	t195 = t166 * t206;
	t194 = -0.2e1 * t166 * t231;
	t193 = t231 * t236;
	t190 = t146 * t156 * t222;
	t188 = t196 * t176;
	t187 = t191 * t176;
	t186 = t149 * t223 - t225;
	t153 = -t167 * t210 + t213;
	t142 = 0.1e1 / t144;
	t141 = t196 * t174 * t156;
	t131 = 0.1e1 / t133;
	t129 = (t234 * t166 * t145 - t174 * t190) * t176;
	t127 = -t167 * t227 + t226 + (t145 * t167 - t146 * t212) * t141;
	t126 = -t196 * t193 + (qJD(1) * t188 + t184 * t236) * t156;
	t123 = t186 * t204 * t217 + (t186 * t167 * t214 + (-t186 * t208 + ((-qJD(5) * t148 - 0.2e1 * t201) * t175 + (-t135 * t175 + (-qJD(5) * t154 + t136) * t173) * t149) * t176) * t166) * t142;
	t122 = (t127 * t229 - t137 * t167) * t176 * t206 + ((-t137 * t208 + (-t127 * t170 - t125) * t228) * t167 + (-t137 * t214 - (-t126 * t146 * t174 - t235 * t145 + (t130 * t227 - t145 * t170 - t146 * t207) * t141) * t203 + (t138 * t208 + t176 * t205) * t127 - ((t126 - t207) * t145 + ((-t141 * t174 + 0.1e1) * t170 + (t141 - t174) * t130) * t146) * t167 * t228) * t166) * t131;
	t1 = [t163 * t176 * t194 + (t170 * t188 - t208 * t218) * t156, 0, t126, t126, 0, 0; (t137 * t195 + (-t137 * t216 + (qJD(1) * t129 + t125) * t229) * t131) * t174 + (t138 * t195 * t129 + (-((t194 - t216 + (t130 * t163 * t219 + t216) * t156) * t145 + (t193 * t222 - t130 * t166 + (-t161 * t200 + (t130 - 0.2e1 * t215) * t166) * t156) * t146) * t203 + (-t138 * t216 + t166 * t205) * t129 + (-t137 + ((-t171 + t172) * t190 + t234 * t202) * t138) * t166 * qJD(1)) * t131) * t176, 0, t122, t122, 0, 0; 0.2e1 * (t148 * t185 + t153 * t224) * t232 + (0.2e1 * t153 * t201 - t192 * t148 * t210 + (t170 * t212 + t187) * t225 + (t153 * t135 + t185 * t136 - t187 * t223 - (t166 * t170 * t175 + t192 * t173) * t154 * t174) * t149) * t142, 0, t123, t123, t204 + 0.2e1 * (-t135 * t142 * t149 + (-t142 * t230 - t149 * t232) * t154) * t154, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:52:23
	% EndTime: 2019-10-10 01:52:25
	% DurationCPUTime: 1.63s
	% Computational Cost: add. (8326->125), mult. (8382->273), div. (1515->15), fcn. (10508->9), ass. (0->119)
	t197 = pkin(10) + qJ(3) + qJ(4);
	t194 = cos(t197);
	t203 = sin(qJ(5));
	t278 = sin(qJ(1));
	t234 = t278 * t203;
	t204 = cos(qJ(5));
	t205 = cos(qJ(1));
	t255 = t205 * t204;
	t180 = t194 * t255 + t234;
	t174 = 0.1e1 / t180 ^ 2;
	t193 = sin(t197);
	t189 = t193 ^ 2;
	t202 = t205 ^ 2;
	t266 = t189 * t202;
	t242 = t174 * t266;
	t169 = 0.1e1 + t242;
	t229 = qJD(1) * t278;
	t198 = qJD(3) + qJD(4);
	t259 = t198 * t205;
	t237 = t193 * t259;
	t215 = t194 * t229 + t237;
	t228 = t278 * qJD(5);
	t256 = t205 * t203;
	t159 = (-qJD(5) * t194 + qJD(1)) * t256 + (t228 - t215) * t204;
	t173 = 0.1e1 / t180;
	t273 = t159 * t173 * t174;
	t223 = t266 * t273;
	t238 = t193 * t198 * t202;
	t281 = (-t223 + (-t189 * t205 * t229 + t194 * t238) * t174) / t169 ^ 2;
	t262 = t193 * t205;
	t176 = t194 * t234 + t255;
	t220 = t203 * t228;
	t251 = qJD(5) * t205;
	t231 = t204 * t251;
	t158 = t176 * qJD(1) - t194 * t231 + t203 * t237 - t220;
	t233 = t278 * t204;
	t179 = t194 * t256 - t233;
	t190 = 0.1e1 / t193;
	t191 = 0.1e1 / t193 ^ 2;
	t200 = 0.1e1 / t203 ^ 2;
	t252 = qJD(5) * t204;
	t232 = t200 * t252;
	t199 = 0.1e1 / t203;
	t260 = t198 * t199;
	t236 = t194 * t260;
	t265 = t190 * t199;
	t280 = (t190 * t232 + t191 * t236) * t179 + t158 * t265;
	t263 = t193 * t203;
	t168 = atan2(-t176, t263);
	t163 = cos(t168);
	t162 = sin(t168);
	t272 = t162 * t176;
	t157 = t163 * t263 - t272;
	t154 = 0.1e1 / t157;
	t155 = 0.1e1 / t157 ^ 2;
	t279 = 0.2e1 * t179;
	t171 = t176 ^ 2;
	t264 = t191 * t200;
	t170 = t171 * t264 + 0.1e1;
	t166 = 0.1e1 / t170;
	t261 = t194 * t198;
	t216 = t193 * t252 + t203 * t261;
	t240 = t176 * t264;
	t221 = t204 * t229;
	t235 = t278 * t193;
	t222 = t198 * t235;
	t254 = qJD(1) * t205;
	t160 = t204 * t228 * t194 - t221 + (t254 * t194 - t222 - t251) * t203;
	t243 = t160 * t265;
	t146 = (t216 * t240 - t243) * t166;
	t213 = -t146 * t176 + t216;
	t142 = (-t146 * t263 - t160) * t162 + t213 * t163;
	t156 = t154 * t155;
	t277 = t142 * t156;
	t192 = t190 / t189;
	t201 = t199 * t200;
	t276 = (t160 * t240 + (-t191 * t201 * t252 - t192 * t200 * t261) * t171) / t170 ^ 2;
	t275 = t155 * t179;
	t274 = t158 * t155;
	t271 = t162 * t179;
	t270 = t162 * t193;
	t269 = t163 * t176;
	t268 = t163 * t179;
	t267 = t163 * t194;
	t258 = t200 * t204;
	t257 = t205 * t154;
	t253 = qJD(5) * t203;
	t172 = t179 ^ 2;
	t152 = t155 * t172 + 0.1e1;
	t250 = 0.2e1 * (-t172 * t277 - t179 * t274) / t152 ^ 2;
	t249 = 0.2e1 * t281;
	t248 = -0.2e1 * t276;
	t247 = t156 * t279;
	t246 = t190 * t276;
	t245 = t155 * t271;
	t241 = t176 * t265;
	t239 = t191 * t194 * t199;
	t218 = t176 * t239 + t278;
	t153 = t218 * t166;
	t230 = t278 - t153;
	t227 = t154 * t250;
	t226 = t155 * t250;
	t225 = t262 * t279;
	t224 = t199 * t246;
	t178 = t194 * t233 - t256;
	t219 = t176 * t258 - t178 * t199;
	t217 = t174 * t178 * t205 - t278 * t173;
	t164 = 0.1e1 / t169;
	t161 = t180 * qJD(1) - t194 * t220 - t204 * t222 - t231;
	t150 = 0.1e1 / t152;
	t149 = t219 * t190 * t166;
	t145 = (-t162 + (t163 * t241 + t162) * t166) * t179;
	t144 = -t153 * t269 + (t230 * t270 + t267) * t203;
	t143 = t163 * t193 * t204 - t162 * t178 + (-t162 * t263 - t269) * t149;
	t141 = t218 * t248 + (t160 * t239 + t254 + (-t190 * t260 + (-t191 * t232 - 0.2e1 * t192 * t236) * t194) * t176) * t166;
	t139 = -0.2e1 * t219 * t246 + (-t219 * t191 * t261 + (t160 * t258 - t161 * t199 + (t178 * t258 + (-0.2e1 * t201 * t204 ^ 2 - t199) * t176) * qJD(5)) * t190) * t166;
	t138 = (t173 * t194 * t205 + t204 * t242) * t249 + (0.2e1 * t204 * t223 + t215 * t173 + ((t159 * t205 - 0.2e1 * t204 * t238) * t194 + (t202 * t253 + 0.2e1 * t205 * t221) * t189) * t174) * t164;
	t137 = t144 * t179 * t226 + (-(-t141 * t269 + (t146 * t272 - t160 * t163) * t153) * t275 + (t142 * t247 + t274) * t144 + (-t193 * t257 - (-t153 * t270 + t162 * t235 + t267) * t275) * t252) * t150 + (t227 * t262 + ((-t198 * t257 - (t230 * t198 - t146) * t245) * t194 + (t154 * t229 + (t205 * t142 - (-t141 + t254) * t271 - (t230 * t146 - t198) * t268) * t155) * t193) * t150) * t203;
	t1 = [t166 * t280 + t224 * t279, 0, t141, t141, t139, 0; t176 * t227 + (-t160 * t154 + (t142 * t176 + t145 * t158) * t155) * t150 + (t145 * t226 + (0.2e1 * t145 * t277 + (t158 * t166 - t158 - (-t146 * t166 * t241 + t248) * t179) * t155 * t162 + (-(-0.2e1 * t176 * t224 - t146) * t275 + (-(t146 + t243) * t179 + t280 * t176) * t155 * t166) * t163) * t150) * t179, 0, t137, t137, (t143 * t275 - t154 * t180) * t250 + (t143 * t274 + t159 * t154 + (t143 * t247 - t155 * t180) * t142 - (-t193 * t253 + t204 * t261 - t139 * t176 - t149 * t160 + (-t149 * t263 - t178) * t146) * t155 * t268 - (-t161 + (-t139 * t203 - t146 * t204) * t193 - t213 * t149) * t245) * t150, 0; t217 * t193 * t249 + (-t217 * t261 + ((qJD(1) * t173 + 0.2e1 * t178 * t273) * t205 + (-t278 * t159 - t161 * t205 + t178 * t229) * t174) * t193) * t164, 0, t138, t138, t174 * t225 * t281 + (t225 * t273 + (t158 * t262 + (t193 * t229 - t194 * t259) * t179) * t174) * t164, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end