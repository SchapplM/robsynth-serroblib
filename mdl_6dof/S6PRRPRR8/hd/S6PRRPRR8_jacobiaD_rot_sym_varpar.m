% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR8
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
%   Wie in S6PRRPRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(12));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(12)) * cos(pkin(6));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (1326->63), mult. (4214->161), div. (275->12), fcn. (5416->13), ass. (0->79)
	t164 = sin(pkin(7));
	t162 = t164 ^ 2;
	t206 = 0.2e1 * t162;
	t166 = cos(pkin(12));
	t163 = sin(pkin(12));
	t170 = sin(qJ(2));
	t168 = cos(pkin(6));
	t172 = cos(qJ(2));
	t190 = t168 * t172;
	t182 = -t163 * t170 + t166 * t190;
	t165 = sin(pkin(6));
	t167 = cos(pkin(7));
	t195 = t165 * t167;
	t146 = t182 * t164 + t166 * t195;
	t196 = t164 * t172;
	t156 = -t165 * t196 + t168 * t167;
	t141 = atan2(t146, t156);
	t136 = sin(t141);
	t137 = cos(t141);
	t123 = t136 * t146 + t137 * t156;
	t120 = 0.1e1 / t123;
	t169 = sin(qJ(3));
	t171 = cos(qJ(3));
	t191 = t168 * t170;
	t180 = t163 * t191 - t166 * t172;
	t181 = t163 * t190 + t166 * t170;
	t187 = t163 * t164 * t165;
	t183 = -t167 * t181 + t187;
	t135 = t183 * t169 - t171 * t180;
	t131 = 0.1e1 / t135;
	t153 = 0.1e1 / t156;
	t121 = 0.1e1 / t123 ^ 2;
	t132 = 0.1e1 / t135 ^ 2;
	t154 = 0.1e1 / t156 ^ 2;
	t147 = t163 * t195 + t164 * t181;
	t145 = t147 ^ 2;
	t119 = t121 * t145 + 0.1e1;
	t152 = t180 * qJD(2);
	t157 = -t163 * t172 - t166 * t191;
	t150 = t157 * qJD(2);
	t194 = t165 * t170;
	t198 = t146 * t154;
	t185 = t194 * t198;
	t144 = t146 ^ 2;
	t140 = t144 * t154 + 0.1e1;
	t138 = 0.1e1 / t140;
	t199 = t138 * t164;
	t115 = (-qJD(2) * t185 + t150 * t153) * t199;
	t184 = -t136 * t156 + t137 * t146;
	t189 = qJD(2) * t165;
	t186 = t170 * t189;
	t112 = (t136 * t150 + t137 * t186) * t164 + t184 * t115;
	t204 = t112 * t120 * t121;
	t205 = (-t121 * t147 * t152 * t164 - t145 * t204) / t119 ^ 2;
	t192 = t167 * t171;
	t197 = t180 * t169;
	t134 = -t171 * t187 + t181 * t192 - t197;
	t130 = t134 ^ 2;
	t127 = t130 * t132 + 0.1e1;
	t151 = t181 * qJD(2);
	t193 = t167 * t169;
	t129 = t152 * t193 - t151 * t171 + (t183 * t171 + t197) * qJD(3);
	t201 = t129 * t131 * t132;
	t128 = t135 * qJD(3) - t151 * t169 - t152 * t192;
	t202 = t128 * t132;
	t203 = (-t130 * t201 + t134 * t202) / t127 ^ 2;
	t143 = -t171 * t181 + t180 * t193;
	t200 = t134 * t143;
	t188 = t154 * t162 * t170;
	t142 = -t169 * t181 - t180 * t192;
	t179 = -t153 * t157 + t185;
	t155 = t153 * t154;
	t149 = t182 * qJD(2);
	t125 = 0.1e1 / t127;
	t117 = 0.1e1 / t119;
	t116 = t179 * t199;
	t113 = (t136 * t157 + t137 * t194) * t164 - t184 * t116;
	t111 = t179 / t140 ^ 2 * (-t144 * t155 * t186 + t150 * t198) * t206 + (-t149 * t153 * t164 + (-t150 * t188 + (-t157 * t188 + (t155 * t165 * t170 ^ 2 * t206 - t154 * t196) * t146) * qJD(2)) * t165) * t138;
	t1 = [0, t111, 0, 0, 0, 0; 0, (-(t123 * t116 * t115 + t184 * t111) * t121 * t117 + 0.2e1 * (t117 * t204 + t121 * t205) * t113) * t147 + (0.2e1 * t180 * t120 * t205 + (-t151 * t120 + (t180 * t112 + t113 * t152 + (-(t115 * t157 - t116 * t150 + t172 * t189) * t137 - (-t149 + (qJD(2) * t116 - t115) * t194) * t136) * t147) * t121) * t117) * t164, 0, 0, 0, 0; 0, 0.2e1 * (-t131 * t142 + t132 * t200) * t203 + ((t143 * qJD(3) - t151 * t192 + t152 * t169) * t131 + 0.2e1 * t200 * t201 + (-t142 * t129 - (-t142 * qJD(3) + t151 * t193 + t152 * t171) * t134 - t143 * t128) * t132) * t125, -0.2e1 * t203 + 0.2e1 * (t125 * t202 + (-t125 * t201 - t132 * t203) * t134) * t134, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:11
	% DurationCPUTime: 1.44s
	% Computational Cost: add. (4346->110), mult. (14031->217), div. (500->12), fcn. (17771->13), ass. (0->108)
	t216 = sin(pkin(7));
	t218 = cos(pkin(7));
	t220 = sin(qJ(2));
	t222 = cos(qJ(2));
	t279 = sin(pkin(12));
	t281 = cos(pkin(6));
	t249 = t281 * t279;
	t280 = cos(pkin(12));
	t235 = t280 * t220 + t222 * t249;
	t217 = sin(pkin(6));
	t257 = t217 * t279;
	t199 = t216 * t235 + t218 * t257;
	t195 = 0.1e1 / t199 ^ 2;
	t236 = t220 * t249 - t280 * t222;
	t207 = t236 * qJD(2);
	t285 = t207 * t195 * t216;
	t250 = t281 * t280;
	t238 = -t279 * t220 + t222 * t250;
	t205 = t238 * qJD(2);
	t253 = t216 * t217 * t280;
	t284 = -qJD(3) * t253 + t205;
	t219 = sin(qJ(3));
	t237 = -t220 * t250 - t279 * t222;
	t232 = qJD(2) * t237;
	t221 = cos(qJ(3));
	t234 = t221 * t238;
	t283 = qJD(3) * t234 + t219 * t232;
	t233 = t238 * t219;
	t229 = t218 * t233 - t221 * t237;
	t267 = t218 * t221;
	t163 = t229 * qJD(3) + t284 * t219 - t232 * t267;
	t270 = t237 * t219;
	t182 = -t218 * t234 + t221 * t253 - t270;
	t179 = t182 ^ 2;
	t263 = t221 * t222;
	t266 = t219 * t220;
	t243 = t218 * t263 - t266;
	t258 = t216 * t281;
	t197 = -t243 * t217 - t221 * t258;
	t192 = 0.1e1 / t197 ^ 2;
	t174 = t179 * t192 + 0.1e1;
	t273 = t182 * t192;
	t264 = t220 * t221;
	t265 = t219 * t222;
	t241 = t218 * t265 + t264;
	t242 = t218 * t264 + t265;
	t251 = qJD(3) * t258;
	t177 = t219 * t251 + (t242 * qJD(2) + t241 * qJD(3)) * t217;
	t191 = 0.1e1 / t197;
	t274 = t177 * t191 * t192;
	t282 = -0.2e1 * (t163 * t273 - t179 * t274) / t174 ^ 2;
	t175 = atan2(-t182, t197);
	t168 = sin(t175);
	t169 = cos(t175);
	t162 = -t168 * t182 + t169 * t197;
	t159 = 0.1e1 / t162;
	t194 = 0.1e1 / t199;
	t160 = 0.1e1 / t162 ^ 2;
	t170 = 0.1e1 / t174;
	t151 = (-t163 * t191 + t177 * t273) * t170;
	t248 = -t168 * t197 - t169 * t182;
	t148 = t248 * t151 - t163 * t168 + t169 * t177;
	t278 = t148 * t159 * t160;
	t252 = t216 * t257;
	t269 = t236 * t219;
	t185 = -t221 * t252 + t235 * t267 - t269;
	t277 = t160 * t185;
	t276 = t168 * t185;
	t275 = t169 * t185;
	t239 = -t218 * t235 + t252;
	t186 = t239 * t219 - t221 * t236;
	t272 = t186 * t195;
	t268 = t218 * t219;
	t180 = t185 ^ 2;
	t157 = t160 * t180 + 0.1e1;
	t206 = t235 * qJD(2);
	t165 = t186 * qJD(3) - t206 * t219 - t207 * t267;
	t262 = 0.2e1 * (t165 * t277 - t180 * t278) / t157 ^ 2;
	t166 = t207 * t268 - t206 * t221 + (t239 * t221 + t269) * qJD(3);
	t181 = t186 ^ 2;
	t176 = t181 * t195 + 0.1e1;
	t260 = t194 * t285;
	t261 = 0.2e1 * (t166 * t272 + t181 * t260) / t176 ^ 2;
	t259 = qJD(3) * t270;
	t255 = -0.2e1 * t182 * t274;
	t254 = 0.2e1 * t185 * t278;
	t184 = -t219 * t253 + t229;
	t198 = t241 * t217 + t219 * t258;
	t245 = -t184 * t191 + t198 * t273;
	t188 = -t237 * t267 + t233;
	t204 = t242 * t217;
	t244 = -t188 * t191 + t204 * t273;
	t189 = -t219 * t235 - t236 * t267;
	t190 = -t221 * t235 + t236 * t268;
	t240 = -t218 * t266 + t263;
	t187 = (t243 * qJD(2) + t240 * qJD(3)) * t217;
	t178 = t221 * t251 + (t240 * qJD(2) + t243 * qJD(3)) * t217;
	t172 = 0.1e1 / t176;
	t167 = t205 * t267 + t218 * t259 + t283;
	t164 = t283 * t218 + t284 * t221 + t259;
	t155 = 0.1e1 / t157;
	t153 = t244 * t170;
	t152 = t245 * t170;
	t150 = t248 * t153 - t168 * t188 + t169 * t204;
	t149 = t248 * t152 - t168 * t184 + t169 * t198;
	t147 = t244 * t282 + (t204 * t255 - t167 * t191 + (t163 * t204 + t177 * t188 + t182 * t187) * t192) * t170;
	t146 = t245 * t282 + (t198 * t255 - t164 * t191 + (t163 * t198 + t177 * t184 + t178 * t182) * t192) * t170;
	t1 = [0, t147, t146, 0, 0, 0; 0, (t150 * t277 - t159 * t189) * t262 + ((t190 * qJD(3) - t206 * t267 + t207 * t219) * t159 + t150 * t254 + (-t189 * t148 - t150 * t165 - (-t147 * t182 - t153 * t163 + t187 + (-t153 * t197 - t188) * t151) * t275 - (-t147 * t197 - t153 * t177 - t167 + (t153 * t182 - t204) * t151) * t276) * t160) * t155, (t149 * t277 - t159 * t186) * t262 + (t149 * t254 + t166 * t159 + (-t186 * t148 - t149 * t165 - (-t146 * t182 - t152 * t163 + t178 + (-t152 * t197 - t184) * t151) * t275 - (-t146 * t197 - t152 * t177 - t164 + (t152 * t182 - t198) * t151) * t276) * t160) * t155, 0, 0, 0; 0, (-t216 * t236 * t272 - t190 * t194) * t261 + ((-t189 * qJD(3) + t206 * t268 + t207 * t221) * t194 + (0.2e1 * t186 * t236 * t260 + (t166 * t236 + t186 * t206 + t190 * t207) * t195) * t216) * t172, t185 * t194 * t261 + (-t165 * t194 - t185 * t285) * t172, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:12
	% DurationCPUTime: 1.74s
	% Computational Cost: add. (5341->132), mult. (16984->251), div. (538->12), fcn. (21480->15), ass. (0->120)
	t286 = cos(pkin(7));
	t289 = sin(qJ(3));
	t334 = t286 * t289;
	t288 = sin(qJ(5));
	t291 = cos(qJ(5));
	t283 = sin(pkin(12));
	t284 = sin(pkin(7));
	t285 = sin(pkin(6));
	t287 = cos(pkin(6));
	t290 = sin(qJ(2));
	t354 = cos(pkin(12));
	t323 = t354 * t290;
	t293 = cos(qJ(2));
	t339 = t283 * t293;
	t302 = t287 * t339 + t323;
	t309 = t283 * t285 * t286 + t284 * t302;
	t292 = cos(qJ(3));
	t338 = t284 * t285;
	t325 = t283 * t338;
	t333 = t286 * t292;
	t322 = t354 * t293;
	t340 = t283 * t290;
	t303 = t287 * t340 - t322;
	t341 = t303 * t289;
	t316 = -t292 * t325 + t302 * t333 - t341;
	t237 = t316 * t288 + t309 * t291;
	t301 = -t287 * t323 - t339;
	t318 = t354 * t338;
	t356 = -t289 * t301 + t292 * t318;
	t274 = t287 * t322 - t340;
	t270 = t274 * qJD(2);
	t271 = t301 * qJD(2);
	t344 = t274 * t292;
	t227 = t271 * t334 + t270 * t292 + (t344 * t286 - t356) * qJD(3);
	t342 = t301 * t292;
	t253 = t274 * t334 - t289 * t318 - t342;
	t250 = t253 ^ 2;
	t330 = t290 * t292;
	t331 = t289 * t293;
	t305 = t286 * t331 + t330;
	t337 = t284 * t287;
	t265 = t305 * t285 + t289 * t337;
	t262 = 0.1e1 / t265 ^ 2;
	t242 = t250 * t262 + 0.1e1;
	t345 = t253 * t262;
	t329 = t292 * t293;
	t332 = t289 * t290;
	t304 = -t286 * t332 + t329;
	t307 = t286 * t329 - t332;
	t324 = qJD(3) * t337;
	t248 = t292 * t324 + (t304 * qJD(2) + t307 * qJD(3)) * t285;
	t261 = 0.1e1 / t265;
	t346 = t248 * t261 * t262;
	t355 = -0.2e1 * (t227 * t345 - t250 * t346) / t242 ^ 2;
	t243 = atan2(-t253, t265);
	t238 = sin(t243);
	t239 = cos(t243);
	t220 = -t238 * t253 + t239 * t265;
	t217 = 0.1e1 / t220;
	t233 = 0.1e1 / t237;
	t218 = 0.1e1 / t220 ^ 2;
	t234 = 0.1e1 / t237 ^ 2;
	t240 = 0.1e1 / t242;
	t210 = (-t227 * t261 + t248 * t345) * t240;
	t315 = -t238 * t265 - t239 * t253;
	t206 = t315 * t210 - t238 * t227 + t239 * t248;
	t353 = t206 * t217 * t218;
	t308 = -t286 * t302 + t325;
	t255 = t308 * t289 - t292 * t303;
	t272 = t302 * qJD(2);
	t273 = t303 * qJD(2);
	t228 = t255 * qJD(3) - t272 * t289 - t273 * t333;
	t336 = t284 * t288;
	t221 = -t237 * qJD(5) + t228 * t291 + t273 * t336;
	t236 = t309 * t288 - t316 * t291;
	t232 = t236 ^ 2;
	t225 = t232 * t234 + 0.1e1;
	t349 = t234 * t236;
	t328 = qJD(5) * t236;
	t335 = t284 * t291;
	t222 = t228 * t288 - t273 * t335 - t328;
	t350 = t222 * t233 * t234;
	t352 = (-t221 * t349 - t232 * t350) / t225 ^ 2;
	t351 = t218 * t255;
	t348 = t238 * t255;
	t347 = t239 * t255;
	t251 = t255 ^ 2;
	t216 = t218 * t251 + 0.1e1;
	t229 = t273 * t334 - t272 * t292 + (t308 * t292 + t341) * qJD(3);
	t327 = 0.2e1 * (t229 * t351 - t251 * t353) / t216 ^ 2;
	t326 = 0.2e1 * t352;
	t321 = 0.2e1 * t236 * t350;
	t320 = -0.2e1 * t253 * t346;
	t319 = 0.2e1 * t255 * t353;
	t312 = t291 * t233 + t288 * t349;
	t252 = -t274 * t333 + t356;
	t264 = t307 * t285 + t292 * t337;
	t311 = t252 * t261 + t264 * t345;
	t257 = t301 * t334 + t344;
	t268 = t304 * t285;
	t310 = -t257 * t261 + t268 * t345;
	t258 = -t289 * t302 - t303 * t333;
	t245 = t258 * t288 - t303 * t335;
	t244 = -t258 * t291 - t303 * t336;
	t259 = -t292 * t302 + t303 * t334;
	t306 = -t286 * t330 - t331;
	t256 = (-t305 * qJD(2) + t306 * qJD(3)) * t285;
	t247 = -t289 * t324 + (t306 * qJD(2) - t305 * qJD(3)) * t285;
	t231 = t259 * qJD(3) - t272 * t333 + t273 * t289;
	t230 = -t270 * t334 + t271 * t292 + (-t274 * t289 + t301 * t333) * qJD(3);
	t226 = -t271 * t333 + t270 * t289 + (-t342 + (t274 * t286 - t318) * t289) * qJD(3);
	t223 = 0.1e1 / t225;
	t214 = 0.1e1 / t216;
	t212 = t310 * t240;
	t211 = t311 * t240;
	t208 = t315 * t212 - t238 * t257 + t239 * t268;
	t207 = t315 * t211 + t238 * t252 + t239 * t264;
	t205 = t310 * t355 + (t268 * t320 - t230 * t261 + (t227 * t268 + t248 * t257 + t253 * t256) * t262) * t240;
	t204 = t311 * t355 + (t264 * t320 + t226 * t261 + (t227 * t264 + t247 * t253 - t248 * t252) * t262) * t240;
	t1 = [0, t205, t204, 0, 0, 0; 0, (t208 * t351 - t217 * t259) * t327 + ((-t258 * qJD(3) + t272 * t334 + t273 * t292) * t217 + t208 * t319 + (-t259 * t206 - t208 * t229 - (-t205 * t253 - t212 * t227 + t256 + (-t212 * t265 - t257) * t210) * t347 - (-t205 * t265 - t212 * t248 - t230 + (t212 * t253 - t268) * t210) * t348) * t218) * t214, (t207 * t351 + t316 * t217) * t327 + (t207 * t319 - t228 * t217 + (t316 * t206 - t207 * t229 - (-t204 * t253 - t211 * t227 + t247 + (-t211 * t265 + t252) * t210) * t347 - (-t204 * t265 - t211 * t248 + t226 + (t211 * t253 - t264) * t210) * t348) * t218) * t214, 0, 0, 0; 0, (-t233 * t244 + t245 * t349) * t326 + ((t245 * qJD(5) - t231 * t291 - t272 * t336) * t233 + t245 * t321 + (-t244 * t222 - (-t244 * qJD(5) + t231 * t288 - t272 * t335) * t236 + t245 * t221) * t234) * t223, t312 * t255 * t326 + (-t312 * t229 + ((qJD(5) * t233 + t321) * t288 + (t221 * t288 + (t222 - t328) * t291) * t234) * t255) * t223, 0, -0.2e1 * t352 + 0.2e1 * (-t221 * t234 * t223 + (-t223 * t350 - t234 * t352) * t236) * t236, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:11
	% EndTime: 2019-10-09 22:39:15
	% DurationCPUTime: 4.25s
	% Computational Cost: add. (14691->203), mult. (44088->385), div. (816->12), fcn. (56595->17), ass. (0->165)
	t395 = sin(qJ(6));
	t483 = sin(pkin(12));
	t485 = cos(pkin(6));
	t434 = t485 * t483;
	t484 = cos(pkin(12));
	t486 = sin(qJ(2));
	t488 = cos(qJ(2));
	t388 = -t486 * t434 + t484 * t488;
	t394 = cos(pkin(7));
	t397 = sin(qJ(3));
	t416 = t488 * t434 + t484 * t486;
	t392 = sin(pkin(7));
	t393 = sin(pkin(6));
	t449 = t393 * t483;
	t439 = t392 * t449;
	t487 = cos(qJ(3));
	t368 = t388 * t487 + (-t394 * t416 + t439) * t397;
	t398 = cos(qJ(6));
	t472 = t368 * t398;
	t396 = sin(qJ(5));
	t399 = cos(qJ(5));
	t422 = -t392 * t416 - t394 * t449;
	t454 = t394 * t487;
	t489 = -t388 * t397 - t416 * t454 + t487 * t439;
	t491 = t396 * t489 + t422 * t399;
	t326 = -t395 * t491 - t472;
	t495 = 0.2e1 * t326;
	t435 = t485 * t484;
	t414 = -t488 * t435 + t483 * t486;
	t415 = t486 * t435 + t483 * t488;
	t450 = t393 * t484;
	t440 = t392 * t450;
	t367 = t415 * t487 + (-t414 * t394 - t440) * t397;
	t382 = t414 * qJD(2);
	t383 = t415 * qJD(2);
	t340 = t367 * qJD(3) - t382 * t397 + t383 * t454;
	t410 = t414 * t487;
	t407 = t394 * t410 + t415 * t397 + t487 * t440;
	t409 = t414 * t392 - t394 * t450;
	t351 = t407 * t396 + t409 * t399;
	t468 = t392 * t396;
	t316 = -t351 * qJD(5) + t340 * t399 - t383 * t468;
	t349 = t409 * t396 - t407 * t399;
	t347 = t349 ^ 2;
	t443 = t488 * t487;
	t452 = t486 * t397;
	t419 = t394 * t443 - t452;
	t451 = t392 * t485;
	t430 = t487 * t451;
	t377 = -t419 * t393 - t430;
	t469 = t392 * t393;
	t444 = t488 * t469;
	t386 = t485 * t394 - t444;
	t432 = t377 * t399 - t386 * t396;
	t365 = 0.1e1 / t432 ^ 2;
	t332 = t347 * t365 + 0.1e1;
	t330 = 0.1e1 / t332;
	t442 = t486 * t487;
	t453 = t488 * t397;
	t417 = t394 * t442 + t453;
	t420 = t394 * t453 + t442;
	t441 = t397 * t451;
	t359 = qJD(3) * t441 + (t417 * qJD(2) + t420 * qJD(3)) * t393;
	t370 = t377 * t396 + t386 * t399;
	t445 = t486 * t469;
	t431 = qJD(2) * t445;
	t334 = t370 * qJD(5) - t359 * t399 + t396 * t431;
	t364 = 0.1e1 / t432;
	t474 = t349 * t365;
	t299 = (-t316 * t364 + t334 * t474) * t330;
	t333 = atan2(-t349, -t432);
	t328 = sin(t333);
	t329 = cos(t333);
	t433 = t328 * t432 - t329 * t349;
	t294 = t433 * t299 + t328 * t316 + t329 * t334;
	t312 = -t328 * t349 - t329 * t432;
	t309 = 0.1e1 / t312;
	t310 = 0.1e1 / t312 ^ 2;
	t494 = t294 * t309 * t310;
	t352 = -t422 * t396 + t399 * t489;
	t448 = 0.2e1 * t352 * t494;
	t475 = t334 * t364 * t365;
	t493 = (-t316 * t474 + t347 * t475) / t332 ^ 2;
	t378 = t420 * t393 + t441;
	t425 = t364 * t367 + t378 * t474;
	t492 = t399 * t425;
	t490 = -0.2e1 * t493;
	t327 = t368 * t395 - t398 * t491;
	t323 = 0.1e1 / t327;
	t324 = 0.1e1 / t327 ^ 2;
	t384 = t416 * qJD(2);
	t385 = t388 * qJD(2);
	t342 = t368 * qJD(3) - t384 * t397 + t385 * t454;
	t467 = t392 * t399;
	t319 = -t352 * qJD(5) + t342 * t396 + t385 * t467;
	t466 = t394 * t397;
	t343 = t489 * qJD(3) - t384 * t487 - t385 * t466;
	t463 = qJD(6) * t326;
	t308 = t319 * t398 + t343 * t395 - t463;
	t482 = t308 * t323 * t324;
	t481 = t310 * t352;
	t322 = t326 ^ 2;
	t315 = t322 * t324 + 0.1e1;
	t313 = 0.1e1 / t315;
	t480 = t313 * t324;
	t307 = t327 * qJD(6) + t319 * t395 - t343 * t398;
	t478 = t324 * t326;
	t479 = 0.1e1 / t315 ^ 2 * (t307 * t478 - t322 * t482);
	t477 = t328 * t352;
	t476 = t329 * t352;
	t473 = t368 * t396;
	t471 = t368 * t399;
	t465 = qJD(5) * t396;
	t464 = qJD(5) * t399;
	t348 = t352 ^ 2;
	t306 = t310 * t348 + 0.1e1;
	t318 = t491 * qJD(5) + t342 * t399 - t385 * t468;
	t462 = 0.2e1 * (-t318 * t481 - t348 * t494) / t306 ^ 2;
	t460 = -0.2e1 * t479;
	t459 = 0.2e1 * t479;
	t458 = t324 * t479;
	t457 = t307 * t480;
	t456 = t326 * t482;
	t455 = t349 * t475;
	t447 = 0.2e1 * t456;
	t446 = 0.2e1 * t455;
	t436 = qJD(6) * t473 + t342;
	t423 = -t388 * t454 + t397 * t416;
	t356 = t388 * t467 - t396 * t423;
	t373 = -t388 * t466 - t416 * t487;
	t339 = t356 * t398 + t373 * t395;
	t338 = t356 * t395 - t373 * t398;
	t428 = t323 * t395 - t398 * t478;
	t427 = t351 * t364 + t370 * t474;
	t411 = t394 * t415;
	t371 = -t414 * t397 + t487 * t411;
	t412 = t392 * t415;
	t354 = t371 * t399 - t396 * t412;
	t381 = t417 * t393;
	t374 = -t381 * t399 + t396 * t445;
	t426 = -t354 * t364 + t374 * t474;
	t355 = t388 * t468 + t399 * t423;
	t418 = -t394 * t452 + t443;
	t413 = qJD(6) * t489 + t343 * t396 + t368 * t464;
	t360 = qJD(3) * t430 + (t418 * qJD(2) + t419 * qJD(3)) * t393;
	t346 = t423 * qJD(3) + t384 * t466 - t385 * t487;
	t345 = t373 * qJD(3) - t384 * t454 - t385 * t397;
	t344 = qJD(2) * t396 * t444 + t445 * t464 - (t419 * qJD(2) + t418 * qJD(3)) * t393 * t399 + t381 * t465;
	t341 = -t407 * qJD(3) - t382 * t487 - t383 * t466;
	t337 = t395 * t489 + t396 * t472;
	t335 = t432 * qJD(5) + t359 * t396 + t399 * t431;
	t321 = -t355 * qJD(5) + t345 * t396 - t384 * t467;
	t320 = t382 * t468 - t412 * t464 + (-t382 * t454 - t383 * t397 + (-t397 * t411 - t410) * qJD(3)) * t399 - t371 * t465;
	t317 = t349 * qJD(5) - t340 * t396 - t383 * t467;
	t304 = 0.1e1 / t306;
	t303 = t330 * t492;
	t302 = t426 * t330;
	t301 = t427 * t330;
	t297 = (t328 * t367 - t329 * t378) * t399 - t433 * t303;
	t296 = t433 * t302 + t328 * t354 + t329 * t374;
	t295 = t433 * t301 - t328 * t351 + t329 * t370;
	t293 = t426 * t490 + (t374 * t446 - t320 * t364 + (-t316 * t374 - t334 * t354 + t344 * t349) * t365) * t330;
	t291 = t427 * t490 + (t370 * t446 - t317 * t364 + (-t316 * t370 + t334 * t351 + t335 * t349) * t365) * t330;
	t290 = 0.2e1 * t492 * t493 + (t425 * t465 + (-0.2e1 * t378 * t455 - t341 * t364 + (t316 * t378 - t334 * t367 - t349 * t360) * t365) * t399) * t330;
	t1 = [0, t293, t290, 0, t291, 0; 0, (t296 * t481 - t309 * t355) * t462 + ((t356 * qJD(5) - t345 * t399 - t384 * t468) * t309 + t296 * t448 + (-t355 * t294 + t296 * t318 - (-t293 * t349 + t302 * t316 + t344 + (t302 * t432 + t354) * t299) * t476 - (t293 * t432 - t302 * t334 + t320 + (t302 * t349 - t374) * t299) * t477) * t310) * t304, (t297 * t481 + t309 * t471) * t462 + ((-t343 * t399 + t368 * t465) * t309 + t297 * t448 + (t297 * t318 + t471 * t294 - (t378 * t465 - t290 * t349 - t303 * t316 - t360 * t399 + (-t303 * t432 + t367 * t399) * t299) * t476 - (-t367 * t465 + t290 * t432 + t303 * t334 + t341 * t399 + (-t303 * t349 + t378 * t399) * t299) * t477) * t310) * t304, 0, (t295 * t481 + t309 * t491) * t462 + (t295 * t448 + t319 * t309 + (t491 * t294 + t295 * t318 - (-t291 * t349 + t301 * t316 + t335 + (t301 * t432 - t351) * t299) * t476 - (t291 * t432 - t301 * t334 + t317 + (t301 * t349 - t370) * t299) * t477) * t310) * t304, 0; 0, (-t323 * t338 + t339 * t478) * t459 + ((t339 * qJD(6) + t321 * t395 - t346 * t398) * t323 + t339 * t447 + (-t338 * t308 - (-t338 * qJD(6) + t321 * t398 + t346 * t395) * t326 - t339 * t307) * t324) * t313, (t458 * t495 - t457) * t337 + (-t308 * t480 + t323 * t460) * (t395 * t473 - t398 * t489) + ((t413 * t395 + t436 * t398) * t323 - (-t436 * t395 + t413 * t398) * t478 + t337 * t447) * t313, 0, t428 * t352 * t459 + (t428 * t318 + ((-qJD(6) * t323 - 0.2e1 * t456) * t398 + (t307 * t398 + (t308 - t463) * t395) * t324) * t352) * t313, t460 + (t457 + (-t313 * t482 - t458) * t326) * t495;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end