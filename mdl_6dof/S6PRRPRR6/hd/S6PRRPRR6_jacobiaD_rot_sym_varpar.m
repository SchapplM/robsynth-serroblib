% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR6
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
%   Wie in S6PRRPRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiaD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
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
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.67s
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
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:21
	% DurationCPUTime: 1.61s
	% Computational Cost: add. (4823->122), mult. (15462->242), div. (514->12), fcn. (19568->15), ass. (0->121)
	t259 = sin(qJ(2));
	t261 = cos(qJ(2));
	t323 = cos(pkin(12));
	t324 = cos(pkin(6));
	t289 = t324 * t323;
	t322 = sin(pkin(12));
	t277 = -t322 * t259 + t261 * t289;
	t242 = t277 * qJD(2);
	t254 = sin(pkin(7));
	t255 = sin(pkin(6));
	t292 = t254 * t255 * t323;
	t327 = -qJD(3) * t292 + t242;
	t258 = sin(qJ(3));
	t276 = -t259 * t289 - t322 * t261;
	t271 = qJD(2) * t276;
	t260 = cos(qJ(3));
	t273 = t260 * t277;
	t326 = qJD(3) * t273 + t258 * t271;
	t257 = cos(pkin(7));
	t272 = t277 * t258;
	t268 = t257 * t272 - t260 * t276;
	t306 = t257 * t260;
	t197 = t268 * qJD(3) + t327 * t258 - t271 * t306;
	t311 = t276 * t258;
	t221 = -t257 * t273 + t260 * t292 - t311;
	t219 = t221 ^ 2;
	t302 = t260 * t261;
	t305 = t258 * t259;
	t282 = t257 * t302 - t305;
	t297 = t254 * t324;
	t233 = -t282 * t255 - t260 * t297;
	t231 = 0.1e1 / t233 ^ 2;
	t213 = t219 * t231 + 0.1e1;
	t312 = t221 * t231;
	t303 = t259 * t260;
	t304 = t258 * t261;
	t280 = t257 * t304 + t303;
	t281 = t257 * t303 + t304;
	t290 = qJD(3) * t297;
	t217 = t258 * t290 + (t281 * qJD(2) + t280 * qJD(3)) * t255;
	t230 = 0.1e1 / t233;
	t313 = t217 * t230 * t231;
	t325 = -0.2e1 * (t197 * t312 - t219 * t313) / t213 ^ 2;
	t214 = atan2(-t221, t233);
	t209 = sin(t214);
	t210 = cos(t214);
	t191 = -t209 * t221 + t210 * t233;
	t188 = 0.1e1 / t191;
	t288 = t324 * t322;
	t275 = t259 * t288 - t323 * t261;
	t274 = t323 * t259 + t261 * t288;
	t296 = t255 * t322;
	t291 = t254 * t296;
	t278 = -t257 * t274 + t291;
	t225 = t278 * t258 - t260 * t275;
	t235 = t254 * t274 + t257 * t296;
	t253 = sin(pkin(13));
	t256 = cos(pkin(13));
	t208 = t225 * t256 + t235 * t253;
	t202 = 0.1e1 / t208;
	t189 = 0.1e1 / t191 ^ 2;
	t203 = 0.1e1 / t208 ^ 2;
	t211 = 0.1e1 / t213;
	t181 = (-t197 * t230 + t217 * t312) * t211;
	t287 = -t209 * t233 - t210 * t221;
	t177 = t287 * t181 - t209 * t197 + t210 * t217;
	t321 = t177 * t188 * t189;
	t310 = t275 * t258;
	t224 = -t260 * t291 + t274 * t306 - t310;
	t320 = t189 * t224;
	t243 = t274 * qJD(2);
	t244 = t275 * qJD(2);
	t307 = t257 * t258;
	t200 = t244 * t307 - t243 * t260 + (t278 * t260 + t310) * qJD(3);
	t309 = t253 * t254;
	t196 = t200 * t256 - t244 * t309;
	t319 = t196 * t202 * t203;
	t318 = t202 * t253;
	t207 = t225 * t253 - t235 * t256;
	t317 = t203 * t207;
	t316 = t207 * t256;
	t315 = t209 * t224;
	t314 = t210 * t224;
	t308 = t254 * t256;
	t220 = t224 ^ 2;
	t187 = t189 * t220 + 0.1e1;
	t199 = t225 * qJD(3) - t243 * t258 - t244 * t306;
	t301 = 0.2e1 * (t199 * t320 - t220 * t321) / t187 ^ 2;
	t201 = t207 ^ 2;
	t194 = t201 * t203 + 0.1e1;
	t195 = t200 * t253 + t244 * t308;
	t300 = 0.2e1 * (t195 * t317 - t201 * t319) / t194 ^ 2;
	t299 = t207 * t319;
	t298 = qJD(3) * t311;
	t294 = 0.2e1 * t224 * t321;
	t293 = -0.2e1 * t221 * t313;
	t223 = -t258 * t292 + t268;
	t234 = t280 * t255 + t258 * t297;
	t284 = -t223 * t230 + t234 * t312;
	t227 = -t276 * t306 + t272;
	t241 = t281 * t255;
	t283 = -t227 * t230 + t241 * t312;
	t228 = -t258 * t274 - t275 * t306;
	t229 = -t260 * t274 + t275 * t307;
	t279 = -t257 * t305 + t302;
	t226 = (t282 * qJD(2) + t279 * qJD(3)) * t255;
	t218 = t260 * t290 + (t279 * qJD(2) + t282 * qJD(3)) * t255;
	t216 = t229 * t256 - t275 * t309;
	t215 = t229 * t253 + t275 * t308;
	t206 = -t228 * qJD(3) + t243 * t307 + t244 * t260;
	t205 = t242 * t306 + t257 * t298 + t326;
	t198 = t326 * t257 + t327 * t260 + t298;
	t192 = 0.1e1 / t194;
	t185 = 0.1e1 / t187;
	t183 = t283 * t211;
	t182 = t284 * t211;
	t179 = t287 * t183 - t209 * t227 + t210 * t241;
	t178 = t287 * t182 - t209 * t223 + t210 * t234;
	t176 = t283 * t325 + (t241 * t293 - t205 * t230 + (t197 * t241 + t217 * t227 + t221 * t226) * t231) * t211;
	t175 = t284 * t325 + (t234 * t293 - t198 * t230 + (t197 * t234 + t217 * t223 + t218 * t221) * t231) * t211;
	t1 = [0, t176, t175, 0, 0, 0; 0, (t179 * t320 - t188 * t228) * t301 + ((t229 * qJD(3) - t243 * t306 + t244 * t258) * t188 + t179 * t294 + (-t228 * t177 - t179 * t199 - (-t176 * t221 - t183 * t197 + t226 + (-t183 * t233 - t227) * t181) * t314 - (-t176 * t233 - t183 * t217 - t205 + (t183 * t221 - t241) * t181) * t315) * t189) * t185, (t178 * t320 - t188 * t225) * t301 + (t178 * t294 + t200 * t188 + (-t225 * t177 - t178 * t199 - (-t175 * t221 - t182 * t197 + t218 + (-t182 * t233 - t223) * t181) * t314 - (-t175 * t233 - t182 * t217 - t198 + (t182 * t221 - t234) * t181) * t315) * t189) * t185, 0, 0, 0; 0, (-t202 * t215 + t216 * t317) * t300 + ((t206 * t253 + t243 * t308) * t202 + 0.2e1 * t216 * t299 + (-t215 * t196 - (t206 * t256 - t243 * t309) * t207 - t216 * t195) * t203) * t192, (-t203 * t316 + t318) * t224 * t300 + (-0.2e1 * t224 * t256 * t299 - t199 * t318 + (t199 * t316 + (t195 * t256 + t196 * t253) * t224) * t203) * t192, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:22
	% DurationCPUTime: 1.75s
	% Computational Cost: add. (5650->131), mult. (16984->252), div. (538->12), fcn. (21480->15), ass. (0->123)
	t297 = sin(qJ(2));
	t299 = cos(qJ(2));
	t363 = cos(pkin(12));
	t364 = cos(pkin(6));
	t329 = t364 * t363;
	t362 = sin(pkin(12));
	t315 = -t362 * t297 + t299 * t329;
	t279 = t315 * qJD(2);
	t293 = sin(pkin(7));
	t294 = sin(pkin(6));
	t332 = t293 * t294 * t363;
	t367 = -qJD(3) * t332 + t279;
	t296 = sin(qJ(3));
	t314 = -t297 * t329 - t362 * t299;
	t309 = qJD(2) * t314;
	t298 = cos(qJ(3));
	t311 = t298 * t315;
	t366 = qJD(3) * t311 + t296 * t309;
	t295 = cos(pkin(7));
	t310 = t315 * t296;
	t306 = t295 * t310 - t298 * t314;
	t347 = t295 * t298;
	t233 = t306 * qJD(3) + t367 * t296 - t309 * t347;
	t352 = t314 * t296;
	t257 = -t295 * t311 + t298 * t332 - t352;
	t255 = t257 ^ 2;
	t343 = t298 * t299;
	t346 = t296 * t297;
	t320 = t295 * t343 - t346;
	t337 = t293 * t364;
	t270 = -t320 * t294 - t298 * t337;
	t268 = 0.1e1 / t270 ^ 2;
	t249 = t255 * t268 + 0.1e1;
	t353 = t257 * t268;
	t344 = t297 * t298;
	t345 = t296 * t299;
	t318 = t295 * t345 + t344;
	t319 = t295 * t344 + t345;
	t330 = qJD(3) * t337;
	t253 = t296 * t330 + (t319 * qJD(2) + t318 * qJD(3)) * t294;
	t267 = 0.1e1 / t270;
	t354 = t253 * t267 * t268;
	t365 = -0.2e1 * (t233 * t353 - t255 * t354) / t249 ^ 2;
	t250 = atan2(-t257, t270);
	t245 = sin(t250);
	t246 = cos(t250);
	t227 = -t245 * t257 + t246 * t270;
	t224 = 0.1e1 / t227;
	t328 = t364 * t362;
	t313 = t297 * t328 - t363 * t299;
	t312 = t363 * t297 + t299 * t328;
	t336 = t294 * t362;
	t331 = t293 * t336;
	t316 = -t295 * t312 + t331;
	t261 = t316 * t296 - t298 * t313;
	t272 = t293 * t312 + t295 * t336;
	t292 = pkin(13) + qJ(5);
	t290 = sin(t292);
	t291 = cos(t292);
	t242 = t261 * t291 + t272 * t290;
	t238 = 0.1e1 / t242;
	t225 = 0.1e1 / t227 ^ 2;
	t239 = 0.1e1 / t242 ^ 2;
	t247 = 0.1e1 / t249;
	t217 = (-t233 * t267 + t253 * t353) * t247;
	t327 = -t245 * t270 - t246 * t257;
	t213 = t327 * t217 - t245 * t233 + t246 * t253;
	t361 = t213 * t224 * t225;
	t280 = t312 * qJD(2);
	t281 = t313 * qJD(2);
	t348 = t295 * t296;
	t351 = t313 * t296;
	t236 = t281 * t348 - t280 * t298 + (t316 * t298 + t351) * qJD(3);
	t349 = t291 * t293;
	t228 = t242 * qJD(5) + t236 * t290 + t281 * t349;
	t241 = t261 * t290 - t272 * t291;
	t237 = t241 ^ 2;
	t232 = t237 * t239 + 0.1e1;
	t357 = t239 * t241;
	t342 = qJD(5) * t241;
	t350 = t290 * t293;
	t229 = t236 * t291 - t281 * t350 - t342;
	t358 = t229 * t238 * t239;
	t360 = (t228 * t357 - t237 * t358) / t232 ^ 2;
	t260 = -t298 * t331 + t312 * t347 - t351;
	t359 = t225 * t260;
	t356 = t245 * t260;
	t355 = t246 * t260;
	t256 = t260 ^ 2;
	t223 = t256 * t225 + 0.1e1;
	t235 = t261 * qJD(3) - t280 * t296 - t281 * t347;
	t341 = 0.2e1 * (t235 * t359 - t256 * t361) / t223 ^ 2;
	t340 = -0.2e1 * t360;
	t339 = t241 * t358;
	t338 = qJD(3) * t352;
	t334 = 0.2e1 * t260 * t361;
	t333 = -0.2e1 * t257 * t354;
	t324 = -t290 * t238 + t291 * t357;
	t259 = -t296 * t332 + t306;
	t271 = t318 * t294 + t296 * t337;
	t323 = -t259 * t267 + t271 * t353;
	t263 = -t314 * t347 + t310;
	t278 = t319 * t294;
	t322 = -t263 * t267 + t278 * t353;
	t265 = -t298 * t312 + t313 * t348;
	t321 = -t265 * t290 - t313 * t349;
	t252 = t265 * t291 - t313 * t350;
	t264 = -t296 * t312 - t313 * t347;
	t317 = -t295 * t346 + t343;
	t262 = (t320 * qJD(2) + t317 * qJD(3)) * t294;
	t254 = t298 * t330 + (t317 * qJD(2) + t320 * qJD(3)) * t294;
	t244 = -t264 * qJD(3) + t280 * t348 + t281 * t298;
	t243 = t279 * t347 + t295 * t338 + t366;
	t234 = t366 * t295 + t367 * t298 + t338;
	t230 = 0.1e1 / t232;
	t221 = 0.1e1 / t223;
	t219 = t322 * t247;
	t218 = t323 * t247;
	t215 = t327 * t219 - t245 * t263 + t246 * t278;
	t214 = t327 * t218 - t245 * t259 + t246 * t271;
	t212 = t322 * t365 + (t278 * t333 - t243 * t267 + (t233 * t278 + t253 * t263 + t257 * t262) * t268) * t247;
	t211 = t323 * t365 + (t271 * t333 - t234 * t267 + (t233 * t271 + t253 * t259 + t254 * t257) * t268) * t247;
	t1 = [0, t212, t211, 0, 0, 0; 0, (t215 * t359 - t224 * t264) * t341 + ((t265 * qJD(3) - t280 * t347 + t281 * t296) * t224 + t215 * t334 + (-t264 * t213 - t215 * t235 - (-t212 * t257 - t219 * t233 + t262 + (-t219 * t270 - t263) * t217) * t355 - (-t212 * t270 - t219 * t253 - t243 + (t219 * t257 - t278) * t217) * t356) * t225) * t221, (t214 * t359 - t224 * t261) * t341 + (t214 * t334 + t236 * t224 + (-t261 * t213 - t214 * t235 - (-t211 * t257 - t218 * t233 + t254 + (-t218 * t270 - t259) * t217) * t355 - (-t211 * t270 - t218 * t253 - t234 + (t218 * t257 - t271) * t217) * t356) * t225) * t221, 0, 0, 0; 0, 0.2e1 * (t238 * t321 + t252 * t357) * t360 + ((t252 * qJD(5) + t244 * t290 + t280 * t349) * t238 + 0.2e1 * t252 * t339 + (t321 * t229 - (t321 * qJD(5) + t244 * t291 - t280 * t350) * t241 - t252 * t228) * t239) * t230, t324 * t260 * t340 + (t324 * t235 + ((-qJD(5) * t238 - 0.2e1 * t339) * t291 + (t228 * t291 + (t229 - t342) * t290) * t239) * t260) * t230, 0, t340 + 0.2e1 * (t228 * t239 * t230 + (-t230 * t358 - t239 * t360) * t241) * t241, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:21
	% EndTime: 2019-10-09 22:35:25
	% DurationCPUTime: 4.19s
	% Computational Cost: add. (19059->203), mult. (44088->383), div. (816->12), fcn. (56595->17), ass. (0->162)
	t404 = sin(qJ(3));
	t407 = cos(qJ(3));
	t400 = sin(pkin(7));
	t402 = cos(pkin(7));
	t405 = sin(qJ(2));
	t408 = cos(qJ(2));
	t484 = cos(pkin(12));
	t485 = cos(pkin(6));
	t437 = t485 * t484;
	t483 = sin(pkin(12));
	t423 = t483 * t405 - t408 * t437;
	t401 = sin(pkin(6));
	t448 = t401 * t484;
	t415 = -t400 * t448 - t423 * t402;
	t422 = -t405 * t437 - t483 * t408;
	t372 = t404 * t422 + t415 * t407;
	t387 = t423 * qJD(2);
	t388 = t422 * qJD(2);
	t465 = t402 * t404;
	t349 = t372 * qJD(3) - t387 * t407 + t388 * t465;
	t373 = t415 * t404 - t407 * t422;
	t399 = pkin(13) + qJ(5);
	t397 = sin(t399);
	t398 = cos(t399);
	t416 = t423 * t400 - t402 * t448;
	t356 = t373 * t398 + t416 * t397;
	t467 = t400 * t398;
	t323 = t356 * qJD(5) + t349 * t397 + t388 * t467;
	t354 = t373 * t397 - t416 * t398;
	t352 = t354 ^ 2;
	t460 = t405 * t407;
	t461 = t404 * t408;
	t426 = t402 * t461 + t460;
	t446 = t485 * t400;
	t383 = t426 * t401 + t404 * t446;
	t466 = t400 * t408;
	t391 = -t401 * t466 + t485 * t402;
	t367 = t383 * t397 - t391 * t398;
	t365 = 0.1e1 / t367 ^ 2;
	t339 = t352 * t365 + 0.1e1;
	t337 = 0.1e1 / t339;
	t368 = t383 * t398 + t391 * t397;
	t458 = t407 * t408;
	t462 = t404 * t405;
	t425 = -t402 * t462 + t458;
	t428 = t402 * t458 - t462;
	t439 = qJD(3) * t446;
	t370 = t407 * t439 + (t425 * qJD(2) + t428 * qJD(3)) * t401;
	t449 = t400 * t401 * t405;
	t441 = qJD(2) * t449;
	t341 = t368 * qJD(5) + t370 * t397 - t398 * t441;
	t364 = 0.1e1 / t367;
	t474 = t354 * t365;
	t306 = (-t323 * t364 + t341 * t474) * t337;
	t340 = atan2(-t354, t367);
	t333 = sin(t340);
	t334 = cos(t340);
	t435 = -t333 * t367 - t334 * t354;
	t301 = t435 * t306 - t333 * t323 + t334 * t341;
	t319 = -t333 * t354 + t334 * t367;
	t316 = 0.1e1 / t319;
	t317 = 0.1e1 / t319 ^ 2;
	t488 = t301 * t316 * t317;
	t436 = t485 * t483;
	t420 = t484 * t405 + t408 * t436;
	t447 = t401 * t483;
	t440 = t400 * t447;
	t417 = -t420 * t402 + t440;
	t421 = t405 * t436 - t484 * t408;
	t375 = t417 * t404 - t407 * t421;
	t418 = t420 * t400 + t402 * t447;
	t357 = t375 * t397 - t418 * t398;
	t444 = 0.2e1 * t357 * t488;
	t382 = t428 * t401 + t407 * t446;
	t430 = -t364 * t372 + t382 * t474;
	t487 = t397 * t430;
	t475 = t341 * t364 * t365;
	t486 = -0.2e1 * (t323 * t474 - t352 * t475) / t339 ^ 2;
	t358 = t375 * t398 + t418 * t397;
	t406 = cos(qJ(6));
	t419 = t420 * t407;
	t469 = t421 * t404;
	t374 = t402 * t419 - t407 * t440 - t469;
	t403 = sin(qJ(6));
	t472 = t374 * t403;
	t336 = t358 * t406 + t472;
	t330 = 0.1e1 / t336;
	t331 = 0.1e1 / t336 ^ 2;
	t389 = t420 * qJD(2);
	t390 = t421 * qJD(2);
	t351 = t390 * t465 - t389 * t407 + (t417 * t407 + t469) * qJD(3);
	t468 = t397 * t400;
	t326 = -t357 * qJD(5) + t351 * t398 - t390 * t468;
	t464 = t402 * t407;
	t350 = t375 * qJD(3) - t389 * t404 - t390 * t464;
	t314 = t336 * qJD(6) + t326 * t403 - t350 * t406;
	t471 = t374 * t406;
	t335 = t358 * t403 - t471;
	t329 = t335 ^ 2;
	t322 = t329 * t331 + 0.1e1;
	t478 = t331 * t335;
	t455 = qJD(6) * t335;
	t315 = t326 * t406 + t350 * t403 - t455;
	t481 = t315 * t330 * t331;
	t482 = (t314 * t478 - t329 * t481) / t322 ^ 2;
	t480 = t317 * t357;
	t325 = t358 * qJD(5) + t351 * t397 + t390 * t467;
	t479 = t325 * t317;
	t477 = t333 * t357;
	t476 = t334 * t357;
	t473 = t374 * t397;
	t463 = t403 * t330;
	t459 = t406 * t335;
	t457 = qJD(5) * t398;
	t456 = qJD(5) * t400;
	t353 = t357 ^ 2;
	t313 = t353 * t317 + 0.1e1;
	t454 = 0.2e1 * (-t353 * t488 + t357 * t479) / t313 ^ 2;
	t453 = -0.2e1 * t482;
	t452 = 0.2e1 * t482;
	t450 = t335 * t481;
	t443 = 0.2e1 * t450;
	t442 = -0.2e1 * t354 * t475;
	t438 = qJD(6) * t374 * t398 + t351;
	t378 = t421 * t465 - t419;
	t363 = t378 * t398 - t421 * t468;
	t377 = -t420 * t404 - t421 * t464;
	t344 = t363 * t406 + t377 * t403;
	t343 = t363 * t403 - t377 * t406;
	t433 = t331 * t459 - t463;
	t432 = -t356 * t364 + t368 * t474;
	t376 = -t423 * t407 + t422 * t465;
	t361 = t376 * t397 + t422 * t467;
	t386 = t425 * t401;
	t379 = t386 * t397 - t398 * t449;
	t431 = -t361 * t364 + t379 * t474;
	t429 = -t378 * t397 - t421 * t467;
	t427 = -t402 * t460 - t461;
	t424 = qJD(5) * t473 + qJD(6) * t375 - t350 * t398;
	t369 = -t404 * t439 + (t427 * qJD(2) - t426 * qJD(3)) * t401;
	t360 = -t377 * qJD(3) + t389 * t465 + t390 * t407;
	t359 = t378 * qJD(3) - t389 * t464 + t390 * t404;
	t348 = -t373 * qJD(3) + t387 * t404 + t388 * t464;
	t347 = t386 * t457 + ((t427 * qJD(3) + t405 * t456) * t397 + (-t426 * t397 - t398 * t466) * qJD(2)) * t401;
	t346 = t375 * t403 - t398 * t471;
	t345 = -t375 * t406 - t398 * t472;
	t342 = -t367 * qJD(5) + t370 * t398 + t397 * t441;
	t328 = t429 * qJD(5) + t360 * t398 - t389 * t468;
	t327 = t376 * t457 + t387 * t467 + (t387 * t465 + t388 * t407 + (t423 * t404 + t422 * t464) * qJD(3) - t422 * t456) * t397;
	t324 = -t354 * qJD(5) + t349 * t398 - t388 * t468;
	t320 = 0.1e1 / t322;
	t311 = 0.1e1 / t313;
	t310 = t337 * t487;
	t309 = t431 * t337;
	t308 = t432 * t337;
	t304 = (-t333 * t372 + t334 * t382) * t397 + t435 * t310;
	t303 = t435 * t309 - t333 * t361 + t334 * t379;
	t302 = t435 * t308 - t333 * t356 + t334 * t368;
	t300 = t431 * t486 + (t379 * t442 - t327 * t364 + (t323 * t379 + t341 * t361 + t347 * t354) * t365) * t337;
	t298 = t432 * t486 + (t368 * t442 - t324 * t364 + (t323 * t368 + t341 * t356 + t342 * t354) * t365) * t337;
	t297 = t486 * t487 + (t430 * t457 + (t382 * t442 - t348 * t364 + (t323 * t382 + t341 * t372 + t354 * t369) * t365) * t397) * t337;
	t1 = [0, t300, t297, 0, t298, 0; 0, (t303 * t480 + t316 * t429) * t454 + ((t363 * qJD(5) + t360 * t397 + t389 * t467) * t316 + t303 * t444 + (t429 * t301 - t303 * t325 - (-t300 * t354 - t309 * t323 + t347 + (-t309 * t367 - t361) * t306) * t476 - (-t300 * t367 - t309 * t341 - t327 + (t309 * t354 - t379) * t306) * t477) * t317) * t311, (t304 * t480 + t316 * t473) * t454 + ((-t350 * t397 - t374 * t457) * t316 + (-t479 + t444) * t304 + (t473 * t301 - (t382 * t457 - t297 * t354 - t310 * t323 + t369 * t397 + (-t310 * t367 - t372 * t397) * t306) * t476 - (-t372 * t457 - t297 * t367 - t310 * t341 - t348 * t397 + (t310 * t354 - t382 * t397) * t306) * t477) * t317) * t311, 0, (t302 * t480 - t316 * t358) * t454 + (t302 * t444 + t326 * t316 + (-t358 * t301 - t302 * t325 - (-t298 * t354 - t308 * t323 + t342 + (-t308 * t367 - t356) * t306) * t476 - (-t298 * t367 - t308 * t341 - t324 + (t308 * t354 - t368) * t306) * t477) * t317) * t311, 0; 0, (-t330 * t343 + t344 * t478) * t452 + ((t344 * qJD(6) + t328 * t403 - t359 * t406) * t330 + t344 * t443 + (-t343 * t315 - (-t343 * qJD(6) + t328 * t406 + t359 * t403) * t335 - t344 * t314) * t331) * t320, (-t330 * t345 + t346 * t478) * t452 + (t346 * t443 - t438 * t330 * t406 + t424 * t463 + (-t335 * t403 * t438 - t346 * t314 - t345 * t315 - t424 * t459) * t331) * t320, 0, t433 * t357 * t453 + (t433 * t325 + ((-qJD(6) * t330 - 0.2e1 * t450) * t406 + (t314 * t406 + (t315 - t455) * t403) * t331) * t357) * t320, t453 + 0.2e1 * (t314 * t331 * t320 + (-t320 * t481 - t331 * t482) * t335) * t335;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end