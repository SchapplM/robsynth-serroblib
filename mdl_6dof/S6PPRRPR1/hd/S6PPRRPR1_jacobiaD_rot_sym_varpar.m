% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR1
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
%   Wie in S6PPRRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiaD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (151->10), mult. (481->28), div. (18->4), fcn. (595->10), ass. (0->20)
	t84 = sin(pkin(12));
	t88 = cos(pkin(12));
	t89 = cos(pkin(11));
	t85 = sin(pkin(11));
	t98 = t85 * cos(pkin(6));
	t82 = -t84 * t98 + t89 * t88;
	t92 = sin(qJ(3));
	t93 = cos(qJ(3));
	t96 = t85 * sin(pkin(7)) * sin(pkin(6)) + (-t89 * t84 - t88 * t98) * cos(pkin(7));
	t78 = t82 * t92 - t96 * t93;
	t79 = t82 * t93 + t96 * t92;
	t76 = 0.1e1 / t79 ^ 2;
	t104 = qJD(3) * t76;
	t101 = t79 * t104;
	t102 = t78 / t79 * t104;
	t75 = t78 ^ 2;
	t72 = t75 * t76 + 0.1e1;
	t103 = (t78 * t101 + t75 * t102) / t72 ^ 2;
	t70 = 0.1e1 / t72;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, -0.2e1 * t103 + 0.2e1 * (t70 * t101 + (t70 * t102 - t76 * t103) * t78) * t78, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:41
	% DurationCPUTime: 0.97s
	% Computational Cost: add. (2591->65), mult. (8150->143), div. (281->12), fcn. (10563->15), ass. (0->78)
	t208 = sin(pkin(7));
	t209 = sin(pkin(6));
	t211 = cos(pkin(6));
	t207 = sin(pkin(12));
	t250 = sin(pkin(11));
	t234 = t250 * t207;
	t210 = cos(pkin(12));
	t251 = cos(pkin(11));
	t235 = t251 * t210;
	t252 = cos(pkin(7));
	t256 = (t211 * t235 - t234) * t252 - t208 * t209 * t251;
	t233 = t250 * t210;
	t236 = t251 * t207;
	t224 = t211 * t233 + t236;
	t238 = t209 * t250;
	t255 = -t208 * t238 + t224 * t252;
	t202 = t211 * t236 + t233;
	t213 = sin(qJ(3));
	t253 = cos(qJ(3));
	t188 = t202 * t253 + t213 * t256;
	t237 = t210 * t252;
	t241 = t208 * t211;
	t254 = (-t207 * t213 + t253 * t237) * t209 + t253 * t241;
	t186 = t202 * t213 - t253 * t256;
	t179 = atan2(-t186, -t254);
	t174 = sin(t179);
	t175 = cos(t179);
	t162 = -t174 * t186 - t175 * t254;
	t159 = 0.1e1 / t162;
	t203 = -t211 * t234 + t235;
	t190 = t203 * t253 - t255 * t213;
	t199 = t224 * t208 + t252 * t238;
	t212 = sin(qJ(4));
	t214 = cos(qJ(4));
	t173 = t190 * t214 + t199 * t212;
	t169 = 0.1e1 / t173;
	t194 = 0.1e1 / t254;
	t160 = 0.1e1 / t162 ^ 2;
	t170 = 0.1e1 / t173 ^ 2;
	t195 = 0.1e1 / t254 ^ 2;
	t184 = t186 ^ 2;
	t178 = t184 * t195 + 0.1e1;
	t176 = 0.1e1 / t178;
	t181 = t188 * qJD(3);
	t198 = t213 * t241 + (t253 * t207 + t213 * t237) * t209;
	t192 = t198 * qJD(3);
	t244 = t186 * t195;
	t153 = (t181 * t194 + t192 * t244) * t176;
	t229 = t174 * t254 - t175 * t186;
	t150 = t229 * t153 - t174 * t181 + t175 * t192;
	t249 = t150 * t159 * t160;
	t172 = t190 * t212 - t199 * t214;
	t168 = t172 ^ 2;
	t165 = t168 * t170 + 0.1e1;
	t189 = t203 * t213 + t255 * t253;
	t182 = t189 * qJD(3);
	t166 = t173 * qJD(4) - t182 * t212;
	t245 = t170 * t172;
	t240 = qJD(4) * t172;
	t167 = -t182 * t214 - t240;
	t246 = t167 * t169 * t170;
	t248 = (t166 * t245 - t168 * t246) / t165 ^ 2;
	t247 = t160 * t189;
	t243 = t186 * t198;
	t242 = t192 * t194 * t195;
	t239 = -0.2e1 * t248;
	t227 = -t169 * t212 + t214 * t245;
	t226 = t188 * t194 + t195 * t243;
	t191 = t254 * qJD(3);
	t185 = t189 ^ 2;
	t183 = t190 * qJD(3);
	t180 = t186 * qJD(3);
	t163 = 0.1e1 / t165;
	t157 = t185 * t160 + 0.1e1;
	t154 = t226 * t176;
	t151 = t229 * t154 - t174 * t188 + t175 * t198;
	t149 = -0.2e1 * t226 / t178 ^ 2 * (t181 * t244 + t184 * t242) + (0.2e1 * t242 * t243 - t180 * t194 + (t181 * t198 + t186 * t191 + t188 * t192) * t195) * t176;
	t1 = [0, 0, t149, 0, 0, 0; 0, 0, 0.2e1 * (t151 * t247 - t159 * t190) / t157 ^ 2 * (t183 * t247 - t185 * t249) + (-t182 * t159 + (-t190 * t150 - t151 * t183) * t160 + (0.2e1 * t151 * t249 + (-(-t149 * t186 - t154 * t181 + t191 + (t154 * t254 - t188) * t153) * t175 - (t149 * t254 - t154 * t192 + t180 + (t154 * t186 - t198) * t153) * t174) * t160) * t189) / t157, 0, 0, 0; 0, 0, t227 * t189 * t239 + (t227 * t183 + ((-qJD(4) * t169 - 0.2e1 * t172 * t246) * t214 + (t166 * t214 + (t167 - t240) * t212) * t170) * t189) * t163, t239 + 0.2e1 * (t163 * t166 * t170 + (-t163 * t246 - t170 * t248) * t172) * t172, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:43
	% DurationCPUTime: 2.31s
	% Computational Cost: add. (8432->112), mult. (25028->236), div. (535->12), fcn. (32914->17), ass. (0->114)
	t318 = sin(pkin(12));
	t319 = sin(pkin(11));
	t290 = t319 * t318;
	t321 = cos(pkin(12));
	t322 = cos(pkin(11));
	t293 = t322 * t321;
	t324 = cos(pkin(6));
	t261 = -t324 * t290 + t293;
	t267 = sin(qJ(3));
	t269 = cos(qJ(3));
	t291 = t319 * t321;
	t292 = t322 * t318;
	t282 = t324 * t291 + t292;
	t264 = sin(pkin(6));
	t297 = t264 * t319;
	t320 = sin(pkin(7));
	t323 = cos(pkin(7));
	t326 = t282 * t323 - t320 * t297;
	t245 = t261 * t267 + t326 * t269;
	t260 = t324 * t292 + t291;
	t281 = -t324 * t293 + t290;
	t298 = t264 * t320;
	t277 = -t281 * t323 - t322 * t298;
	t244 = t260 * t269 + t277 * t267;
	t266 = sin(qJ(4));
	t268 = cos(qJ(4));
	t276 = -t322 * t264 * t323 + t281 * t320;
	t236 = t244 * t268 + t276 * t266;
	t243 = -t260 * t267 + t277 * t269;
	t239 = t243 * qJD(3);
	t212 = t236 * qJD(4) + t239 * t266;
	t234 = t244 * t266 - t276 * t268;
	t232 = t234 ^ 2;
	t294 = t323 * t321;
	t295 = t324 * t320;
	t257 = (t267 * t294 + t318 * t269) * t264 + t267 * t295;
	t259 = -t321 * t298 + t324 * t323;
	t250 = t257 * t266 - t259 * t268;
	t248 = 0.1e1 / t250 ^ 2;
	t226 = t232 * t248 + 0.1e1;
	t224 = 0.1e1 / t226;
	t251 = t257 * t268 + t259 * t266;
	t256 = t269 * t295 + (-t318 * t267 + t269 * t294) * t264;
	t252 = t256 * qJD(3);
	t230 = t251 * qJD(4) + t252 * t266;
	t247 = 0.1e1 / t250;
	t308 = t234 * t248;
	t196 = (-t212 * t247 + t230 * t308) * t224;
	t227 = atan2(-t234, t250);
	t222 = sin(t227);
	t223 = cos(t227);
	t289 = -t222 * t250 - t223 * t234;
	t192 = t289 * t196 - t222 * t212 + t223 * t230;
	t206 = -t222 * t234 + t223 * t250;
	t203 = 0.1e1 / t206;
	t204 = 0.1e1 / t206 ^ 2;
	t329 = t192 * t203 * t204;
	t246 = t261 * t269 - t267 * t326;
	t278 = t282 * t320 + t323 * t297;
	t237 = t246 * t266 - t278 * t268;
	t328 = 0.2e1 * t237 * t329;
	t285 = -t243 * t247 + t256 * t308;
	t327 = t266 * t285;
	t309 = t230 * t247 * t248;
	t325 = -0.2e1 * (t212 * t308 - t232 * t309) / t226 ^ 2;
	t238 = t246 * t268 + t278 * t266;
	t263 = sin(pkin(13));
	t265 = cos(pkin(13));
	t221 = t238 * t265 + t245 * t263;
	t217 = 0.1e1 / t221;
	t218 = 0.1e1 / t221 ^ 2;
	t317 = t204 * t237;
	t241 = t245 * qJD(3);
	t215 = -t237 * qJD(4) - t241 * t268;
	t242 = t246 * qJD(3);
	t211 = t215 * t265 + t242 * t263;
	t316 = t211 * t217 * t218;
	t214 = t238 * qJD(4) - t241 * t266;
	t315 = t214 * t204;
	t314 = t217 * t263;
	t220 = t238 * t263 - t245 * t265;
	t313 = t218 * t220;
	t312 = t220 * t265;
	t311 = t222 * t237;
	t310 = t223 * t237;
	t307 = t245 * t266;
	t306 = t245 * t268;
	t303 = qJD(4) * t268;
	t233 = t237 ^ 2;
	t202 = t233 * t204 + 0.1e1;
	t302 = 0.2e1 * (-t233 * t329 + t237 * t315) / t202 ^ 2;
	t216 = t220 ^ 2;
	t209 = t216 * t218 + 0.1e1;
	t210 = t215 * t263 - t242 * t265;
	t301 = 0.2e1 * (t210 * t313 - t216 * t316) / t209 ^ 2;
	t299 = t220 * t316;
	t296 = -0.2e1 * t234 * t309;
	t286 = -t236 * t247 + t251 * t308;
	t284 = qJD(4) * t307 - t242 * t268;
	t253 = t257 * qJD(3);
	t240 = t244 * qJD(3);
	t231 = -t250 * qJD(4) + t252 * t268;
	t229 = t246 * t263 - t265 * t306;
	t228 = -t246 * t265 - t263 * t306;
	t213 = -t234 * qJD(4) + t239 * t268;
	t207 = 0.1e1 / t209;
	t200 = 0.1e1 / t202;
	t198 = t224 * t327;
	t197 = t286 * t224;
	t194 = (-t222 * t243 + t223 * t256) * t266 + t289 * t198;
	t193 = t289 * t197 - t222 * t236 + t223 * t251;
	t190 = t286 * t325 + (t251 * t296 - t213 * t247 + (t212 * t251 + t230 * t236 + t231 * t234) * t248) * t224;
	t189 = t325 * t327 + (t285 * t303 + (t256 * t296 + t240 * t247 + (t212 * t256 + t230 * t243 - t234 * t253) * t248) * t266) * t224;
	t1 = [0, 0, t189, t190, 0, 0; 0, 0, (t194 * t317 + t203 * t307) * t302 + ((-t242 * t266 - t245 * t303) * t203 + (-t315 + t328) * t194 + (t307 * t192 - (t256 * t303 - t189 * t234 - t198 * t212 - t253 * t266 + (-t198 * t250 - t243 * t266) * t196) * t310 - (-t243 * t303 - t189 * t250 - t198 * t230 + t240 * t266 + (t198 * t234 - t256 * t266) * t196) * t311) * t204) * t200, (t193 * t317 - t203 * t238) * t302 + (t193 * t328 + t215 * t203 + (-t238 * t192 - t193 * t214 - (-t190 * t234 - t197 * t212 + t231 + (-t197 * t250 - t236) * t196) * t310 - (-t190 * t250 - t197 * t230 - t213 + (t197 * t234 - t251) * t196) * t311) * t204) * t200, 0, 0; 0, 0, (-t217 * t228 + t229 * t313) * t301 + ((t241 * t265 + t284 * t263) * t217 + 0.2e1 * t229 * t299 + (-t228 * t211 - (-t241 * t263 + t284 * t265) * t220 - t229 * t210) * t218) * t207, (-t218 * t312 + t314) * t237 * t301 + (-0.2e1 * t237 * t265 * t299 - t214 * t314 + (t214 * t312 + (t210 * t265 + t211 * t263) * t237) * t218) * t207, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:41
	% EndTime: 2019-10-09 21:10:43
	% DurationCPUTime: 2.53s
	% Computational Cost: add. (9604->120), mult. (27507->244), div. (559->12), fcn. (36133->17), ass. (0->120)
	t380 = sin(pkin(12));
	t381 = sin(pkin(11));
	t349 = t381 * t380;
	t384 = cos(pkin(12));
	t385 = cos(pkin(11));
	t355 = t385 * t384;
	t387 = cos(pkin(6));
	t319 = -t387 * t349 + t355;
	t325 = sin(qJ(3));
	t327 = cos(qJ(3));
	t351 = t381 * t384;
	t354 = t385 * t380;
	t341 = t387 * t351 + t354;
	t383 = sin(pkin(6));
	t350 = t381 * t383;
	t382 = sin(pkin(7));
	t386 = cos(pkin(7));
	t389 = t341 * t386 - t382 * t350;
	t303 = t319 * t325 + t389 * t327;
	t318 = t387 * t354 + t351;
	t340 = -t387 * t355 + t349;
	t353 = t383 * t382;
	t335 = -t340 * t386 - t385 * t353;
	t302 = t318 * t327 + t335 * t325;
	t324 = sin(qJ(4));
	t326 = cos(qJ(4));
	t356 = t386 * t383;
	t334 = t340 * t382 - t385 * t356;
	t293 = t302 * t326 + t334 * t324;
	t301 = -t318 * t325 + t335 * t327;
	t297 = t301 * qJD(3);
	t275 = t293 * qJD(4) + t297 * t324;
	t291 = t302 * t324 - t334 * t326;
	t289 = t291 ^ 2;
	t339 = t384 * t356 + t387 * t382;
	t352 = t383 * t380;
	t315 = t339 * t325 + t327 * t352;
	t317 = -t384 * t353 + t387 * t386;
	t308 = t315 * t324 - t317 * t326;
	t306 = 0.1e1 / t308 ^ 2;
	t283 = t289 * t306 + 0.1e1;
	t281 = 0.1e1 / t283;
	t309 = t315 * t326 + t317 * t324;
	t314 = -t325 * t352 + t339 * t327;
	t310 = t314 * qJD(3);
	t287 = t309 * qJD(4) + t310 * t324;
	t305 = 0.1e1 / t308;
	t369 = t291 * t306;
	t253 = (-t275 * t305 + t287 * t369) * t281;
	t284 = atan2(-t291, t308);
	t279 = sin(t284);
	t280 = cos(t284);
	t348 = -t279 * t308 - t280 * t291;
	t249 = t348 * t253 - t279 * t275 + t280 * t287;
	t263 = -t279 * t291 + t280 * t308;
	t260 = 0.1e1 / t263;
	t261 = 0.1e1 / t263 ^ 2;
	t392 = t249 * t260 * t261;
	t304 = t319 * t327 - t325 * t389;
	t336 = t341 * t382 + t386 * t350;
	t294 = t304 * t324 - t336 * t326;
	t391 = 0.2e1 * t294 * t392;
	t344 = -t301 * t305 + t314 * t369;
	t390 = t324 * t344;
	t370 = t287 * t305 * t306;
	t388 = -0.2e1 * (t275 * t369 - t289 * t370) / t283 ^ 2;
	t295 = t304 * t326 + t336 * t324;
	t323 = pkin(13) + qJ(6);
	t321 = sin(t323);
	t322 = cos(t323);
	t274 = t295 * t322 + t303 * t321;
	t270 = 0.1e1 / t274;
	t271 = 0.1e1 / t274 ^ 2;
	t299 = t303 * qJD(3);
	t278 = -t294 * qJD(4) - t299 * t326;
	t300 = t304 * qJD(3);
	t264 = t274 * qJD(6) + t278 * t321 - t300 * t322;
	t273 = t295 * t321 - t303 * t322;
	t269 = t273 ^ 2;
	t268 = t269 * t271 + 0.1e1;
	t375 = t271 * t273;
	t363 = qJD(6) * t273;
	t265 = t278 * t322 + t300 * t321 - t363;
	t377 = t265 * t270 * t271;
	t379 = (t264 * t375 - t269 * t377) / t268 ^ 2;
	t378 = t261 * t294;
	t376 = t270 * t321;
	t374 = t273 * t322;
	t277 = t295 * qJD(4) - t299 * t324;
	t373 = t277 * t261;
	t372 = t279 * t294;
	t371 = t280 * t294;
	t368 = t303 * t324;
	t367 = t303 * t326;
	t364 = qJD(4) * t326;
	t290 = t294 ^ 2;
	t259 = t290 * t261 + 0.1e1;
	t362 = 0.2e1 * (-t290 * t392 + t294 * t373) / t259 ^ 2;
	t361 = -0.2e1 * t379;
	t359 = t273 * t377;
	t358 = -0.2e1 * t291 * t370;
	t357 = qJD(6) * t367 - t299;
	t346 = t271 * t374 - t376;
	t345 = -t293 * t305 + t309 * t369;
	t342 = qJD(4) * t368 + qJD(6) * t304 - t300 * t326;
	t311 = t315 * qJD(3);
	t298 = t302 * qJD(3);
	t288 = -t308 * qJD(4) + t310 * t326;
	t286 = t304 * t321 - t322 * t367;
	t285 = -t304 * t322 - t321 * t367;
	t276 = -t291 * qJD(4) + t297 * t326;
	t266 = 0.1e1 / t268;
	t257 = 0.1e1 / t259;
	t255 = t281 * t390;
	t254 = t345 * t281;
	t251 = (-t279 * t301 + t280 * t314) * t324 + t348 * t255;
	t250 = t348 * t254 - t279 * t293 + t280 * t309;
	t247 = t345 * t388 + (t309 * t358 - t276 * t305 + (t275 * t309 + t287 * t293 + t288 * t291) * t306) * t281;
	t246 = t388 * t390 + (t344 * t364 + (t314 * t358 + t298 * t305 + (t275 * t314 + t287 * t301 - t291 * t311) * t306) * t324) * t281;
	t1 = [0, 0, t246, t247, 0, 0; 0, 0, (t251 * t378 + t260 * t368) * t362 + ((-t300 * t324 - t303 * t364) * t260 + (-t373 + t391) * t251 + (t368 * t249 - (t314 * t364 - t246 * t291 - t255 * t275 - t311 * t324 + (-t255 * t308 - t301 * t324) * t253) * t371 - (-t301 * t364 - t246 * t308 - t255 * t287 + t298 * t324 + (t255 * t291 - t314 * t324) * t253) * t372) * t261) * t257, (t250 * t378 - t260 * t295) * t362 + (t250 * t391 + t278 * t260 + (-t295 * t249 - t250 * t277 - (-t247 * t291 - t254 * t275 + t288 + (-t254 * t308 - t293) * t253) * t371 - (-t247 * t308 - t254 * t287 - t276 + (t254 * t291 - t309) * t253) * t372) * t261) * t257, 0, 0; 0, 0, 0.2e1 * (-t270 * t285 + t286 * t375) * t379 + (0.2e1 * t286 * t359 - t357 * t270 * t322 + t342 * t376 + (-t357 * t273 * t321 - t286 * t264 - t285 * t265 - t342 * t374) * t271) * t266, t346 * t294 * t361 + (t346 * t277 + ((-qJD(6) * t270 - 0.2e1 * t359) * t322 + (t264 * t322 + (t265 - t363) * t321) * t271) * t294) * t266, 0, t361 + 0.2e1 * (t264 * t271 * t266 + (-t266 * t377 - t271 * t379) * t273) * t273;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end