% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR2
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
%   Wie in S6PPRRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
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
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:35
	% DurationCPUTime: 0.96s
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
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:36
	% DurationCPUTime: 2.04s
	% Computational Cost: add. (7669->100), mult. (22804->213), div. (524->12), fcn. (30024->15), ass. (0->101)
	t240 = cos(pkin(6));
	t237 = sin(pkin(12));
	t286 = sin(pkin(11));
	t265 = t286 * t237;
	t239 = cos(pkin(12));
	t288 = cos(pkin(11));
	t267 = t288 * t239;
	t235 = -t240 * t265 + t267;
	t242 = sin(qJ(3));
	t244 = cos(qJ(3));
	t264 = t286 * t239;
	t268 = t288 * t237;
	t256 = t240 * t264 + t268;
	t238 = sin(pkin(6));
	t270 = t238 * t286;
	t287 = sin(pkin(7));
	t289 = cos(pkin(7));
	t291 = t256 * t289 - t287 * t270;
	t218 = t235 * t242 + t291 * t244;
	t234 = t240 * t268 + t264;
	t255 = -t240 * t267 + t265;
	t271 = t238 * t287;
	t252 = -t255 * t289 - t288 * t271;
	t217 = t234 * t244 + t252 * t242;
	t241 = sin(qJ(4));
	t243 = cos(qJ(4));
	t251 = -t288 * t238 * t289 + t255 * t287;
	t206 = t217 * t243 + t251 * t241;
	t216 = -t234 * t242 + t252 * t244;
	t209 = t216 * qJD(3);
	t186 = t206 * qJD(4) + t209 * t241;
	t204 = t217 * t241 - t251 * t243;
	t201 = t204 ^ 2;
	t266 = t287 * t240;
	t269 = t239 * t289;
	t230 = (t237 * t244 + t242 * t269) * t238 + t242 * t266;
	t233 = -t239 * t271 + t240 * t289;
	t223 = t230 * t241 - t233 * t243;
	t221 = 0.1e1 / t223 ^ 2;
	t197 = t201 * t221 + 0.1e1;
	t195 = 0.1e1 / t197;
	t224 = t230 * t243 + t233 * t241;
	t229 = t244 * t266 + (-t237 * t242 + t244 * t269) * t238;
	t225 = t229 * qJD(3);
	t199 = t224 * qJD(4) + t225 * t241;
	t220 = 0.1e1 / t223;
	t280 = t204 * t221;
	t174 = (-t186 * t220 + t199 * t280) * t195;
	t198 = atan2(-t204, t223);
	t192 = sin(t198);
	t193 = cos(t198);
	t261 = -t192 * t223 - t193 * t204;
	t171 = t261 * t174 - t192 * t186 + t193 * t199;
	t185 = -t192 * t204 + t193 * t223;
	t182 = 0.1e1 / t185;
	t183 = 0.1e1 / t185 ^ 2;
	t294 = t171 * t182 * t183;
	t219 = t235 * t244 - t242 * t291;
	t231 = t256 * t287 + t289 * t270;
	t207 = t219 * t241 - t231 * t243;
	t293 = 0.2e1 * t207 * t294;
	t258 = -t216 * t220 + t229 * t280;
	t292 = t241 * t258;
	t281 = t199 * t220 * t221;
	t290 = -0.2e1 * (t186 * t280 - t201 * t281) / t197 ^ 2;
	t213 = 0.1e1 / t218;
	t214 = 0.1e1 / t218 ^ 2;
	t285 = t183 * t207;
	t208 = t219 * t243 + t231 * t241;
	t211 = t218 * qJD(3);
	t188 = t208 * qJD(4) - t211 * t241;
	t284 = t188 * t183;
	t283 = t192 * t207;
	t282 = t193 * t207;
	t279 = t208 * t219;
	t278 = t218 * t241;
	t275 = qJD(4) * t243;
	t202 = t207 ^ 2;
	t180 = t202 * t183 + 0.1e1;
	t274 = 0.2e1 * (-t202 * t294 + t207 * t284) / t180 ^ 2;
	t189 = -t207 * qJD(4) - t211 * t243;
	t203 = t208 ^ 2;
	t194 = t203 * t214 + 0.1e1;
	t212 = t219 * qJD(3);
	t215 = t213 * t214;
	t273 = 0.2e1 * (t208 * t214 * t189 - t203 * t215 * t212) / t194 ^ 2;
	t263 = -0.2e1 * t204 * t281;
	t259 = -t206 * t220 + t224 * t280;
	t226 = t230 * qJD(3);
	t210 = t217 * qJD(3);
	t200 = -t223 * qJD(4) + t225 * t243;
	t190 = 0.1e1 / t194;
	t187 = -t204 * qJD(4) + t209 * t243;
	t178 = 0.1e1 / t180;
	t176 = t195 * t292;
	t175 = t259 * t195;
	t173 = (-t192 * t216 + t193 * t229) * t241 + t261 * t176;
	t172 = t261 * t175 - t192 * t206 + t193 * t224;
	t169 = t259 * t290 + (t224 * t263 - t187 * t220 + (t186 * t224 + t199 * t206 + t200 * t204) * t221) * t195;
	t168 = t290 * t292 + (t258 * t275 + (t229 * t263 + t210 * t220 + (t186 * t229 + t199 * t216 - t204 * t226) * t221) * t241) * t195;
	t1 = [0, 0, t168, t169, 0, 0; 0, 0, (t173 * t285 + t182 * t278) * t274 + ((-t212 * t241 - t218 * t275) * t182 + (-t284 + t293) * t173 + (t278 * t171 - (t229 * t275 - t168 * t204 - t176 * t186 - t226 * t241 + (-t176 * t223 - t216 * t241) * t174) * t282 - (-t216 * t275 - t168 * t223 - t176 * t199 + t210 * t241 + (t176 * t204 - t229 * t241) * t174) * t283) * t183) * t178, (t172 * t285 - t182 * t208) * t274 + (t172 * t293 + t189 * t182 + (-t208 * t171 - t172 * t188 - (-t169 * t204 - t175 * t186 + t200 + (-t175 * t223 - t206) * t174) * t282 - (-t169 * t223 - t175 * t199 - t187 + (t175 * t204 - t224) * t174) * t283) * t183) * t178, 0, 0; 0, 0, (t213 * t218 * t243 + t214 * t279) * t273 + (qJD(4) * t213 * t278 + (-t189 * t219 + t208 * t211) * t214 + (0.2e1 * t215 * t279 + (t214 * t218 - t213) * t243) * t212) * t190, t207 * t213 * t273 + (t207 * t212 * t214 - t188 * t213) * t190, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:37
	% DurationCPUTime: 2.56s
	% Computational Cost: add. (9293->121), mult. (27507->247), div. (559->12), fcn. (36133->17), ass. (0->119)
	t319 = cos(pkin(12));
	t320 = cos(pkin(11));
	t322 = cos(pkin(6));
	t316 = sin(pkin(12));
	t380 = sin(pkin(11));
	t353 = t380 * t316;
	t314 = t320 * t319 - t322 * t353;
	t325 = sin(qJ(3));
	t328 = cos(qJ(3));
	t317 = sin(pkin(7));
	t321 = cos(pkin(7));
	t352 = t380 * t319;
	t340 = t320 * t316 + t322 * t352;
	t318 = sin(pkin(6));
	t354 = t318 * t380;
	t383 = -t317 * t354 + t340 * t321;
	t295 = t314 * t325 + t383 * t328;
	t360 = t320 * t322;
	t313 = t316 * t360 + t352;
	t312 = t319 * t360 - t353;
	t362 = t318 * t320;
	t342 = -t312 * t321 + t317 * t362;
	t294 = -t313 * t325 - t342 * t328;
	t290 = t294 * qJD(3);
	t308 = -t312 * t317 - t321 * t362;
	t327 = cos(qJ(4));
	t324 = sin(qJ(4));
	t381 = -t313 * t328 + t342 * t325;
	t336 = t381 * t324;
	t262 = qJD(4) * t336 + (t308 * qJD(4) + t290) * t327;
	t285 = t308 * t324 - t327 * t381;
	t282 = t285 ^ 2;
	t361 = t319 * t321;
	t363 = t317 * t322;
	t307 = (t316 * t328 + t325 * t361) * t318 + t325 * t363;
	t311 = -t318 * t319 * t317 + t322 * t321;
	t301 = t307 * t327 + t311 * t324;
	t298 = 0.1e1 / t301 ^ 2;
	t275 = t282 * t298 + 0.1e1;
	t273 = 0.1e1 / t275;
	t300 = -t307 * t324 + t311 * t327;
	t306 = t328 * t363 + (-t316 * t325 + t328 * t361) * t318;
	t302 = t306 * qJD(3);
	t281 = t300 * qJD(4) + t302 * t327;
	t297 = 0.1e1 / t301;
	t369 = t285 * t298;
	t245 = (-t262 * t297 + t281 * t369) * t273;
	t276 = atan2(-t285, t301);
	t271 = sin(t276);
	t272 = cos(t276);
	t347 = -t271 * t301 - t272 * t285;
	t241 = t347 * t245 - t271 * t262 + t272 * t281;
	t255 = -t271 * t285 + t272 * t301;
	t252 = 0.1e1 / t255;
	t253 = 0.1e1 / t255 ^ 2;
	t387 = t241 * t252 * t253;
	t296 = t314 * t328 - t325 * t383;
	t309 = t340 * t317 + t321 * t354;
	t287 = t296 * t324 - t309 * t327;
	t323 = sin(qJ(6));
	t326 = cos(qJ(6));
	t346 = t287 * t326 - t295 * t323;
	t386 = t346 * qJD(6);
	t288 = t296 * t327 + t309 * t324;
	t385 = 0.2e1 * t288 * t387;
	t343 = -t294 * t297 + t306 * t369;
	t384 = t327 * t343;
	t370 = t281 * t297 * t298;
	t382 = -0.2e1 * (t262 * t369 - t282 * t370) / t275 ^ 2;
	t367 = t295 * t326;
	t270 = t287 * t323 + t367;
	t266 = 0.1e1 / t270;
	t267 = 0.1e1 / t270 ^ 2;
	t292 = t295 * qJD(3);
	t263 = t288 * qJD(4) - t292 * t324;
	t293 = t296 * qJD(3);
	t256 = t270 * qJD(6) - t263 * t326 + t293 * t323;
	t265 = t346 ^ 2;
	t260 = t265 * t267 + 0.1e1;
	t374 = t267 * t346;
	t257 = t263 * t323 + t293 * t326 + t386;
	t377 = t257 * t266 * t267;
	t379 = (-t256 * t374 - t265 * t377) / t260 ^ 2;
	t378 = t253 * t288;
	t264 = -t287 * qJD(4) - t292 * t327;
	t376 = t264 * t253;
	t375 = t266 * t326;
	t373 = t346 * t323;
	t372 = t271 * t288;
	t371 = t272 * t288;
	t368 = t295 * t324;
	t366 = t295 * t327;
	t359 = qJD(4) * t324;
	t283 = t288 ^ 2;
	t251 = t283 * t253 + 0.1e1;
	t358 = 0.2e1 * (-t283 * t387 + t288 * t376) / t251 ^ 2;
	t357 = 0.2e1 * t379;
	t351 = -0.2e1 * t346 * t377;
	t350 = -0.2e1 * t285 * t370;
	t348 = -qJD(6) * t368 - t292;
	t345 = -t267 * t373 + t375;
	t284 = -t308 * t327 - t336;
	t344 = t284 * t297 + t300 * t369;
	t339 = qJD(4) * t366 + qJD(6) * t296 + t293 * t324;
	t303 = t307 * qJD(3);
	t291 = t381 * qJD(3);
	t280 = -t301 * qJD(4) - t302 * t324;
	t278 = t296 * t326 - t323 * t368;
	t277 = t296 * t323 + t324 * t367;
	t261 = t285 * qJD(4) + t290 * t324;
	t258 = 0.1e1 / t260;
	t249 = 0.1e1 / t251;
	t247 = t273 * t384;
	t246 = t344 * t273;
	t243 = (-t271 * t294 + t272 * t306) * t327 + t347 * t247;
	t242 = t347 * t246 + t271 * t284 + t272 * t300;
	t239 = t344 * t382 + (t300 * t350 + t261 * t297 + (t262 * t300 + t280 * t285 - t281 * t284) * t298) * t273;
	t238 = t382 * t384 + (-t343 * t359 + (t306 * t350 - t291 * t297 + (t262 * t306 + t281 * t294 - t285 * t303) * t298) * t327) * t273;
	t1 = [0, 0, t238, t239, 0, 0; 0, 0, (t243 * t378 + t252 * t366) * t358 + ((-t293 * t327 + t295 * t359) * t252 + (-t376 + t385) * t243 + (t366 * t241 - (-t306 * t359 - t238 * t285 - t247 * t262 - t303 * t327 + (-t247 * t301 - t294 * t327) * t245) * t371 - (t294 * t359 - t238 * t301 - t247 * t281 - t291 * t327 + (t247 * t285 - t306 * t327) * t245) * t372) * t253) * t249, (t242 * t378 + t252 * t287) * t358 + (t242 * t385 - t263 * t252 + (t287 * t241 - t242 * t264 - (-t239 * t285 - t246 * t262 + t280 + (-t246 * t301 + t284) * t245) * t371 - (-t239 * t301 - t246 * t281 + t261 + (t246 * t285 - t300) * t245) * t372) * t253) * t249, 0, 0; 0, 0, (-t266 * t277 - t278 * t374) * t357 + (t278 * t351 + t348 * t266 * t323 + t339 * t375 + (t326 * t346 * t348 - t278 * t256 - t277 * t257 - t339 * t373) * t267) * t258, t345 * t288 * t357 + (-t345 * t264 + ((qJD(6) * t266 + t351) * t323 + (-t256 * t323 + (t257 + t386) * t326) * t267) * t288) * t258, 0, -0.2e1 * t379 - 0.2e1 * (t256 * t267 * t258 - (-t258 * t377 - t267 * t379) * t346) * t346;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end