% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR9
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
%   Wie in S6RPRRPR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiaD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (137->30), mult. (614->95), div. (108->12), fcn. (792->9), ass. (0->49)
	t86 = sin(pkin(6));
	t79 = t86 ^ 2;
	t88 = cos(pkin(6));
	t81 = 0.1e1 / t88 ^ 2;
	t90 = cos(qJ(1));
	t84 = t90 ^ 2;
	t77 = t79 * t81 * t84 + 0.1e1;
	t89 = sin(qJ(1));
	t83 = t89 ^ 2;
	t108 = 0.1e1 / t77 ^ 2 * t83;
	t112 = t108 * t81;
	t103 = t90 * t86;
	t76 = atan2(t103, t88);
	t72 = sin(t76);
	t73 = cos(t76);
	t58 = t72 * t103 + t73 * t88;
	t55 = 0.1e1 / t58;
	t85 = sin(pkin(12));
	t105 = t89 * t85;
	t87 = cos(pkin(12));
	t99 = t88 * t105 - t87 * t90;
	t65 = 0.1e1 / t99;
	t80 = 0.1e1 / t88;
	t56 = 0.1e1 / t58 ^ 2;
	t66 = 0.1e1 / t99 ^ 2;
	t111 = t56 * t89;
	t104 = t89 * t87;
	t70 = t88 * t104 + t85 * t90;
	t110 = t66 * t70;
	t106 = t88 * t90;
	t69 = -t85 * t106 - t104;
	t109 = t69 * t70;
	t107 = t79 * t80;
	t102 = qJD(1) * t90;
	t74 = 0.1e1 / t77;
	t101 = (t74 - 0.1e1) * t86;
	t100 = -0.2e1 * t80 * t112;
	t68 = t87 * t106 - t105;
	t51 = (-t73 * t74 * t90 * t107 + t72 * t101) * t89;
	t78 = t86 * t79;
	t67 = t65 * t66;
	t64 = t70 ^ 2;
	t63 = t69 * qJD(1);
	t62 = t68 * qJD(1);
	t61 = t64 * t66 + 0.1e1;
	t57 = t55 * t56;
	t54 = t56 * t79 * t83 + 0.1e1;
	t50 = qJD(1) * t51;
	t1 = [(-t74 * t80 * t86 + t78 * t100) * t102, 0, 0, 0, 0, 0; (0.2e1 * (t51 * t111 - t55 * t90) / t54 ^ 2 * (-t50 * t57 * t83 + t102 * t111) * t79 + ((0.2e1 * t51 * t57 * t89 - t56 * t90) * t50 + (-t89 * t55 + ((-t51 + (-t78 * t112 - t101) * t89 * t72) * t90 - (t84 * t79 ^ 2 * t100 + (-t108 + (0.2e1 * t83 - t84) * t74) * t107) * t89 * t73) * t56) * qJD(1)) / t54) * t86, 0, 0, 0, 0, 0; 0.2e1 * (t66 * t109 + t65 * t68) / t61 ^ 2 * (t63 * t64 * t67 + t62 * t110) + (-t69 * t62 * t66 + (-0.2e1 * t67 * t109 - t68 * t66) * t63 + (-t99 * t110 + t70 * t65) * qJD(1)) / t61, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (994->58), mult. (3107->139), div. (132->12), fcn. (4021->13), ass. (0->80)
	t183 = sin(pkin(7));
	t184 = sin(pkin(6));
	t185 = cos(pkin(12));
	t186 = cos(pkin(7));
	t187 = cos(pkin(6));
	t173 = -t184 * t185 * t183 + t187 * t186;
	t191 = cos(qJ(1));
	t211 = t191 * t185;
	t182 = sin(pkin(12));
	t189 = sin(qJ(1));
	t214 = t189 * t182;
	t174 = -t187 * t211 + t214;
	t216 = t184 * t191;
	t204 = -t174 * t183 + t186 * t216;
	t157 = atan2(t204, t173);
	t152 = sin(t157);
	t153 = cos(t157);
	t138 = t152 * t204 + t153 * t173;
	t135 = 0.1e1 / t138;
	t188 = sin(qJ(3));
	t190 = cos(qJ(3));
	t200 = t187 * t214 - t211;
	t212 = t191 * t182;
	t213 = t189 * t185;
	t201 = t187 * t213 + t212;
	t217 = t184 * t189;
	t209 = t183 * t217;
	t202 = -t186 * t201 + t209;
	t151 = t202 * t188 - t190 * t200;
	t145 = 0.1e1 / t151;
	t230 = t204 ^ 2;
	t170 = 0.1e1 / t173;
	t136 = 0.1e1 / t138 ^ 2;
	t146 = 0.1e1 / t151 ^ 2;
	t171 = 0.1e1 / t173 ^ 2;
	t229 = -0.2e1 * t170 * t171;
	t175 = -t187 * t212 - t213;
	t167 = t175 * qJD(1);
	t166 = t174 * qJD(1);
	t210 = qJD(1) * t184;
	t206 = t191 * t210;
	t199 = t166 * t186 + t183 * t206;
	t139 = t151 * qJD(3) + t167 * t188 - t199 * t190;
	t215 = t186 * t190;
	t218 = t200 * t188;
	t150 = -t190 * t209 + t201 * t215 - t218;
	t144 = t150 ^ 2;
	t143 = t144 * t146 + 0.1e1;
	t225 = t146 * t150;
	t140 = t167 * t190 + t199 * t188 + (t202 * t190 + t218) * qJD(3);
	t226 = t140 * t145 * t146;
	t228 = (t139 * t225 - t144 * t226) / t143 ^ 2;
	t168 = t201 * qJD(1);
	t207 = t189 * t210;
	t159 = -t168 * t183 - t186 * t207;
	t156 = t230 * t171 + 0.1e1;
	t154 = 0.1e1 / t156;
	t198 = t152 + (t153 * t170 * t204 - t152) * t154;
	t130 = t198 * t159;
	t227 = t130 * t135 * t136;
	t208 = t183 * t216;
	t203 = t174 * t186 + t208;
	t149 = t175 * t190 + t203 * t188;
	t224 = t149 * t150;
	t223 = t154 * t170;
	t155 = 0.1e1 / t156 ^ 2;
	t222 = t155 * t204;
	t158 = t166 * t183 - t186 * t206;
	t221 = t158 * t136;
	t163 = -t183 * t201 - t186 * t217;
	t220 = t159 * t163;
	t219 = t175 * t188;
	t205 = t183 * t207;
	t169 = t200 * qJD(1);
	t160 = t163 ^ 2;
	t148 = -t203 * t190 + t219;
	t141 = 0.1e1 / t143;
	t134 = t160 * t136 + 0.1e1;
	t131 = t198 * t163;
	t1 = [t220 * t222 * t229 + t158 * t223, 0, 0, 0, 0, 0; 0.2e1 * (-t131 * t136 * t163 - t135 * t204) / t134 ^ 2 * (-t160 * t227 + t163 * t221) + (t159 * t135 + (-t130 * t204 + t131 * t158) * t136 + (-0.2e1 * t131 * t227 + t198 * t221 + (t152 * t171 * t222 + (0.2e1 * t223 + (t230 * t229 - t170) * t155) * t153) * t136 * t220) * t163) / t134, 0, 0, 0, 0, 0; 0.2e1 * (-t145 * t148 + t146 * t224) * t228 + ((-t168 * t215 + t169 * t188 + t190 * t205) * t145 + 0.2e1 * t224 * t226 + (-t148 * t140 - (t169 * t190 + (t168 * t186 - t205) * t188) * t150 - t149 * t139) * t146 + (t149 * t145 - (t174 * t215 + t190 * t208 - t219) * t225) * qJD(3)) * t141, 0, -0.2e1 * t228 + 0.2e1 * (t139 * t146 * t141 + (-t141 * t226 - t146 * t228) * t150) * t150, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:34
	% DurationCPUTime: 1.94s
	% Computational Cost: add. (4837->111), mult. (15324->224), div. (448->12), fcn. (19446->15), ass. (0->108)
	t347 = sin(pkin(12));
	t352 = cos(pkin(6));
	t320 = t352 * t347;
	t350 = cos(pkin(12));
	t353 = sin(qJ(1));
	t355 = cos(qJ(1));
	t276 = t320 * t355 + t350 * t353;
	t287 = sin(qJ(3));
	t348 = sin(pkin(7));
	t349 = sin(pkin(6));
	t318 = t349 * t348;
	t311 = t355 * t318;
	t283 = t287 * t311;
	t322 = t352 * t350;
	t275 = -t322 * t355 + t347 * t353;
	t351 = cos(pkin(7));
	t326 = t275 * t351;
	t354 = cos(qJ(3));
	t362 = -t276 * t354 + t287 * t326 + t283;
	t304 = t354 * t311;
	t324 = t351 * t354;
	t361 = -t275 * t324 - t304;
	t249 = t276 * t287 - t361;
	t317 = t349 * t347;
	t319 = t351 * t349;
	t359 = t350 * t319 + t352 * t348;
	t295 = -t287 * t317 + t359 * t354;
	t246 = atan2(-t249, -t295);
	t241 = sin(t246);
	t242 = cos(t246);
	t222 = -t241 * t249 - t242 * t295;
	t220 = 0.1e1 / t222 ^ 2;
	t277 = -t353 * t320 + t355 * t350;
	t300 = t322 * t353 + t347 * t355;
	t298 = t300 * t351;
	t308 = t353 * t318;
	t303 = t354 * t308;
	t358 = -t277 * t287 - t354 * t298 + t303;
	t248 = t358 ^ 2;
	t218 = t248 * t220 + 0.1e1;
	t255 = t277 * t354 + (-t298 + t308) * t287;
	t272 = t276 * qJD(1);
	t297 = qJD(1) * t326;
	t228 = -qJD(1) * t304 + qJD(3) * t255 - t272 * t287 - t297 * t354;
	t339 = t228 * t220;
	t219 = 0.1e1 / t222;
	t247 = t249 ^ 2;
	t262 = 0.1e1 / t295 ^ 2;
	t245 = t247 * t262 + 0.1e1;
	t243 = 0.1e1 / t245;
	t261 = 0.1e1 / t295;
	t265 = t359 * t287 + t354 * t317;
	t257 = t265 * qJD(3);
	t335 = t257 * t262;
	t273 = t300 * qJD(1);
	t274 = t277 * qJD(1);
	t357 = qJD(1) * t303 + t362 * qJD(3) - t273 * t324 - t274 * t287;
	t314 = t249 * t335 - t261 * t357;
	t213 = t314 * t243;
	t316 = t241 * t295 - t242 * t249;
	t209 = t213 * t316 + t241 * t357 + t242 * t257;
	t360 = t209 * t220;
	t345 = t219 * t360;
	t331 = 0.2e1 * (-t248 * t345 - t339 * t358) / t218 ^ 2;
	t231 = t287 * (qJD(1) * t308 - qJD(3) * t276 - t273 * t351) + t274 * t354 + t361 * qJD(3);
	t309 = t353 * t319;
	t267 = t300 * t348 + t309;
	t286 = sin(qJ(4));
	t288 = cos(qJ(4));
	t240 = t255 * t288 + t267 * t286;
	t234 = 0.1e1 / t240;
	t235 = 0.1e1 / t240 ^ 2;
	t356 = -0.2e1 * t358;
	t229 = qJD(1) * t283 + t358 * qJD(3) - t272 * t354 + t287 * t297;
	t266 = -t275 * t348 + t319 * t355;
	t258 = t266 * qJD(1);
	t223 = qJD(4) * t240 + t229 * t286 - t258 * t288;
	t239 = t255 * t286 - t267 * t288;
	t233 = t239 ^ 2;
	t227 = t233 * t235 + 0.1e1;
	t338 = t235 * t239;
	t332 = qJD(4) * t239;
	t224 = t229 * t288 + t258 * t286 - t332;
	t340 = t224 * t234 * t235;
	t344 = (t223 * t338 - t233 * t340) / t227 ^ 2;
	t334 = t261 * t335;
	t342 = (-t249 * t262 * t357 + t247 * t334) / t245 ^ 2;
	t341 = t220 * t358;
	t337 = t249 * t261;
	t336 = t249 * t265;
	t330 = -0.2e1 * t344;
	t329 = -0.2e1 * t342;
	t328 = t261 * t342;
	t327 = t239 * t340;
	t325 = t345 * t356;
	t238 = t266 * t286 + t288 * t362;
	t237 = -t266 * t288 + t286 * t362;
	t313 = -t286 * t234 + t288 * t338;
	t312 = -t261 * t362 + t262 * t336;
	t306 = -t241 + (-t242 * t337 + t241) * t243;
	t259 = -qJD(1) * t309 - t273 * t348;
	t256 = t295 * qJD(3);
	t225 = 0.1e1 / t227;
	t216 = 0.1e1 / t218;
	t214 = t312 * t243;
	t210 = t214 * t316 + t241 * t362 + t242 * t265;
	t208 = t312 * t329 + (0.2e1 * t334 * t336 + t231 * t261 + (t249 * t256 - t257 * t362 - t265 * t357) * t262) * t243;
	t1 = [-t328 * t356 + (t228 * t261 - t335 * t358) * t243, 0, t208, 0, 0, 0; t249 * t219 * t331 + (t357 * t219 + t249 * t360 + (t306 * t228 - ((t213 * t243 * t337 + t329) * t241 + (0.2e1 * t249 * t328 - t213 + (t213 - t314) * t243) * t242) * t358) * t341) * t216 - (-t341 * t331 + (-t339 + t325) * t216) * t306 * t358, 0, (-t210 * t341 - t219 * t255) * t331 + (t210 * t325 + t229 * t219 + (-t255 * t209 - t210 * t228 - (-(-t208 * t249 + t214 * t357 + t256 + (t214 * t295 + t362) * t213) * t242 - (t208 * t295 - t214 * t257 - t231 + (t214 * t249 - t265) * t213) * t241) * t358) * t220) * t216, 0, 0, 0; 0.2e1 * (-t234 * t237 + t238 * t338) * t344 + ((qJD(4) * t238 - t231 * t286 - t259 * t288) * t234 + 0.2e1 * t238 * t327 + (-t237 * t224 - (-qJD(4) * t237 - t231 * t288 + t259 * t286) * t239 - t238 * t223) * t235) * t225, 0, -t313 * t358 * t330 + (t313 * t228 - ((-qJD(4) * t234 - 0.2e1 * t327) * t288 + (t223 * t288 + (t224 - t332) * t286) * t235) * t358) * t225, t330 + 0.2e1 * (t223 * t235 * t225 + (-t225 * t340 - t235 * t344) * t239) * t239, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:34
	% DurationCPUTime: 1.97s
	% Computational Cost: add. (5146->112), mult. (15324->224), div. (448->12), fcn. (19446->15), ass. (0->109)
	t359 = sin(pkin(12));
	t364 = cos(pkin(6));
	t332 = t364 * t359;
	t362 = cos(pkin(12));
	t365 = sin(qJ(1));
	t367 = cos(qJ(1));
	t287 = t367 * t332 + t365 * t362;
	t300 = sin(qJ(3));
	t360 = sin(pkin(7));
	t361 = sin(pkin(6));
	t330 = t361 * t360;
	t323 = t367 * t330;
	t294 = t300 * t323;
	t334 = t364 * t362;
	t286 = -t367 * t334 + t365 * t359;
	t363 = cos(pkin(7));
	t338 = t286 * t363;
	t366 = cos(qJ(3));
	t374 = -t287 * t366 + t300 * t338 + t294;
	t316 = t366 * t323;
	t336 = t363 * t366;
	t373 = -t286 * t336 - t316;
	t260 = t287 * t300 - t373;
	t329 = t361 * t359;
	t331 = t363 * t361;
	t371 = t362 * t331 + t364 * t360;
	t307 = -t300 * t329 + t371 * t366;
	t257 = atan2(-t260, -t307);
	t252 = sin(t257);
	t253 = cos(t257);
	t233 = -t252 * t260 - t253 * t307;
	t231 = 0.1e1 / t233 ^ 2;
	t288 = -t365 * t332 + t367 * t362;
	t312 = t365 * t334 + t367 * t359;
	t310 = t312 * t363;
	t320 = t365 * t330;
	t315 = t366 * t320;
	t370 = -t288 * t300 - t366 * t310 + t315;
	t259 = t370 ^ 2;
	t229 = t259 * t231 + 0.1e1;
	t266 = t288 * t366 + (-t310 + t320) * t300;
	t283 = t287 * qJD(1);
	t309 = qJD(1) * t338;
	t239 = -qJD(1) * t316 + t266 * qJD(3) - t283 * t300 - t366 * t309;
	t351 = t239 * t231;
	t230 = 0.1e1 / t233;
	t258 = t260 ^ 2;
	t273 = 0.1e1 / t307 ^ 2;
	t256 = t258 * t273 + 0.1e1;
	t254 = 0.1e1 / t256;
	t272 = 0.1e1 / t307;
	t276 = t371 * t300 + t366 * t329;
	t268 = t276 * qJD(3);
	t347 = t268 * t273;
	t284 = t312 * qJD(1);
	t285 = t288 * qJD(1);
	t369 = qJD(1) * t315 + t374 * qJD(3) - t284 * t336 - t285 * t300;
	t326 = t260 * t347 - t272 * t369;
	t224 = t326 * t254;
	t328 = t252 * t307 - t253 * t260;
	t220 = t328 * t224 + t252 * t369 + t253 * t268;
	t372 = t220 * t231;
	t357 = t230 * t372;
	t343 = 0.2e1 * (-t259 * t357 - t351 * t370) / t229 ^ 2;
	t242 = (qJD(1) * t320 - qJD(3) * t287 - t363 * t284) * t300 + t285 * t366 + t373 * qJD(3);
	t321 = t365 * t331;
	t278 = t312 * t360 + t321;
	t299 = qJ(4) + pkin(13);
	t297 = sin(t299);
	t298 = cos(t299);
	t251 = t266 * t298 + t278 * t297;
	t245 = 0.1e1 / t251;
	t246 = 0.1e1 / t251 ^ 2;
	t368 = -0.2e1 * t370;
	t240 = qJD(1) * t294 + t370 * qJD(3) - t283 * t366 + t300 * t309;
	t277 = -t286 * t360 + t367 * t331;
	t270 = t277 * qJD(1);
	t234 = t251 * qJD(4) + t240 * t297 - t270 * t298;
	t250 = t266 * t297 - t278 * t298;
	t244 = t250 ^ 2;
	t238 = t244 * t246 + 0.1e1;
	t350 = t246 * t250;
	t344 = qJD(4) * t250;
	t235 = t240 * t298 + t270 * t297 - t344;
	t352 = t235 * t245 * t246;
	t356 = (t234 * t350 - t244 * t352) / t238 ^ 2;
	t346 = t272 * t347;
	t354 = (-t260 * t273 * t369 + t258 * t346) / t256 ^ 2;
	t353 = t231 * t370;
	t349 = t260 * t272;
	t348 = t260 * t276;
	t342 = -0.2e1 * t356;
	t341 = -0.2e1 * t354;
	t340 = t272 * t354;
	t339 = t250 * t352;
	t337 = t357 * t368;
	t249 = t277 * t297 + t298 * t374;
	t248 = -t277 * t298 + t297 * t374;
	t325 = -t297 * t245 + t298 * t350;
	t324 = -t272 * t374 + t273 * t348;
	t318 = -t252 + (-t253 * t349 + t252) * t254;
	t271 = -qJD(1) * t321 - t284 * t360;
	t267 = t307 * qJD(3);
	t236 = 0.1e1 / t238;
	t227 = 0.1e1 / t229;
	t225 = t324 * t254;
	t221 = t328 * t225 + t252 * t374 + t253 * t276;
	t219 = t324 * t341 + (0.2e1 * t346 * t348 + t242 * t272 + (t260 * t267 - t268 * t374 - t276 * t369) * t273) * t254;
	t1 = [-t340 * t368 + (t239 * t272 - t347 * t370) * t254, 0, t219, 0, 0, 0; t260 * t230 * t343 + (t369 * t230 + t260 * t372 + (t318 * t239 - ((t224 * t254 * t349 + t341) * t252 + (0.2e1 * t260 * t340 - t224 + (t224 - t326) * t254) * t253) * t370) * t353) * t227 - (-t353 * t343 + (-t351 + t337) * t227) * t318 * t370, 0, (-t221 * t353 - t230 * t266) * t343 + (t221 * t337 + t240 * t230 + (-t266 * t220 - t221 * t239 - (-(-t219 * t260 + t225 * t369 + t267 + (t225 * t307 + t374) * t224) * t253 - (t219 * t307 - t225 * t268 - t242 + (t225 * t260 - t276) * t224) * t252) * t370) * t231) * t227, 0, 0, 0; 0.2e1 * (-t245 * t248 + t249 * t350) * t356 + ((t249 * qJD(4) - t242 * t297 - t271 * t298) * t245 + 0.2e1 * t249 * t339 + (-t248 * t235 - (-t248 * qJD(4) - t242 * t298 + t271 * t297) * t250 - t249 * t234) * t246) * t236, 0, -t325 * t370 * t342 + (t325 * t239 - ((-qJD(4) * t245 - 0.2e1 * t339) * t298 + (t234 * t298 + (t235 - t344) * t297) * t246) * t370) * t236, t342 + 0.2e1 * (t234 * t246 * t236 + (-t236 * t352 - t246 * t356) * t250) * t250, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:33
	% EndTime: 2019-10-10 01:37:38
	% DurationCPUTime: 5.27s
	% Computational Cost: add. (17692->166), mult. (41237->320), div. (726->12), fcn. (53008->17), ass. (0->147)
	t494 = sin(pkin(12));
	t499 = cos(pkin(6));
	t454 = t499 * t494;
	t497 = cos(pkin(12));
	t500 = sin(qJ(1));
	t501 = cos(qJ(1));
	t393 = -t500 * t454 + t501 * t497;
	t390 = t393 * qJD(1);
	t402 = sin(qJ(3));
	t404 = cos(qJ(3));
	t392 = t501 * t454 + t500 * t497;
	t495 = sin(pkin(7));
	t496 = sin(pkin(6));
	t452 = t496 * t495;
	t440 = t501 * t452;
	t498 = cos(pkin(7));
	t455 = t499 * t497;
	t508 = -t501 * t455 + t500 * t494;
	t509 = t508 * t498;
	t414 = t509 + t440;
	t504 = t392 * t402 + t414 * t404;
	t423 = t500 * t455 + t501 * t494;
	t389 = t423 * qJD(1);
	t437 = t500 * t452;
	t514 = qJD(1) * t437 - t389 * t498;
	t352 = t504 * qJD(3) - t390 * t404 - t514 * t402;
	t471 = qJ(4) + pkin(13);
	t400 = sin(t471);
	t453 = t498 * t496;
	t438 = t500 * t453;
	t424 = qJD(1) * t438 + t389 * t495;
	t460 = cos(t471);
	t475 = t392 * t404;
	t373 = t414 * t402 - t475;
	t415 = -t501 * t453 + t495 * t508;
	t523 = -t373 * t400 - t415 * t460;
	t330 = qJD(4) * t523 + t352 * t460 - t424 * t400;
	t361 = t373 * t460 - t415 * t400;
	t525 = t361 * qJD(4) + t352 * t400 + t424 * t460;
	t355 = t523 ^ 2;
	t420 = t497 * t453 + t495 * t499;
	t451 = t496 * t494;
	t384 = t420 * t402 + t404 * t451;
	t391 = -t497 * t452 + t499 * t498;
	t430 = -t384 * t400 + t391 * t460;
	t365 = 0.1e1 / t430 ^ 2;
	t339 = t355 * t365 + 0.1e1;
	t333 = 0.1e1 / t339;
	t368 = t384 * t460 + t391 * t400;
	t383 = -t402 * t451 + t420 * t404;
	t378 = t383 * qJD(3);
	t353 = t368 * qJD(4) + t378 * t400;
	t364 = 0.1e1 / t430;
	t480 = t523 * t365;
	t446 = t353 * t480 - t364 * t525;
	t310 = t446 * t333;
	t340 = atan2(-t523, -t430);
	t331 = sin(t340);
	t332 = cos(t340);
	t448 = t331 * t430 - t332 * t523;
	t305 = t448 * t310 + t331 * t525 + t332 * t353;
	t322 = -t331 * t523 - t332 * t430;
	t320 = 0.1e1 / t322 ^ 2;
	t522 = t305 * t320;
	t519 = -t423 * t498 + t437;
	t375 = t393 * t404 + t519 * t402;
	t413 = t423 * t495 + t438;
	t363 = t375 * t460 + t413 * t400;
	t506 = t519 * t404;
	t374 = t393 * t402 - t506;
	t401 = sin(qJ(6));
	t403 = cos(qJ(6));
	t343 = t363 * t401 - t374 * t403;
	t518 = 0.2e1 * t343;
	t319 = 0.1e1 / t322;
	t517 = t319 * t522;
	t411 = t413 * t460;
	t362 = t375 * t400 - t411;
	t502 = 0.2e1 * t362;
	t457 = t502 * t517;
	t388 = t392 * qJD(1);
	t474 = qJD(3) * t402;
	t507 = t414 * qJD(1);
	t348 = t506 * qJD(3) - t388 * t404 - t393 * t474 + t507 * t402;
	t412 = t415 * qJD(1);
	t326 = t363 * qJD(4) + t348 * t400 + t460 * t412;
	t486 = t326 * t320;
	t513 = -t486 + t457;
	t512 = (t402 * t509 - t475) * qJD(3) - t390 * t402 + t514 * t404 + t440 * t474;
	t356 = t362 ^ 2;
	t316 = t356 * t320 + 0.1e1;
	t470 = 0.2e1 * (-t356 * t517 + t362 * t486) / t316 ^ 2;
	t511 = t353 * t365;
	t443 = -t364 * t504 + t383 * t480;
	t510 = t400 * t443;
	t449 = qJD(4) * t460;
	t344 = t363 * t403 + t374 * t401;
	t336 = 0.1e1 / t344;
	t337 = 0.1e1 / t344 ^ 2;
	t503 = -0.2e1 * t523;
	t473 = qJD(4) * t400;
	t327 = qJD(4) * t411 + t348 * t460 - t375 * t473 - t400 * t412;
	t347 = t375 * qJD(3) - t388 * t402 - t507 * t404;
	t317 = t344 * qJD(6) + t327 * t401 - t347 * t403;
	t335 = t343 ^ 2;
	t325 = t335 * t337 + 0.1e1;
	t483 = t337 * t343;
	t472 = qJD(6) * t343;
	t318 = t327 * t403 + t347 * t401 - t472;
	t489 = t318 * t336 * t337;
	t492 = (t317 * t483 - t335 * t489) / t325 ^ 2;
	t482 = t364 * t511;
	t490 = (t355 * t482 - t480 * t525) / t339 ^ 2;
	t488 = t320 * t362;
	t323 = 0.1e1 / t325;
	t487 = t323 * t337;
	t485 = t331 * t362;
	t484 = t332 * t362;
	t481 = t523 * t364;
	t478 = t374 * t400;
	t469 = -0.2e1 * t492;
	t468 = -0.2e1 * t490;
	t466 = t337 * t492;
	t465 = t364 * t490;
	t464 = t317 * t487;
	t463 = t343 * t489;
	t459 = 0.2e1 * t463;
	t458 = t482 * t503;
	t450 = t374 * t460;
	t342 = t361 * t403 - t401 * t504;
	t341 = t361 * t401 + t403 * t504;
	t445 = -t401 * t336 + t403 * t483;
	t444 = -t361 * t364 + t368 * t480;
	t436 = -t331 + (-t332 * t481 + t331) * t333;
	t435 = qJD(6) * t450 + t348;
	t425 = qJD(6) * t375 - t460 * t347 + t374 * t473;
	t379 = t384 * qJD(3);
	t354 = t430 * qJD(4) + t378 * t460;
	t346 = t375 * t401 - t403 * t450;
	t314 = 0.1e1 / t316;
	t313 = t333 * t510;
	t311 = t444 * t333;
	t307 = (t331 * t504 + t332 * t383) * t400 + t448 * t313;
	t306 = t448 * t311 + t331 * t361 + t332 * t368;
	t303 = t444 * t468 + (-t368 * t458 - t330 * t364 + (-t353 * t361 + t354 * t523 - t368 * t525) * t365) * t333;
	t302 = t468 * t510 + ((-t383 * t458 + t512 * t364 + (-t353 * t504 - t379 * t523 - t383 * t525) * t365) * t400 + t443 * t449) * t333;
	t1 = [-t465 * t502 + (t326 * t364 + t362 * t511) * t333, 0, t302, t303, 0, 0; t523 * t319 * t470 + (t525 * t319 + t523 * t522 - (t436 * t326 + ((t310 * t333 * t481 + t468) * t331 + (-t465 * t503 - t310 + (t310 - t446) * t333) * t332) * t362) * t488) * t314 + (t513 * t314 + t488 * t470) * t436 * t362, 0, (t307 * t488 + t319 * t478) * t470 + ((-t347 * t400 - t374 * t449) * t319 + t513 * t307 + (t478 * t305 - (t383 * t449 - t302 * t523 + t313 * t525 - t379 * t400 + (t313 * t430 + t400 * t504) * t310) * t484 - (t504 * t449 + t302 * t430 - t313 * t353 - t512 * t400 + (t313 * t523 - t383 * t400) * t310) * t485) * t320) * t314, (t306 * t488 - t319 * t363) * t470 + (t306 * t457 + t327 * t319 + (-t363 * t305 - t306 * t326 - (-t303 * t523 + t311 * t525 + t354 + (t311 * t430 + t361) * t310) * t484 - (t303 * t430 - t311 * t353 + t330 + (t311 * t523 - t368) * t310) * t485) * t320) * t314, 0, 0; 0.2e1 * (-t336 * t341 + t342 * t483) * t492 + ((t342 * qJD(6) + t330 * t401 - t403 * t512) * t336 + t342 * t459 + (-t341 * t318 - (-t341 * qJD(6) + t330 * t403 + t401 * t512) * t343 - t342 * t317) * t337) * t323, 0, (t466 * t518 - t464) * t346 + (-t318 * t487 + t336 * t469) * (-t375 * t403 - t401 * t450) + ((t425 * t401 - t435 * t403) * t336 - (t435 * t401 + t425 * t403) * t483 + t346 * t459) * t323, t445 * t362 * t469 + (t445 * t326 + ((-qJD(6) * t336 - 0.2e1 * t463) * t403 + (t317 * t403 + (t318 - t472) * t401) * t337) * t362) * t323, 0, t469 + (t464 + (-t323 * t489 - t466) * t343) * t518;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end