% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR11
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
%   Wie in S6RPRPRR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR11_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR11_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiaD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.32s
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
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:01
	% DurationCPUTime: 0.60s
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
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:02
	% DurationCPUTime: 1.78s
	% Computational Cost: add. (4301->104), mult. (13743->215), div. (424->12), fcn. (17471->15), ass. (0->107)
	t308 = sin(pkin(12));
	t312 = cos(pkin(6));
	t279 = t312 * t308;
	t310 = cos(pkin(12));
	t313 = sin(qJ(1));
	t315 = cos(qJ(1));
	t237 = t315 * t279 + t313 * t310;
	t250 = sin(qJ(3));
	t248 = sin(pkin(6));
	t309 = sin(pkin(7));
	t285 = t248 * t309;
	t275 = t315 * t285;
	t244 = t250 * t275;
	t281 = t312 * t310;
	t236 = -t315 * t281 + t313 * t308;
	t311 = cos(pkin(7));
	t287 = t236 * t311;
	t314 = cos(qJ(3));
	t323 = -t237 * t314 + t250 * t287 + t244;
	t238 = -t313 * t279 + t315 * t310;
	t261 = t313 * t281 + t315 * t308;
	t259 = t261 * t311;
	t273 = t313 * t285;
	t217 = t238 * t314 + (-t259 + t273) * t250;
	t286 = t248 * t311;
	t274 = t313 * t286;
	t228 = t261 * t309 + t274;
	t247 = sin(pkin(13));
	t249 = cos(pkin(13));
	t201 = t217 * t247 - t228 * t249;
	t233 = t237 * qJD(1);
	t258 = qJD(1) * t287;
	t267 = t314 * t273;
	t319 = -t238 * t250 - t314 * t259 + t267;
	t191 = qJD(1) * t244 + t319 * qJD(3) - t233 * t314 + t250 * t258;
	t227 = -t236 * t309 + t315 * t286;
	t220 = t227 * qJD(1);
	t189 = t191 * t249 + t220 * t247;
	t202 = t217 * t249 + t228 * t247;
	t196 = 0.1e1 / t202;
	t197 = 0.1e1 / t202 ^ 2;
	t302 = t189 * t196 * t197;
	t322 = 0.2e1 * t201 * t302;
	t268 = t314 * t275;
	t283 = t311 * t314;
	t321 = -t236 * t283 - t268;
	t211 = t237 * t250 - t321;
	t278 = t311 * t310;
	t280 = t312 * t309;
	t318 = (-t308 * t250 + t314 * t278) * t248 + t314 * t280;
	t208 = atan2(-t211, -t318);
	t203 = sin(t208);
	t204 = cos(t208);
	t184 = -t203 * t211 - t204 * t318;
	t182 = 0.1e1 / t184 ^ 2;
	t210 = t319 ^ 2;
	t180 = t210 * t182 + 0.1e1;
	t190 = -qJD(1) * t268 + t217 * qJD(3) - t233 * t250 - t314 * t258;
	t301 = t190 * t182;
	t181 = 0.1e1 / t184;
	t209 = t211 ^ 2;
	t223 = 0.1e1 / t318 ^ 2;
	t207 = t209 * t223 + 0.1e1;
	t205 = 0.1e1 / t207;
	t222 = 0.1e1 / t318;
	t226 = t250 * t280 + (t250 * t278 + t308 * t314) * t248;
	t219 = t226 * qJD(3);
	t297 = t219 * t223;
	t234 = t261 * qJD(1);
	t235 = t238 * qJD(1);
	t317 = qJD(1) * t267 + t323 * qJD(3) - t234 * t283 - t235 * t250;
	t272 = t211 * t297 - t222 * t317;
	t175 = t272 * t205;
	t277 = t203 * t318 - t204 * t211;
	t171 = t277 * t175 + t203 * t317 + t204 * t219;
	t320 = t171 * t182;
	t306 = t181 * t320;
	t292 = 0.2e1 * (-t210 * t306 - t301 * t319) / t180 ^ 2;
	t193 = (qJD(1) * t273 - qJD(3) * t237 - t311 * t234) * t250 + t235 * t314 + t321 * qJD(3);
	t316 = -0.2e1 * t319;
	t296 = t222 * t297;
	t304 = (-t211 * t223 * t317 + t209 * t296) / t207 ^ 2;
	t303 = t182 * t319;
	t300 = t197 * t201;
	t299 = t211 * t222;
	t298 = t211 * t226;
	t294 = t247 * t196;
	t293 = t249 * t201;
	t195 = t201 ^ 2;
	t187 = t195 * t197 + 0.1e1;
	t188 = t191 * t247 - t220 * t249;
	t291 = 0.2e1 * (t188 * t300 - t195 * t302) / t187 ^ 2;
	t290 = -0.2e1 * t304;
	t289 = t222 * t304;
	t284 = t306 * t316;
	t271 = -t222 * t323 + t223 * t298;
	t266 = -t203 + (-t204 * t299 + t203) * t205;
	t221 = -qJD(1) * t274 - t234 * t309;
	t218 = t318 * qJD(3);
	t200 = t227 * t247 + t249 * t323;
	t199 = -t227 * t249 + t247 * t323;
	t185 = 0.1e1 / t187;
	t178 = 0.1e1 / t180;
	t176 = t271 * t205;
	t172 = t277 * t176 + t203 * t323 + t204 * t226;
	t170 = t271 * t290 + (0.2e1 * t296 * t298 + t193 * t222 + (t211 * t218 - t219 * t323 - t226 * t317) * t223) * t205;
	t1 = [-t289 * t316 + (t190 * t222 - t297 * t319) * t205, 0, t170, 0, 0, 0; t211 * t181 * t292 + (t317 * t181 + t211 * t320 + (t266 * t190 - ((t175 * t205 * t299 + t290) * t203 + (0.2e1 * t211 * t289 - t175 + (t175 - t272) * t205) * t204) * t319) * t303) * t178 - (-t303 * t292 + (-t301 + t284) * t178) * t266 * t319, 0, (-t172 * t303 - t181 * t217) * t292 + (t172 * t284 + t191 * t181 + (-t217 * t171 - t172 * t190 - (-(-t170 * t211 + t176 * t317 + t218 + (t176 * t318 + t323) * t175) * t204 - (t170 * t318 - t176 * t219 - t193 + (t176 * t211 - t226) * t175) * t203) * t319) * t182) * t178, 0, 0, 0; (-t196 * t199 + t200 * t300) * t291 + ((-t193 * t247 - t221 * t249) * t196 + t200 * t322 + (-t199 * t189 - (-t193 * t249 + t221 * t247) * t201 - t200 * t188) * t197) * t185, 0, -(-t197 * t293 + t294) * t319 * t291 + (t319 * t249 * t322 - t190 * t294 + (t190 * t293 - (t188 * t249 + t189 * t247) * t319) * t197) * t185, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:01
	% EndTime: 2019-10-10 01:04:03
	% DurationCPUTime: 2.06s
	% Computational Cost: add. (5146->112), mult. (15324->224), div. (448->12), fcn. (19446->15), ass. (0->109)
	t356 = sin(pkin(12));
	t361 = cos(pkin(6));
	t329 = t361 * t356;
	t359 = cos(pkin(12));
	t362 = sin(qJ(1));
	t364 = cos(qJ(1));
	t284 = t364 * t329 + t362 * t359;
	t297 = sin(qJ(3));
	t357 = sin(pkin(7));
	t358 = sin(pkin(6));
	t327 = t358 * t357;
	t320 = t364 * t327;
	t291 = t297 * t320;
	t331 = t361 * t359;
	t283 = -t364 * t331 + t362 * t356;
	t360 = cos(pkin(7));
	t335 = t283 * t360;
	t363 = cos(qJ(3));
	t371 = -t284 * t363 + t297 * t335 + t291;
	t313 = t363 * t320;
	t333 = t360 * t363;
	t370 = -t283 * t333 - t313;
	t257 = t284 * t297 - t370;
	t326 = t358 * t356;
	t328 = t360 * t358;
	t368 = t359 * t328 + t361 * t357;
	t304 = -t297 * t326 + t368 * t363;
	t254 = atan2(-t257, -t304);
	t249 = sin(t254);
	t250 = cos(t254);
	t230 = -t249 * t257 - t250 * t304;
	t228 = 0.1e1 / t230 ^ 2;
	t285 = -t362 * t329 + t364 * t359;
	t309 = t362 * t331 + t364 * t356;
	t307 = t309 * t360;
	t317 = t362 * t327;
	t312 = t363 * t317;
	t367 = -t285 * t297 - t363 * t307 + t312;
	t256 = t367 ^ 2;
	t226 = t256 * t228 + 0.1e1;
	t263 = t285 * t363 + (-t307 + t317) * t297;
	t280 = t284 * qJD(1);
	t306 = qJD(1) * t335;
	t236 = -qJD(1) * t313 + t263 * qJD(3) - t280 * t297 - t363 * t306;
	t348 = t236 * t228;
	t227 = 0.1e1 / t230;
	t255 = t257 ^ 2;
	t270 = 0.1e1 / t304 ^ 2;
	t253 = t255 * t270 + 0.1e1;
	t251 = 0.1e1 / t253;
	t269 = 0.1e1 / t304;
	t273 = t368 * t297 + t363 * t326;
	t265 = t273 * qJD(3);
	t344 = t265 * t270;
	t281 = t309 * qJD(1);
	t282 = t285 * qJD(1);
	t366 = qJD(1) * t312 + t371 * qJD(3) - t281 * t333 - t282 * t297;
	t323 = t257 * t344 - t269 * t366;
	t221 = t323 * t251;
	t325 = t249 * t304 - t250 * t257;
	t217 = t221 * t325 + t249 * t366 + t250 * t265;
	t369 = t217 * t228;
	t354 = t227 * t369;
	t340 = 0.2e1 * (-t256 * t354 - t348 * t367) / t226 ^ 2;
	t239 = t297 * (qJD(1) * t317 - qJD(3) * t284 - t360 * t281) + t282 * t363 + t370 * qJD(3);
	t318 = t362 * t328;
	t275 = t309 * t357 + t318;
	t296 = pkin(13) + qJ(5);
	t294 = sin(t296);
	t295 = cos(t296);
	t248 = t263 * t295 + t275 * t294;
	t242 = 0.1e1 / t248;
	t243 = 0.1e1 / t248 ^ 2;
	t365 = -0.2e1 * t367;
	t237 = qJD(1) * t291 + t367 * qJD(3) - t280 * t363 + t297 * t306;
	t274 = -t283 * t357 + t364 * t328;
	t267 = t274 * qJD(1);
	t231 = qJD(5) * t248 + t237 * t294 - t267 * t295;
	t247 = t263 * t294 - t275 * t295;
	t241 = t247 ^ 2;
	t235 = t241 * t243 + 0.1e1;
	t347 = t243 * t247;
	t341 = qJD(5) * t247;
	t232 = t237 * t295 + t267 * t294 - t341;
	t349 = t232 * t242 * t243;
	t353 = (t231 * t347 - t241 * t349) / t235 ^ 2;
	t343 = t269 * t344;
	t351 = (-t257 * t270 * t366 + t255 * t343) / t253 ^ 2;
	t350 = t228 * t367;
	t346 = t257 * t269;
	t345 = t257 * t273;
	t339 = -0.2e1 * t353;
	t338 = -0.2e1 * t351;
	t337 = t269 * t351;
	t336 = t247 * t349;
	t334 = t354 * t365;
	t246 = t274 * t294 + t295 * t371;
	t245 = -t274 * t295 + t294 * t371;
	t322 = -t294 * t242 + t295 * t347;
	t321 = -t269 * t371 + t270 * t345;
	t315 = -t249 + (-t250 * t346 + t249) * t251;
	t268 = -qJD(1) * t318 - t281 * t357;
	t264 = t304 * qJD(3);
	t233 = 0.1e1 / t235;
	t224 = 0.1e1 / t226;
	t222 = t321 * t251;
	t218 = t222 * t325 + t249 * t371 + t250 * t273;
	t216 = t321 * t338 + (0.2e1 * t343 * t345 + t239 * t269 + (t257 * t264 - t265 * t371 - t273 * t366) * t270) * t251;
	t1 = [-t337 * t365 + (t236 * t269 - t344 * t367) * t251, 0, t216, 0, 0, 0; t257 * t227 * t340 + (t366 * t227 + t257 * t369 + (t315 * t236 - ((t221 * t251 * t346 + t338) * t249 + (0.2e1 * t257 * t337 - t221 + (t221 - t323) * t251) * t250) * t367) * t350) * t224 - (-t350 * t340 + (-t348 + t334) * t224) * t315 * t367, 0, (-t218 * t350 - t227 * t263) * t340 + (t218 * t334 + t237 * t227 + (-t263 * t217 - t218 * t236 - (-(-t216 * t257 + t222 * t366 + t264 + (t222 * t304 + t371) * t221) * t250 - (t216 * t304 - t222 * t265 - t239 + (t222 * t257 - t273) * t221) * t249) * t367) * t228) * t224, 0, 0, 0; 0.2e1 * (-t242 * t245 + t246 * t347) * t353 + ((qJD(5) * t246 - t239 * t294 - t268 * t295) * t242 + 0.2e1 * t246 * t336 + (-t245 * t232 - (-qJD(5) * t245 - t239 * t295 + t268 * t294) * t247 - t246 * t231) * t243) * t233, 0, -t322 * t367 * t339 + (t322 * t236 - ((-qJD(5) * t242 - 0.2e1 * t336) * t295 + (t231 * t295 + (t232 - t341) * t294) * t243) * t367) * t233, 0, t339 + 0.2e1 * (t231 * t243 * t233 + (-t233 * t349 - t243 * t353) * t247) * t247, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:01
	% EndTime: 2019-10-10 01:04:06
	% DurationCPUTime: 5.26s
	% Computational Cost: add. (17692->166), mult. (41237->320), div. (726->12), fcn. (53008->17), ass. (0->147)
	t495 = sin(pkin(12));
	t500 = cos(pkin(6));
	t455 = t500 * t495;
	t498 = cos(pkin(12));
	t501 = sin(qJ(1));
	t502 = cos(qJ(1));
	t394 = -t501 * t455 + t502 * t498;
	t391 = t394 * qJD(1);
	t403 = sin(qJ(3));
	t405 = cos(qJ(3));
	t393 = t502 * t455 + t501 * t498;
	t496 = sin(pkin(7));
	t497 = sin(pkin(6));
	t453 = t497 * t496;
	t441 = t502 * t453;
	t499 = cos(pkin(7));
	t456 = t500 * t498;
	t509 = -t502 * t456 + t501 * t495;
	t510 = t509 * t499;
	t415 = t510 + t441;
	t505 = t393 * t403 + t415 * t405;
	t424 = t501 * t456 + t502 * t495;
	t390 = t424 * qJD(1);
	t438 = t501 * t453;
	t515 = qJD(1) * t438 - t390 * t499;
	t353 = t505 * qJD(3) - t391 * t405 - t515 * t403;
	t472 = pkin(13) + qJ(5);
	t401 = sin(t472);
	t454 = t499 * t497;
	t439 = t501 * t454;
	t425 = qJD(1) * t439 + t390 * t496;
	t461 = cos(t472);
	t476 = t393 * t405;
	t374 = t415 * t403 - t476;
	t416 = -t502 * t454 + t496 * t509;
	t524 = -t374 * t401 - t416 * t461;
	t331 = qJD(5) * t524 + t353 * t461 - t425 * t401;
	t362 = t374 * t461 - t416 * t401;
	t526 = t362 * qJD(5) + t353 * t401 + t425 * t461;
	t356 = t524 ^ 2;
	t421 = t498 * t454 + t500 * t496;
	t452 = t497 * t495;
	t385 = t421 * t403 + t405 * t452;
	t392 = -t498 * t453 + t500 * t499;
	t431 = -t385 * t401 + t392 * t461;
	t366 = 0.1e1 / t431 ^ 2;
	t340 = t356 * t366 + 0.1e1;
	t334 = 0.1e1 / t340;
	t369 = t385 * t461 + t392 * t401;
	t384 = -t403 * t452 + t421 * t405;
	t379 = t384 * qJD(3);
	t354 = t369 * qJD(5) + t379 * t401;
	t365 = 0.1e1 / t431;
	t481 = t524 * t366;
	t447 = t354 * t481 - t365 * t526;
	t311 = t447 * t334;
	t341 = atan2(-t524, -t431);
	t332 = sin(t341);
	t333 = cos(t341);
	t449 = t332 * t431 - t333 * t524;
	t306 = t449 * t311 + t332 * t526 + t333 * t354;
	t323 = -t332 * t524 - t333 * t431;
	t321 = 0.1e1 / t323 ^ 2;
	t523 = t306 * t321;
	t520 = -t424 * t499 + t438;
	t376 = t394 * t405 + t520 * t403;
	t414 = t424 * t496 + t439;
	t364 = t376 * t461 + t414 * t401;
	t507 = t520 * t405;
	t375 = t394 * t403 - t507;
	t402 = sin(qJ(6));
	t404 = cos(qJ(6));
	t344 = t364 * t402 - t375 * t404;
	t519 = 0.2e1 * t344;
	t320 = 0.1e1 / t323;
	t518 = t320 * t523;
	t412 = t414 * t461;
	t363 = t376 * t401 - t412;
	t503 = 0.2e1 * t363;
	t458 = t503 * t518;
	t389 = t393 * qJD(1);
	t475 = qJD(3) * t403;
	t508 = t415 * qJD(1);
	t349 = t507 * qJD(3) - t389 * t405 - t394 * t475 + t508 * t403;
	t413 = t416 * qJD(1);
	t327 = t364 * qJD(5) + t349 * t401 + t461 * t413;
	t487 = t327 * t321;
	t514 = -t487 + t458;
	t513 = (t403 * t510 - t476) * qJD(3) - t391 * t403 + t515 * t405 + t441 * t475;
	t357 = t363 ^ 2;
	t317 = t357 * t321 + 0.1e1;
	t471 = 0.2e1 * (-t357 * t518 + t363 * t487) / t317 ^ 2;
	t512 = t354 * t366;
	t444 = -t365 * t505 + t384 * t481;
	t511 = t401 * t444;
	t450 = qJD(5) * t461;
	t345 = t364 * t404 + t375 * t402;
	t337 = 0.1e1 / t345;
	t338 = 0.1e1 / t345 ^ 2;
	t504 = -0.2e1 * t524;
	t474 = qJD(5) * t401;
	t328 = qJD(5) * t412 + t349 * t461 - t376 * t474 - t401 * t413;
	t348 = t376 * qJD(3) - t389 * t403 - t508 * t405;
	t318 = t345 * qJD(6) + t328 * t402 - t348 * t404;
	t336 = t344 ^ 2;
	t326 = t336 * t338 + 0.1e1;
	t484 = t338 * t344;
	t473 = qJD(6) * t344;
	t319 = t328 * t404 + t348 * t402 - t473;
	t490 = t319 * t337 * t338;
	t493 = (t318 * t484 - t336 * t490) / t326 ^ 2;
	t483 = t365 * t512;
	t491 = (t356 * t483 - t481 * t526) / t340 ^ 2;
	t489 = t321 * t363;
	t324 = 0.1e1 / t326;
	t488 = t324 * t338;
	t486 = t332 * t363;
	t485 = t333 * t363;
	t482 = t524 * t365;
	t479 = t375 * t401;
	t470 = -0.2e1 * t493;
	t469 = -0.2e1 * t491;
	t467 = t338 * t493;
	t466 = t365 * t491;
	t465 = t318 * t488;
	t464 = t344 * t490;
	t460 = 0.2e1 * t464;
	t459 = t483 * t504;
	t451 = t375 * t461;
	t343 = t362 * t404 - t402 * t505;
	t342 = t362 * t402 + t404 * t505;
	t446 = -t402 * t337 + t404 * t484;
	t445 = -t362 * t365 + t369 * t481;
	t437 = -t332 + (-t333 * t482 + t332) * t334;
	t436 = qJD(6) * t451 + t349;
	t426 = qJD(6) * t376 - t461 * t348 + t375 * t474;
	t380 = t385 * qJD(3);
	t355 = t431 * qJD(5) + t379 * t461;
	t347 = t376 * t402 - t404 * t451;
	t315 = 0.1e1 / t317;
	t314 = t334 * t511;
	t312 = t445 * t334;
	t308 = (t332 * t505 + t333 * t384) * t401 + t449 * t314;
	t307 = t449 * t312 + t332 * t362 + t333 * t369;
	t304 = t445 * t469 + (-t369 * t459 - t331 * t365 + (-t354 * t362 + t355 * t524 - t369 * t526) * t366) * t334;
	t303 = t469 * t511 + ((-t384 * t459 + t513 * t365 + (-t354 * t505 - t380 * t524 - t384 * t526) * t366) * t401 + t444 * t450) * t334;
	t1 = [-t466 * t503 + (t327 * t365 + t363 * t512) * t334, 0, t303, 0, t304, 0; t524 * t320 * t471 + (t526 * t320 + t524 * t523 - (t437 * t327 + ((t311 * t334 * t482 + t469) * t332 + (-t466 * t504 - t311 + (t311 - t447) * t334) * t333) * t363) * t489) * t315 + (t514 * t315 + t489 * t471) * t437 * t363, 0, (t308 * t489 + t320 * t479) * t471 + ((-t348 * t401 - t375 * t450) * t320 + t514 * t308 + (t479 * t306 - (t384 * t450 - t303 * t524 + t314 * t526 - t380 * t401 + (t314 * t431 + t401 * t505) * t311) * t485 - (t505 * t450 + t303 * t431 - t314 * t354 - t513 * t401 + (t314 * t524 - t384 * t401) * t311) * t486) * t321) * t315, 0, (t307 * t489 - t320 * t364) * t471 + (t307 * t458 + t328 * t320 + (-t364 * t306 - t307 * t327 - (-t304 * t524 + t312 * t526 + t355 + (t312 * t431 + t362) * t311) * t485 - (t304 * t431 - t312 * t354 + t331 + (t312 * t524 - t369) * t311) * t486) * t321) * t315, 0; 0.2e1 * (-t337 * t342 + t343 * t484) * t493 + ((t343 * qJD(6) + t331 * t402 - t404 * t513) * t337 + t343 * t460 + (-t342 * t319 - (-t342 * qJD(6) + t331 * t404 + t402 * t513) * t344 - t343 * t318) * t338) * t324, 0, (t467 * t519 - t465) * t347 + (-t319 * t488 + t337 * t470) * (-t376 * t404 - t402 * t451) + ((t426 * t402 - t436 * t404) * t337 - (t436 * t402 + t426 * t404) * t484 + t347 * t460) * t324, 0, t446 * t363 * t470 + (t446 * t327 + ((-qJD(6) * t337 - 0.2e1 * t464) * t404 + (t318 * t404 + (t319 - t473) * t402) * t338) * t363) * t324, t470 + (t465 + (-t324 * t490 - t467) * t344) * t519;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end