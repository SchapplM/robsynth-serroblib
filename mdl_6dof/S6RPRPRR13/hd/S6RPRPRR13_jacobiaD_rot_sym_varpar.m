% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR13
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
%   Wie in S6RPRPRR13_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR13_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR13_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
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
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:45
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
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:46
	% DurationCPUTime: 1.58s
	% Computational Cost: add. (3797->92), mult. (12220->188), div. (410->12), fcn. (15578->13), ass. (0->99)
	t222 = sin(pkin(7));
	t225 = cos(pkin(7));
	t226 = cos(pkin(6));
	t228 = cos(qJ(1));
	t289 = sin(pkin(12));
	t259 = t228 * t289;
	t224 = cos(pkin(12));
	t290 = sin(qJ(1));
	t260 = t290 * t224;
	t241 = t226 * t260 + t259;
	t223 = sin(pkin(6));
	t262 = t223 * t290;
	t203 = t241 * t222 + t225 * t262;
	t200 = 0.1e1 / t203 ^ 2;
	t253 = t290 * t289;
	t270 = t228 * t224;
	t244 = t226 * t270 - t253;
	t272 = t223 * t228;
	t202 = t244 * t222 + t225 * t272;
	t299 = t200 * t202;
	t227 = sin(qJ(3));
	t273 = t222 * t227;
	t219 = t272 * t273;
	t243 = -t226 * t259 - t260;
	t271 = t225 * t227;
	t291 = cos(qJ(3));
	t298 = t243 * t291 - t244 * t271 + t219;
	t263 = t222 * t291;
	t252 = t263 * t272;
	t261 = t225 * t291;
	t297 = t244 * t261 - t252;
	t184 = -t227 * t243 - t297;
	t294 = (t224 * t261 - t289 * t227) * t223 + t226 * t263;
	t177 = atan2(-t184, -t294);
	t172 = sin(t177);
	t173 = cos(t177);
	t167 = -t172 * t184 - t173 * t294;
	t165 = 0.1e1 / t167 ^ 2;
	t213 = -t226 * t253 + t270;
	t239 = t241 * t225;
	t255 = t222 * t262;
	t247 = t291 * t255;
	t295 = -t213 * t227 - t291 * t239 + t247;
	t182 = t295 ^ 2;
	t163 = t182 * t165 + 0.1e1;
	t190 = t213 * t291 + (-t239 + t255) * t227;
	t208 = t243 * qJD(1);
	t236 = qJD(1) * t225 * t244;
	t168 = -qJD(1) * t252 + t190 * qJD(3) + t208 * t227 + t291 * t236;
	t282 = t168 * t165;
	t164 = 0.1e1 / t167;
	t181 = t184 ^ 2;
	t195 = 0.1e1 / t294 ^ 2;
	t176 = t181 * t195 + 0.1e1;
	t174 = 0.1e1 / t176;
	t194 = 0.1e1 / t294;
	t198 = t226 * t273 + (t224 * t271 + t289 * t291) * t223;
	t192 = t198 * qJD(3);
	t278 = t192 * t195;
	t209 = t241 * qJD(1);
	t210 = t213 * qJD(1);
	t293 = qJD(1) * t247 + t298 * qJD(3) - t209 * t261 - t210 * t227;
	t249 = t184 * t278 - t194 * t293;
	t157 = t249 * t174;
	t251 = t172 * t294 - t173 * t184;
	t154 = t251 * t157 + t172 * t293 + t173 * t192;
	t296 = t154 * t165;
	t287 = t164 * t296;
	t269 = 0.2e1 * (-t182 * t287 - t282 * t295) / t163 ^ 2;
	t276 = qJD(1) * t299;
	t254 = qJD(1) * t262;
	t171 = (qJD(3) * t243 - t209 * t225 + t222 * t254) * t227 + t210 * t291 + t297 * qJD(3);
	t199 = 0.1e1 / t203;
	t292 = -0.2e1 * t295;
	t277 = t194 * t278;
	t285 = (-t184 * t195 * t293 + t181 * t277) / t176 ^ 2;
	t169 = qJD(1) * t219 + t295 * qJD(3) + t208 * t291 - t227 * t236;
	t183 = t190 ^ 2;
	t180 = t183 * t200 + 0.1e1;
	t275 = t199 * t276;
	t279 = t190 * t200;
	t284 = (t169 * t279 - t183 * t275) / t180 ^ 2;
	t283 = t165 * t295;
	t281 = t184 * t194;
	t280 = t184 * t198;
	t268 = -0.2e1 * t285;
	t267 = 0.2e1 * t190 * t202;
	t266 = t194 * t285;
	t265 = t199 * t284;
	t258 = t287 * t292;
	t248 = -t194 * t298 + t195 * t280;
	t245 = -t172 + (-t173 * t281 + t172) * t174;
	t191 = t294 * qJD(3);
	t178 = 0.1e1 / t180;
	t161 = 0.1e1 / t163;
	t158 = t248 * t174;
	t155 = t251 * t158 + t172 * t298 + t173 * t198;
	t153 = t248 * t268 + (0.2e1 * t277 * t280 + t171 * t194 + (t184 * t191 - t192 * t298 - t198 * t293) * t195) * t174;
	t1 = [-t266 * t292 + (t168 * t194 - t278 * t295) * t174, 0, t153, 0, 0, 0; t184 * t164 * t269 + (t293 * t164 + t184 * t296 + (t245 * t168 - ((t157 * t174 * t281 + t268) * t172 + (0.2e1 * t184 * t266 - t157 + (t157 - t249) * t174) * t173) * t295) * t283) * t161 - (-t283 * t269 + (-t282 + t258) * t161) * t245 * t295, 0, (-t155 * t283 - t164 * t190) * t269 + (t155 * t258 + t169 * t164 + (-t190 * t154 - t155 * t168 - (-(-t153 * t184 + t158 * t293 + t191 + (t158 * t294 + t298) * t157) * t173 - (t153 * t294 - t158 * t192 - t171 + (t158 * t184 - t198) * t157) * t172) * t295) * t165) * t161, 0, 0, 0; t200 * t267 * t284 - 0.2e1 * t298 * t265 + (-t171 * t199 - t298 * t276 - (-t209 * t222 - t225 * t254) * t279 - t169 * t299 + t267 * t275) * t178, 0, t265 * t292 + (-t168 * t199 - t276 * t295) * t178, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:47
	% DurationCPUTime: 1.76s
	% Computational Cost: add. (4837->121), mult. (15324->242), div. (448->12), fcn. (19446->15), ass. (0->119)
	t291 = cos(pkin(6));
	t286 = sin(pkin(12));
	t294 = sin(qJ(1));
	t337 = t294 * t286;
	t322 = t291 * t337;
	t289 = cos(pkin(12));
	t297 = cos(qJ(1));
	t334 = t297 * t289;
	t278 = -t322 + t334;
	t293 = sin(qJ(3));
	t296 = cos(qJ(3));
	t335 = t297 * t286;
	t336 = t294 * t289;
	t307 = t291 * t336 + t335;
	t287 = sin(pkin(7));
	t288 = sin(pkin(6));
	t341 = t288 * t294;
	t325 = t287 * t341;
	t290 = cos(pkin(7));
	t338 = t290 * t296;
	t254 = t278 * t293 - t296 * t325 + t307 * t338;
	t266 = t287 * t307 + t290 * t341;
	t292 = sin(qJ(5));
	t295 = cos(qJ(5));
	t314 = t254 * t295 - t266 * t292;
	t359 = t314 * qJD(5);
	t339 = t290 * t293;
	t342 = t287 * t291;
	t264 = (t286 * t296 + t289 * t339) * t288 + t293 * t342;
	t260 = 0.1e1 / t264;
	t273 = t307 * qJD(1);
	t333 = qJD(1) * t297;
	t274 = -qJD(1) * t322 + t289 * t333;
	t276 = t291 * t335 + t336;
	t340 = t288 * t297;
	t324 = t287 * t340;
	t317 = qJD(3) * t324;
	t280 = t296 * t317;
	t321 = qJD(1) * t341;
	t318 = t287 * t321;
	t331 = qJD(3) * t293;
	t275 = -t291 * t334 + t337;
	t332 = qJD(3) * t275;
	t230 = t293 * t318 + t274 * t296 - t276 * t331 - t280 + (-t273 * t293 - t296 * t332) * t290;
	t310 = t275 * t339 - t276 * t296;
	t253 = t293 * t324 + t310;
	t247 = t253 ^ 2;
	t261 = 0.1e1 / t264 ^ 2;
	t244 = t247 * t261 + 0.1e1;
	t263 = t296 * t342 + (-t286 * t293 + t289 * t338) * t288;
	t256 = t263 * qJD(3);
	t345 = t256 * t261;
	t344 = t260 * t345;
	t353 = (-t230 * t253 * t261 - t247 * t344) / t244 ^ 2;
	t358 = 0.2e1 * t260 * t353;
	t231 = t310 * qJD(3) - t273 * t338 + t296 * t318 + (-t274 + t317) * t293;
	t245 = atan2(t253, t264);
	t240 = sin(t245);
	t241 = cos(t245);
	t221 = t240 * t253 + t241 * t264;
	t218 = 0.1e1 / t221;
	t239 = t254 * t292 + t266 * t295;
	t233 = 0.1e1 / t239;
	t219 = 0.1e1 / t221 ^ 2;
	t234 = 0.1e1 / t239 ^ 2;
	t309 = -t290 * t307 + t325;
	t255 = t278 * t296 + t309 * t293;
	t248 = t255 ^ 2;
	t217 = t248 * t219 + 0.1e1;
	t272 = t276 * qJD(1);
	t271 = t275 * qJD(1);
	t320 = t288 * t333;
	t306 = t271 * t290 + t287 * t320;
	t228 = -t278 * t331 + t306 * t293 + (t309 * qJD(3) - t272) * t296;
	t351 = t228 * t219;
	t242 = 0.1e1 / t244;
	t313 = -t230 * t260 - t253 * t345;
	t212 = t313 * t242;
	t316 = -t240 * t264 + t241 * t253;
	t208 = t316 * t212 - t240 * t230 + t241 * t256;
	t356 = t208 * t218 * t219;
	t357 = (-t248 * t356 + t255 * t351) / t217 ^ 2;
	t343 = t276 * t293;
	t249 = t275 * t338 + t296 * t324 + t343;
	t346 = t253 * t263;
	t311 = t249 * t260 - t261 * t346;
	t213 = t311 * t242;
	t209 = t316 * t213 + t240 * t249 + t241 * t263;
	t355 = t209 * t255;
	t227 = t255 * qJD(3) - t272 * t293 - t306 * t296;
	t258 = -t271 * t287 + t290 * t320;
	t222 = t239 * qJD(5) - t227 * t295 + t258 * t292;
	t232 = t314 ^ 2;
	t226 = t232 * t234 + 0.1e1;
	t350 = t234 * t314;
	t223 = t227 * t292 + t258 * t295 + t359;
	t352 = t223 * t233 * t234;
	t354 = (-t222 * t350 - t232 * t352) / t226 ^ 2;
	t349 = t240 * t255;
	t348 = t241 * t255;
	t347 = t253 * t260;
	t330 = 0.2e1 * t357;
	t329 = 0.2e1 * t356;
	t328 = 0.2e1 * t354;
	t327 = -0.2e1 * t353;
	t319 = -0.2e1 * t314 * t352;
	t252 = -t343 + (-t275 * t290 - t324) * t296;
	t265 = -t275 * t287 + t290 * t340;
	t315 = t252 * t295 - t265 * t292;
	t237 = t252 * t292 + t265 * t295;
	t312 = t295 * t233 - t292 * t350;
	t305 = -t240 + (-t241 * t347 + t240) * t242;
	t259 = -t273 * t287 - t290 * t321;
	t257 = t264 * qJD(3);
	t224 = 0.1e1 / t226;
	t215 = 0.1e1 / t217;
	t211 = t305 * t255;
	t207 = t311 * t327 + (0.2e1 * t344 * t346 - t231 * t260 + (t230 * t263 - t249 * t256 + t253 * t257) * t261) * t242;
	t1 = [t255 * t358 + (-t228 * t260 + t255 * t345) * t242, 0, t207, 0, 0, 0; -0.2e1 * t253 * t218 * t357 + ((t280 + (t290 * t332 - t274) * t296 + (qJD(3) * t276 + t273 * t290 - t318) * t293) * t218 + (-t253 * t208 - t211 * t228) * t219) * t215 + ((t211 * t329 - t305 * t351) * t215 + (t211 * t330 + (-(t212 * t242 * t347 + t327) * t349 - (t253 * t358 - t212 + (t212 - t313) * t242) * t348) * t215) * t219) * t255, 0, (t218 * t254 + t219 * t355) * t330 + (t329 * t355 - t227 * t218 + (t254 * t208 - t209 * t228 - (t207 * t253 - t213 * t230 - t257 + (-t213 * t264 + t249) * t212) * t348 - (-t207 * t264 - t213 * t256 - t231 + (-t213 * t253 - t263) * t212) * t349) * t219) * t215, 0, 0, 0; (t233 * t315 - t237 * t350) * t328 + ((t237 * qJD(5) - t231 * t295 + t259 * t292) * t233 + t237 * t319 + (t315 * t223 + (t315 * qJD(5) + t231 * t292 + t259 * t295) * t314 - t237 * t222) * t234) * t224, 0, t312 * t255 * t328 + (-t312 * t228 + ((qJD(5) * t233 + t319) * t292 + (-t222 * t292 + (t223 + t359) * t295) * t234) * t255) * t224, 0, -0.2e1 * t354 - 0.2e1 * (t222 * t234 * t224 - (-t224 * t352 - t234 * t354) * t314) * t314, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:46
	% EndTime: 2019-10-10 01:07:51
	% DurationCPUTime: 5.23s
	% Computational Cost: add. (13786->169), mult. (41237->335), div. (726->12), fcn. (53008->17), ass. (0->145)
	t389 = sin(qJ(3));
	t479 = sin(pkin(12));
	t483 = cos(pkin(6));
	t437 = t483 * t479;
	t481 = cos(pkin(12));
	t484 = sin(qJ(1));
	t486 = cos(qJ(1));
	t412 = t486 * t437 + t484 * t481;
	t485 = cos(qJ(3));
	t407 = t412 * t485;
	t386 = sin(pkin(6));
	t480 = sin(pkin(7));
	t442 = t386 * t480;
	t431 = t486 * t442;
	t482 = cos(pkin(7));
	t439 = t483 * t481;
	t494 = -t486 * t439 + t484 * t479;
	t495 = t494 * t482;
	t365 = -t407 + (t495 + t431) * t389;
	t414 = -t484 * t437 + t486 * t481;
	t380 = t414 * qJD(1);
	t429 = t484 * t442;
	t383 = t485 * t429;
	t413 = t484 * t439 + t486 * t479;
	t379 = t413 * qJD(1);
	t444 = t379 * t482;
	t340 = qJD(1) * t383 + t365 * qJD(3) - t380 * t389 - t485 * t444;
	t443 = t386 * t482;
	t430 = t484 * t443;
	t371 = qJD(1) * t430 + t379 * t480;
	t388 = sin(qJ(5));
	t391 = cos(qJ(5));
	t375 = -t486 * t443 + t480 * t494;
	t499 = -t412 * t389 - t485 * t431;
	t400 = t485 * t495 - t499;
	t505 = t375 * t388 - t400 * t391;
	t319 = qJD(5) * t505 + t340 * t388 - t371 * t391;
	t352 = t375 * t391 + t400 * t388;
	t320 = t352 * qJD(5) + t340 * t391 + t371 * t388;
	t341 = t400 * qJD(3) - (qJD(1) * t429 - t444) * t389 - t380 * t485;
	t376 = t413 * t480 + t430;
	t406 = t413 * t482;
	t404 = t414 * t389 + t485 * t406 - t383;
	t354 = t376 * t391 + t404 * t388;
	t366 = t414 * t485 + (-t406 + t429) * t389;
	t387 = sin(qJ(6));
	t390 = cos(qJ(6));
	t328 = t354 * t387 - t366 * t390;
	t502 = 0.2e1 * t328;
	t381 = -t481 * t442 + t483 * t482;
	t436 = t482 * t481;
	t438 = t483 * t480;
	t489 = (-t479 * t389 + t485 * t436) * t386 + t485 * t438;
	t433 = -t381 * t388 - t391 * t489;
	t335 = atan2(-t505, -t433);
	t330 = sin(t335);
	t331 = cos(t335);
	t313 = -t330 * t505 - t331 * t433;
	t311 = 0.1e1 / t313 ^ 2;
	t402 = t404 * t391;
	t353 = t376 * t388 - t402;
	t347 = t353 ^ 2;
	t307 = t347 * t311 + 0.1e1;
	t369 = t375 * qJD(1);
	t405 = qJD(1) * t495;
	t398 = t499 * qJD(1) + t366 * qJD(3) - t485 * t405;
	t317 = t354 * qJD(5) - t369 * t388 - t398 * t391;
	t471 = t317 * t311;
	t346 = t505 ^ 2;
	t359 = 0.1e1 / t433 ^ 2;
	t334 = t346 * t359 + 0.1e1;
	t332 = 0.1e1 / t334;
	t362 = t381 * t391 - t388 * t489;
	t373 = t389 * t438 + (t389 * t436 + t479 * t485) * t386;
	t368 = t373 * qJD(3);
	t345 = t362 * qJD(5) - t368 * t391;
	t358 = 0.1e1 / t433;
	t465 = t505 * t359;
	t428 = t320 * t358 + t345 * t465;
	t301 = t428 * t332;
	t435 = t330 * t433 - t331 * t505;
	t296 = t435 * t301 - t330 * t320 + t331 * t345;
	t310 = 0.1e1 / t313;
	t312 = t310 * t311;
	t477 = t296 * t312;
	t458 = 0.2e1 * (-t347 * t477 + t353 * t471) / t307 ^ 2;
	t497 = t345 * t359;
	t425 = -t358 * t365 + t373 * t465;
	t496 = t391 * t425;
	t493 = t366 * t388 * qJD(6) + t398;
	t464 = t366 * t391;
	t491 = qJD(5) * t464 - t404 * qJD(6);
	t329 = t354 * t390 + t366 * t387;
	t323 = 0.1e1 / t329;
	t324 = 0.1e1 / t329 ^ 2;
	t487 = 0.2e1 * t353;
	t460 = qJD(5) * t388;
	t318 = qJD(5) * t402 - t369 * t391 - t376 * t460 + t398 * t388;
	t338 = t389 * t405 + (t389 * t431 - t407) * qJD(1) - t404 * qJD(3);
	t308 = t329 * qJD(6) + t318 * t387 - t338 * t390;
	t322 = t328 ^ 2;
	t316 = t322 * t324 + 0.1e1;
	t470 = t324 * t328;
	t459 = qJD(6) * t328;
	t309 = t318 * t390 + t338 * t387 - t459;
	t474 = t309 * t323 * t324;
	t476 = (t308 * t470 - t322 * t474) / t316 ^ 2;
	t467 = t358 * t497;
	t475 = (t320 * t465 + t346 * t467) / t334 ^ 2;
	t473 = t311 * t353;
	t314 = 0.1e1 / t316;
	t472 = t314 * t324;
	t469 = t330 * t353;
	t468 = t331 * t353;
	t466 = t505 * t358;
	t463 = t388 * t387;
	t462 = t388 * t390;
	t457 = -0.2e1 * t476;
	t456 = -0.2e1 * t475;
	t455 = t312 * t487;
	t454 = t324 * t476;
	t453 = t358 * t475;
	t452 = t308 * t472;
	t451 = t311 * t469;
	t450 = t311 * t468;
	t449 = t328 * t474;
	t448 = t505 * t467;
	t441 = 0.2e1 * t449;
	t327 = -t352 * t390 + t365 * t387;
	t326 = -t352 * t387 - t365 * t390;
	t427 = -t387 * t323 + t390 * t470;
	t426 = t352 * t358 + t362 * t465;
	t420 = -t330 + (-t331 * t466 + t330) * t332;
	t367 = t489 * qJD(3);
	t344 = t433 * qJD(5) + t368 * t388;
	t343 = t366 * t462 - t404 * t387;
	t305 = 0.1e1 / t307;
	t304 = t332 * t496;
	t302 = t426 * t332;
	t300 = t420 * t353;
	t298 = (-t330 * t365 - t331 * t373) * t391 - t435 * t304;
	t297 = t435 * t302 - t330 * t352 + t331 * t362;
	t294 = t426 * t456 + (0.2e1 * t362 * t448 - t319 * t358 + (t320 * t362 + t344 * t505 + t345 * t352) * t359) * t332;
	t293 = 0.2e1 * t475 * t496 + (t425 * t460 + (-0.2e1 * t373 * t448 + t341 * t358 + (-t320 * t373 + t345 * t365 - t367 * t505) * t359) * t391) * t332;
	t1 = [-t453 * t487 + (t317 * t358 + t353 * t497) * t332, 0, t293, 0, t294, 0; t505 * t310 * t458 + (-t320 * t310 + (t296 * t505 - t300 * t317) * t311) * t305 + (t300 * t311 * t458 + (0.2e1 * t300 * t477 - (t301 * t332 * t466 + t456) * t451 - (0.2e1 * t505 * t453 - t301 + (t301 - t428) * t332) * t450 - t420 * t471) * t305) * t353, 0, (t298 * t473 + t310 * t464) * t458 + (-t298 * t471 + (-t338 * t391 + t366 * t460) * t310 + (t298 * t455 + t311 * t464) * t296 - (t373 * t460 - t293 * t505 + t304 * t320 - t367 * t391 + (-t304 * t433 - t365 * t391) * t301) * t450 - (t365 * t460 + t293 * t433 + t304 * t345 - t341 * t391 + (-t304 * t505 + t373 * t391) * t301) * t451) * t305, 0, (t297 * t473 - t310 * t354) * t458 + (t297 * t296 * t455 + t318 * t310 + (-t354 * t296 - t297 * t317 - (-t294 * t505 - t302 * t320 + t344 + (t302 * t433 - t352) * t301) * t468 - (t294 * t433 - t302 * t345 + t319 + (t302 * t505 - t362) * t301) * t469) * t311) * t305, 0; 0.2e1 * (-t323 * t326 + t327 * t470) * t476 + ((t327 * qJD(6) + t319 * t387 - t341 * t390) * t323 + t327 * t441 + (-t326 * t309 - (-t326 * qJD(6) + t319 * t390 + t341 * t387) * t328 - t327 * t308) * t324) * t314, 0, (t454 * t502 - t452) * t343 + (-t309 * t472 + t323 * t457) * (t366 * t463 + t404 * t390) + (t343 * t441 + (t463 * t323 - t462 * t470) * t338 + (t493 * t323 - t491 * t470) * t390 + (t491 * t323 + t493 * t470) * t387) * t314, 0, t427 * t353 * t457 + (t427 * t317 + ((-qJD(6) * t323 - 0.2e1 * t449) * t390 + (t308 * t390 + (t309 - t459) * t387) * t324) * t353) * t314, t457 + (t452 + (-t314 * t474 - t454) * t328) * t502;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end