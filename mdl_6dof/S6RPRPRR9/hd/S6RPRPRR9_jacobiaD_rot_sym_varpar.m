% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR9
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
%   Wie in S6RPRPRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiaD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.36s
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
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:11
	% DurationCPUTime: 0.66s
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
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:11
	% DurationCPUTime: 0.86s
	% Computational Cost: add. (1382->69), mult. (4220->158), div. (139->12), fcn. (5452->15), ass. (0->87)
	t217 = sin(pkin(13));
	t221 = cos(pkin(13));
	t225 = sin(qJ(3));
	t227 = cos(qJ(3));
	t209 = t225 * t217 - t227 * t221;
	t267 = t209 * qJD(3);
	t219 = sin(pkin(7));
	t191 = t209 * t219;
	t223 = cos(pkin(7));
	t193 = t209 * t223;
	t224 = cos(pkin(6));
	t222 = cos(pkin(12));
	t228 = cos(qJ(1));
	t248 = t228 * t222;
	t218 = sin(pkin(12));
	t226 = sin(qJ(1));
	t252 = t226 * t218;
	t237 = t224 * t252 - t248;
	t249 = t228 * t218;
	t251 = t226 * t222;
	t238 = t224 * t251 + t249;
	t240 = t227 * t217 + t225 * t221;
	t220 = sin(pkin(6));
	t255 = t220 * t226;
	t171 = -t191 * t255 + t193 * t238 + t237 * t240;
	t164 = t171 ^ 2;
	t192 = t240 * t219;
	t194 = t240 * t223;
	t235 = t192 * t255 - t194 * t238 + t209 * t237;
	t166 = 0.1e1 / t235 ^ 2;
	t266 = t164 * t166;
	t202 = -t220 * t222 * t219 + t224 * t223;
	t203 = -t224 * t248 + t252;
	t254 = t220 * t228;
	t239 = -t203 * t219 + t223 * t254;
	t179 = atan2(t239, t202);
	t174 = sin(t179);
	t175 = cos(t179);
	t163 = t174 * t239 + t175 * t202;
	t160 = 0.1e1 / t163;
	t165 = 0.1e1 / t235;
	t264 = t239 ^ 2;
	t199 = 0.1e1 / t202;
	t161 = 0.1e1 / t163 ^ 2;
	t200 = 0.1e1 / t202 ^ 2;
	t263 = -0.2e1 * t199 * t200;
	t197 = t238 * qJD(1);
	t243 = qJD(1) * t220 * t223;
	t181 = -t197 * t219 - t226 * t243;
	t178 = t264 * t200 + 0.1e1;
	t176 = 0.1e1 / t178;
	t236 = t174 + (t175 * t199 * t239 - t174) * t176;
	t150 = t236 * t181;
	t262 = t150 * t160 * t161;
	t187 = t267 * t219;
	t189 = t267 * t223;
	t195 = t203 * qJD(1);
	t204 = -t224 * t249 - t251;
	t196 = t204 * qJD(1);
	t208 = t240 * qJD(3);
	t246 = qJD(1) * t228;
	t153 = t238 * t189 + t195 * t194 - t196 * t209 + t237 * t208 + (-t187 * t226 + t192 * t246) * t220;
	t167 = t165 * t166;
	t261 = t167 * t153;
	t169 = t192 * t254 + t203 * t194 - t204 * t209;
	t260 = t169 * t171;
	t259 = t176 * t199;
	t177 = 0.1e1 / t178 ^ 2;
	t258 = t177 * t239;
	t180 = t195 * t219 - t228 * t243;
	t257 = t180 * t161;
	t185 = -t219 * t238 - t223 * t255;
	t256 = t181 * t185;
	t247 = qJD(1) * t226;
	t156 = 0.1e1 + t266;
	t188 = t219 * t208;
	t190 = t223 * t208;
	t152 = t238 * t190 - t195 * t193 - t196 * t240 - t237 * t267 + (-t188 * t226 - t191 * t246) * t220;
	t244 = t171 * t166 * t152;
	t245 = 0.2e1 * (-t164 * t261 + t244) / t156 ^ 2;
	t198 = t237 * qJD(1);
	t182 = t185 ^ 2;
	t168 = t191 * t254 + t203 * t193 + t204 * t240;
	t159 = t182 * t161 + 0.1e1;
	t154 = 0.1e1 / t156;
	t151 = t236 * t185;
	t1 = [t256 * t258 * t263 + t180 * t259, 0, 0, 0, 0, 0; 0.2e1 * (-t151 * t161 * t185 - t160 * t239) / t159 ^ 2 * (-t182 * t262 + t185 * t257) + (t181 * t160 + (-t150 * t239 + t151 * t180) * t161 + (-0.2e1 * t151 * t262 + t236 * t257 + (t174 * t200 * t258 + (0.2e1 * t259 + (t264 * t263 - t199) * t177) * t175) * t161 * t256) * t185) / t159, 0, 0, 0, 0, 0; (-t165 * t168 - t166 * t260) * t245 + ((t203 * t190 + t197 * t193 + t198 * t240 - t204 * t267 + (t188 * t228 - t191 * t247) * t220) * t165 - 0.2e1 * t260 * t261 + (-t168 * t153 + (-t203 * t189 + t197 * t194 - t198 * t209 - t204 * t208 + (-t187 * t228 - t192 * t247) * t220) * t171 + t169 * t152) * t166) * t154, 0, (-t165 * t235 - t266) * t245 + (0.2e1 * t244 + (-0.2e1 * t164 * t167 - t166 * t235 + t165) * t153) * t154, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:13
	% DurationCPUTime: 2.29s
	% Computational Cost: add. (8530->127), mult. (25589->256), div. (448->12), fcn. (32987->17), ass. (0->121)
	t325 = sin(pkin(13));
	t329 = cos(pkin(13));
	t334 = sin(qJ(3));
	t337 = cos(qJ(3));
	t316 = t325 * t334 - t329 * t337;
	t331 = cos(pkin(7));
	t304 = t316 * t331;
	t332 = cos(pkin(6));
	t330 = cos(pkin(12));
	t338 = cos(qJ(1));
	t374 = t338 * t330;
	t326 = sin(pkin(12));
	t335 = sin(qJ(1));
	t378 = t335 * t326;
	t310 = -t332 * t374 + t378;
	t375 = t338 * t326;
	t376 = t335 * t330;
	t311 = t332 * t375 + t376;
	t328 = sin(pkin(6));
	t327 = sin(pkin(7));
	t352 = t316 * t327;
	t350 = t328 * t352;
	t358 = t325 * t337 + t329 * t334;
	t353 = t304 * t310 - t311 * t358 + t338 * t350;
	t274 = t353 ^ 2;
	t290 = (t304 * t330 + t326 * t358) * t328 + t332 * t352;
	t287 = 0.1e1 / t290 ^ 2;
	t264 = t274 * t287 + 0.1e1;
	t262 = 0.1e1 / t264;
	t351 = qJD(3) * t358;
	t300 = t327 * t351;
	t302 = t331 * t351;
	t354 = t332 * t376 + t375;
	t308 = t354 * qJD(1);
	t366 = t332 * t378;
	t373 = qJD(1) * t338;
	t309 = -qJD(1) * t366 + t330 * t373;
	t314 = t316 * qJD(3);
	t346 = qJD(1) * t350;
	t379 = t328 * t338;
	t258 = t300 * t379 + t310 * t302 + t308 * t304 - t309 * t358 + t311 * t314 - t335 * t346;
	t286 = 0.1e1 / t290;
	t284 = t332 * t300 + (t302 * t330 - t314 * t326) * t328;
	t381 = t284 * t287;
	t357 = t258 * t286 - t353 * t381;
	t240 = t357 * t262;
	t265 = atan2(t353, t290);
	t260 = sin(t265);
	t261 = cos(t265);
	t359 = -t260 * t290 + t261 * t353;
	t236 = t240 * t359 + t260 * t258 + t261 * t284;
	t249 = t260 * t353 + t261 * t290;
	t247 = 0.1e1 / t249 ^ 2;
	t394 = t236 * t247;
	t371 = qJD(3) * t337;
	t372 = qJD(3) * t334;
	t393 = t325 * t372 - t329 * t371;
	t246 = 0.1e1 / t249;
	t377 = t335 * t328;
	t295 = t327 * t354 + t331 * t377;
	t333 = sin(qJ(5));
	t336 = cos(qJ(5));
	t303 = t358 * t327;
	t305 = t358 * t331;
	t313 = -t366 + t374;
	t347 = t303 * t377 - t305 * t354 - t313 * t316;
	t273 = t295 * t333 + t336 * t347;
	t267 = 0.1e1 / t273;
	t268 = 0.1e1 / t273 ^ 2;
	t281 = t304 * t354 - t313 * t358 - t335 * t350;
	t275 = t281 ^ 2;
	t245 = t247 * t275 + 0.1e1;
	t306 = t310 * qJD(1);
	t307 = t311 * qJD(1);
	t255 = -t300 * t377 + t302 * t354 - t306 * t304 + t307 * t358 + t313 * t314 - t338 * t346;
	t385 = t255 * t247;
	t391 = t246 * t394;
	t392 = (-t275 * t391 + t281 * t385) / t245 ^ 2;
	t299 = t393 * t327;
	t301 = t393 * t331;
	t315 = -t325 * t371 - t329 * t372;
	t256 = t354 * t301 + t306 * t305 + t307 * t316 + t313 * t315 + (-t299 * t335 + t303 * t373) * t328;
	t365 = qJD(1) * t328 * t331;
	t291 = -t306 * t327 + t338 * t365;
	t250 = qJD(5) * t273 + t256 * t333 - t291 * t336;
	t272 = -t295 * t336 + t333 * t347;
	t266 = t272 ^ 2;
	t254 = t266 * t268 + 0.1e1;
	t384 = t268 * t272;
	t370 = qJD(5) * t272;
	t251 = t256 * t336 + t291 * t333 - t370;
	t386 = t251 * t267 * t268;
	t390 = (t250 * t384 - t266 * t386) / t254 ^ 2;
	t380 = t286 * t381;
	t388 = (t258 * t287 * t353 - t274 * t380) / t264 ^ 2;
	t387 = t247 * t281;
	t383 = t353 * t286;
	t289 = t332 * t303 + (t305 * t330 - t316 * t326) * t328;
	t382 = t353 * t289;
	t368 = 0.2e1 * t390;
	t367 = 0.2e1 * t388;
	t362 = -0.2e1 * t286 * t388;
	t361 = 0.2e1 * t272 * t386;
	t360 = -0.2e1 * t281 * t391;
	t294 = -t310 * t327 + t331 * t379;
	t348 = t303 * t379 + t305 * t310 + t311 * t316;
	t271 = t294 * t333 + t336 * t348;
	t270 = -t294 * t336 + t333 * t348;
	t356 = -t267 * t333 + t336 * t384;
	t355 = -t286 * t348 + t287 * t382;
	t349 = t260 + (t261 * t383 - t260) * t262;
	t345 = -t310 * t301 + t308 * t305 + t309 * t316 - t311 * t315 + (-qJD(1) * t303 * t335 - t299 * t338) * t328;
	t292 = -t308 * t327 - t335 * t365;
	t285 = -t332 * t299 + (-t301 * t330 + t315 * t326) * t328;
	t252 = 0.1e1 / t254;
	t243 = 0.1e1 / t245;
	t241 = t355 * t262;
	t239 = t349 * t281;
	t237 = -t241 * t359 + t260 * t348 + t261 * t289;
	t235 = t355 * t367 + (0.2e1 * t380 * t382 + t345 * t286 + (-t258 * t289 - t284 * t348 - t285 * t353) * t287) * t262;
	t1 = [t281 * t362 + (t255 * t286 - t281 * t381) * t262, 0, t235, 0, 0, 0; -0.2e1 * (t239 * t387 + t246 * t353) * t392 + ((t385 + t360) * t239 + t258 * t246 - t353 * t394 + (t349 * t255 + ((-t240 * t262 * t383 + t367) * t260 + (t353 * t362 + t240 + (-t240 + t357) * t262) * t261) * t281) * t387) * t243, 0, 0.2e1 * (-t237 * t387 - t246 * t347) * t392 + (t256 * t246 + t237 * t360 + (-t347 * t236 + t237 * t255 + ((t235 * t353 - t241 * t258 + t285 + (t241 * t290 + t348) * t240) * t261 + (-t235 * t290 + t241 * t284 + t345 + (t241 * t353 - t289) * t240) * t260) * t281) * t247) * t243, 0, 0, 0; (-t267 * t270 + t271 * t384) * t368 + ((qJD(5) * t271 - t292 * t336 + t333 * t345) * t267 + t271 * t361 + (-t270 * t251 - (-qJD(5) * t270 + t292 * t333 + t336 * t345) * t272 - t271 * t250) * t268) * t252, 0, t356 * t281 * t368 + (-t356 * t255 + ((qJD(5) * t267 + t361) * t336 + (-t250 * t336 + (-t251 + t370) * t333) * t268) * t281) * t252, 0, -0.2e1 * t390 + 0.2e1 * (t250 * t268 * t252 + (-t252 * t386 - t268 * t390) * t272) * t272, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:11
	% EndTime: 2019-10-10 01:00:17
	% DurationCPUTime: 6.05s
	% Computational Cost: add. (20977->187), mult. (61154->360), div. (726->12), fcn. (79375->19), ass. (0->160)
	t435 = sin(qJ(3));
	t527 = sin(pkin(7));
	t528 = cos(pkin(13));
	t476 = t528 * t527;
	t429 = sin(pkin(13));
	t488 = t429 * t527;
	t532 = cos(qJ(3));
	t410 = -t435 * t488 + t532 * t476;
	t406 = t410 * qJD(3);
	t529 = cos(pkin(7));
	t477 = t529 * t528;
	t489 = t429 * t529;
	t412 = -t435 * t489 + t532 * t477;
	t408 = t412 * qJD(3);
	t432 = cos(pkin(6));
	t437 = cos(qJ(1));
	t526 = sin(pkin(12));
	t486 = t437 * t526;
	t431 = cos(pkin(12));
	t530 = sin(qJ(1));
	t492 = t530 * t431;
	t456 = t432 * t492 + t486;
	t415 = t456 * qJD(1);
	t427 = t530 * t526;
	t473 = t432 * t427;
	t507 = qJD(1) * t437;
	t416 = -qJD(1) * t473 + t431 * t507;
	t423 = -t532 * t429 - t435 * t528;
	t421 = t423 * qJD(3);
	t422 = t435 * t429 - t532 * t528;
	t430 = sin(pkin(6));
	t450 = -t435 * t476 - t532 * t488;
	t451 = -t435 * t477 - t532 * t489;
	t457 = -t432 * t486 - t492;
	t508 = t437 * t431;
	t458 = -t432 * t508 + t427;
	t368 = t458 * t408 - t415 * t451 + t416 * t422 + t457 * t421 - (-qJD(1) * t450 * t530 - t437 * t406) * t430;
	t509 = t430 * t437;
	t391 = -t422 * t457 - t450 * t509 - t451 * t458;
	t487 = t430 * t529;
	t404 = t437 * t487 - t458 * t527;
	t434 = sin(qJ(5));
	t531 = cos(qJ(5));
	t379 = t391 * t531 + t404 * t434;
	t474 = t530 * t487;
	t453 = qJD(1) * t474 + t415 * t527;
	t548 = t379 * qJD(5) + t368 * t434 + t453 * t531;
	t544 = -t391 * t434 + t404 * t531;
	t347 = qJD(5) * t544 + t368 * t531 - t453 * t434;
	t373 = t544 ^ 2;
	t398 = -t432 * t450 + (-t526 * t422 - t431 * t451) * t430;
	t417 = -t430 * t431 * t527 + t432 * t529;
	t464 = -t398 * t434 + t417 * t531;
	t384 = 0.1e1 / t464 ^ 2;
	t360 = t373 * t384 + 0.1e1;
	t358 = 0.1e1 / t360;
	t387 = t398 * t531 + t417 * t434;
	t396 = t432 * t406 + (t408 * t431 + t526 * t421) * t430;
	t371 = qJD(5) * t387 + t396 * t434;
	t383 = 0.1e1 / t464;
	t512 = t544 * t384;
	t471 = t371 * t512 - t383 * t548;
	t327 = t471 * t358;
	t361 = atan2(-t544, -t464);
	t356 = sin(t361);
	t357 = cos(t361);
	t475 = t356 * t464 - t357 * t544;
	t322 = t475 * t327 + t356 * t548 + t357 * t371;
	t339 = -t356 * t544 - t357 * t464;
	t337 = 0.1e1 / t339 ^ 2;
	t547 = t322 * t337;
	t447 = t456 * t527 + t474;
	t419 = -t473 + t508;
	t493 = t430 * t530;
	t448 = -t419 * t422 - t450 * t493 + t451 * t456;
	t381 = t447 * t434 + t448 * t531;
	t482 = t410 * t493;
	t393 = -t456 * t412 + t419 * t423 + t482;
	t433 = sin(qJ(6));
	t436 = cos(qJ(6));
	t354 = t381 * t433 + t393 * t436;
	t543 = 0.2e1 * t354;
	t336 = 0.1e1 / t339;
	t542 = t336 * t547;
	t445 = t447 * t531;
	t380 = t434 * t448 - t445;
	t533 = 0.2e1 * t380;
	t483 = t533 * t542;
	t414 = t457 * qJD(1);
	t454 = t458 * qJD(1);
	t491 = t430 * t507;
	t444 = t406 * t493 - t456 * t408 - t414 * t422 + t419 * t421 - t450 * t491 - t451 * t454;
	t446 = t404 * qJD(1);
	t343 = t381 * qJD(5) + t434 * t444 - t531 * t446;
	t518 = t343 * t337;
	t539 = -t518 + t483;
	t374 = t380 ^ 2;
	t335 = t374 * t337 + 0.1e1;
	t504 = 0.2e1 * (-t374 * t542 + t380 * t518) / t335 ^ 2;
	t538 = t371 * t384;
	t397 = t432 * t410 + (t412 * t431 + t526 * t423) * t430;
	t535 = -t410 * t509 - t412 * t458 - t423 * t457;
	t468 = t383 * t535 + t397 * t512;
	t537 = t434 * t468;
	t490 = qJD(5) * t531;
	t407 = t450 * qJD(3);
	t409 = t451 * qJD(3);
	t420 = t422 * qJD(3);
	t536 = qJD(1) * t482 - t407 * t509 - t409 * t458 - t415 * t412 + t416 * t423 - t420 * t457;
	t355 = t381 * t436 - t393 * t433;
	t349 = 0.1e1 / t355;
	t350 = 0.1e1 / t355 ^ 2;
	t534 = -0.2e1 * t544;
	t506 = qJD(5) * t434;
	t344 = qJD(5) * t445 + t434 * t446 + t444 * t531 - t448 * t506;
	t363 = t407 * t493 - t456 * t409 + t410 * t491 + t412 * t454 + t414 * t423 + t419 * t420;
	t331 = qJD(6) * t355 + t344 * t433 + t363 * t436;
	t348 = t354 ^ 2;
	t342 = t348 * t350 + 0.1e1;
	t517 = t350 * t354;
	t505 = qJD(6) * t354;
	t332 = t344 * t436 - t363 * t433 - t505;
	t521 = t332 * t349 * t350;
	t524 = (t331 * t517 - t348 * t521) / t342 ^ 2;
	t514 = t383 * t538;
	t522 = (t373 * t514 - t512 * t548) / t360 ^ 2;
	t520 = t337 * t380;
	t340 = 0.1e1 / t342;
	t519 = t340 * t350;
	t516 = t356 * t380;
	t515 = t357 * t380;
	t513 = t544 * t383;
	t510 = t393 * t434;
	t503 = -0.2e1 * t524;
	t502 = -0.2e1 * t522;
	t500 = t350 * t524;
	t499 = t383 * t522;
	t498 = t331 * t519;
	t497 = t354 * t521;
	t494 = t393 * t531;
	t485 = 0.2e1 * t497;
	t484 = t514 * t534;
	t353 = t379 * t436 + t433 * t535;
	t352 = t379 * t433 - t436 * t535;
	t472 = qJD(6) * t494 - t444;
	t470 = -t433 * t349 + t436 * t517;
	t469 = -t379 * t383 + t387 * t512;
	t460 = -t356 + (-t357 * t513 + t356) * t358;
	t455 = qJD(6) * t448 + t531 * t363 - t393 * t506;
	t395 = t432 * t407 + (t409 * t431 + t526 * t420) * t430;
	t372 = t464 * qJD(5) + t396 * t531;
	t370 = t433 * t448 + t436 * t494;
	t333 = 0.1e1 / t335;
	t330 = t358 * t537;
	t329 = t469 * t358;
	t324 = (-t356 * t535 + t357 * t397) * t434 + t475 * t330;
	t323 = t475 * t329 + t356 * t379 + t357 * t387;
	t320 = t469 * t502 + (-t387 * t484 - t347 * t383 + (-t371 * t379 + t372 * t544 - t387 * t548) * t384) * t358;
	t319 = t502 * t537 + ((-t397 * t484 + t536 * t383 + (t371 * t535 + t395 * t544 - t397 * t548) * t384) * t434 + t468 * t490) * t358;
	t1 = [-t499 * t533 + (t343 * t383 + t380 * t538) * t358, 0, t319, 0, t320, 0; t544 * t336 * t504 + (t548 * t336 + t544 * t547 - (t460 * t343 + ((t327 * t358 * t513 + t502) * t356 + (-t499 * t534 - t327 + (t327 - t471) * t358) * t357) * t380) * t520) * t333 + (t539 * t333 + t520 * t504) * t460 * t380, 0, (t324 * t520 - t336 * t510) * t504 + ((t363 * t434 + t393 * t490) * t336 + t539 * t324 + (-t510 * t322 - (t397 * t490 - t319 * t544 + t330 * t548 + t395 * t434 + (t330 * t464 - t434 * t535) * t327) * t515 - (-t535 * t490 + t319 * t464 - t330 * t371 - t536 * t434 + (t330 * t544 - t397 * t434) * t327) * t516) * t337) * t333, 0, (t323 * t520 - t336 * t381) * t504 + (t323 * t483 + t344 * t336 + (-t381 * t322 - t323 * t343 - (-t320 * t544 + t329 * t548 + t372 + (t329 * t464 + t379) * t327) * t515 - (t320 * t464 - t329 * t371 + t347 + (t329 * t544 - t387) * t327) * t516) * t337) * t333, 0; 0.2e1 * (-t349 * t352 + t353 * t517) * t524 + ((qJD(6) * t353 + t347 * t433 - t436 * t536) * t349 + t353 * t485 + (-t352 * t332 - (-qJD(6) * t352 + t347 * t436 + t433 * t536) * t354 - t353 * t331) * t350) * t340, 0, (t500 * t543 - t498) * t370 + (-t332 * t519 + t349 * t503) * (t433 * t494 - t436 * t448) + ((t433 * t455 + t436 * t472) * t349 - (-t433 * t472 + t436 * t455) * t517 + t370 * t485) * t340, 0, t470 * t380 * t503 + (t470 * t343 + ((-qJD(6) * t349 - 0.2e1 * t497) * t436 + (t331 * t436 + (t332 - t505) * t433) * t350) * t380) * t340, t503 + (t498 + (-t340 * t521 - t500) * t354) * t543;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end