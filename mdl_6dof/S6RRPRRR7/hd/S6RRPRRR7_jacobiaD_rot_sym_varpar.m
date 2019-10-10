% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR7
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
%   Wie in S6RRPRRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:59
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:39
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (776->71), mult. (1835->159), div. (470->13), fcn. (2177->7), ass. (0->71)
	t82 = sin(qJ(1));
	t113 = qJD(1) * t82;
	t133 = 0.2e1 * t82;
	t73 = t82 ^ 2;
	t84 = cos(qJ(1));
	t77 = t84 ^ 2;
	t78 = 0.1e1 / t84;
	t131 = (t73 / t77 + 0.1e1) * t78 * t113;
	t81 = sin(qJ(2));
	t114 = t82 * t81;
	t83 = cos(qJ(2));
	t63 = atan2(-t114, -t83);
	t61 = sin(t63);
	t105 = t61 * t114;
	t62 = cos(t63);
	t57 = -t62 * t83 - t105;
	t54 = 0.1e1 / t57;
	t74 = 0.1e1 / t83;
	t55 = 0.1e1 / t57 ^ 2;
	t75 = 0.1e1 / t83 ^ 2;
	t130 = -0.2e1 * t81;
	t71 = t81 ^ 2;
	t118 = t71 * t75;
	t68 = t73 * t118 + 0.1e1;
	t64 = 0.1e1 / t68;
	t129 = t64 - 0.1e1;
	t112 = qJD(1) * t84;
	t103 = t81 * t112;
	t111 = qJD(2) * t82;
	t120 = t62 * t81;
	t110 = qJD(2) * t83;
	t50 = (-(-t82 * t110 - t103) * t74 + t111 * t118) * t64;
	t46 = (-t50 * t82 + qJD(2)) * t120 + (-t103 + (t50 - t111) * t83) * t61;
	t128 = t46 * t54 * t55;
	t127 = t50 * t61;
	t126 = t50 * t81;
	t125 = t55 * t81;
	t124 = t55 * t84;
	t115 = t74 * t81;
	t70 = t81 * t71;
	t76 = t74 * t75;
	t92 = qJD(2) * (t70 * t76 + t115);
	t96 = t71 * t82 * t112;
	t123 = (t73 * t92 + t75 * t96) / t68 ^ 2;
	t102 = 0.1e1 + t118;
	t60 = t102 * t82 * t64;
	t122 = t60 * t82;
	t121 = t61 * t83;
	t119 = t71 * t74;
	t117 = t71 * t77;
	t116 = t73 / t84 ^ 2;
	t53 = t55 * t117 + 0.1e1;
	t109 = 0.2e1 * (-t117 * t128 + (t77 * t81 * t110 - t96) * t55) / t53 ^ 2;
	t108 = 0.2e1 * t128;
	t69 = t75 * t116 + 0.1e1;
	t107 = 0.2e1 * (t76 * qJD(2) * t81 * t116 + t75 * t131) / t69 ^ 2;
	t106 = t81 * t124;
	t104 = t64 * t119;
	t101 = 0.1e1 + t116;
	t100 = t81 * t109;
	t99 = t123 * t130;
	t98 = t123 * t133;
	t97 = t82 * t104;
	t95 = t102 * t84;
	t93 = t101 * t81 * t75;
	t66 = 0.1e1 / t69;
	t51 = 0.1e1 / t53;
	t49 = (t129 * t81 * t61 - t62 * t97) * t84;
	t48 = -t82 * t121 + t120 + (-t62 * t114 + t121) * t60;
	t47 = -t102 * t98 + (qJD(1) * t95 + t92 * t133) * t64;
	t1 = [t84 * t74 * t99 + (qJD(2) * t95 - t113 * t115) * t64, t47, 0, 0, 0, 0; (t54 * t100 + (-t54 * t110 + (qJD(1) * t49 + t46) * t125) * t51) * t82 + (t55 * t100 * t49 + (-((t129 * t110 + t50 * t97 + t99) * t61 + (t98 * t119 - t126 + (t126 + (-t70 * t75 + t130) * t111) * t64) * t62) * t106 + (t81 * t108 - t55 * t110) * t49 + (-t54 + ((-t73 + t77) * t62 * t104 + t129 * t105) * t55) * t81 * qJD(1)) * t51) * t84, (t48 * t125 - t54 * t83) * t84 * t109 + ((-t54 * t113 + (-qJD(2) * t48 - t46) * t124) * t83 + (-t84 * qJD(2) * t54 - (-t47 * t62 * t82 + t61 * t111 + t122 * t127 - t127 + (-qJD(2) * t61 - t112 * t62) * t60) * t106 + (t84 * t108 + t55 * t113) * t48 - ((t47 - t112) * t61 + ((0.1e1 - t122) * qJD(2) + (t60 - t82) * t50) * t62) * t83 * t124) * t81) * t51, 0, 0, 0, 0; t101 * t74 * t107 + (-qJD(2) * t93 - 0.2e1 * t74 * t131) * t66, t78 * t75 * t107 * t114 + ((-0.2e1 * t71 * t76 - t74) * t78 * t111 - qJD(1) * t93) * t66, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:40
	% DurationCPUTime: 1.88s
	% Computational Cost: add. (3233->93), mult. (9702->187), div. (711->12), fcn. (12089->11), ass. (0->94)
	t280 = sin(qJ(4));
	t281 = sin(qJ(2));
	t283 = cos(qJ(4));
	t284 = cos(qJ(2));
	t201 = -t284 * t280 + t281 * t283;
	t282 = sin(qJ(1));
	t192 = t201 * t282;
	t200 = t281 * t280 + t284 * t283;
	t179 = atan2(t192, t200);
	t174 = sin(t179);
	t175 = cos(t179);
	t194 = t200 * t282;
	t244 = -t174 * t200 + t175 * t192;
	t190 = t192 ^ 2;
	t198 = 0.1e1 / t200 ^ 2;
	t178 = t190 * t198 + 0.1e1;
	t176 = 0.1e1 / t178;
	t197 = 0.1e1 / t200;
	t263 = t192 * t201;
	t240 = -t194 * t197 - t198 * t263;
	t292 = t240 * t176;
	t151 = -t174 * t194 + t175 * t201 + t244 * t292;
	t285 = cos(qJ(1));
	t196 = t201 * t285;
	t311 = t151 * t196;
	t166 = t174 * t192 + t175 * t200;
	t163 = 0.1e1 / t166;
	t164 = 0.1e1 / t166 ^ 2;
	t195 = t200 * t285;
	t191 = t196 ^ 2;
	t160 = t191 * t164 + 0.1e1;
	t300 = qJD(2) - qJD(4);
	t168 = -t192 * qJD(1) + t300 * t195;
	t271 = t168 * t164;
	t169 = -t196 * qJD(1) - t300 * t194;
	t181 = t300 * t201;
	t264 = t192 * t198;
	t242 = -t169 * t197 + t181 * t264;
	t154 = t242 * t176;
	t149 = t244 * t154 - t174 * t169 - t175 * t181;
	t279 = t149 * t163 * t164;
	t261 = 0.2e1 * (-t191 * t279 + t196 * t271) / t160 ^ 2;
	t310 = (t163 * t195 + t164 * t311) * t261;
	t167 = t194 * qJD(1) + t300 * t196;
	t219 = sin(qJ(5));
	t220 = cos(qJ(5));
	t189 = t195 * t220 - t282 * t219;
	t254 = qJD(1) * t285;
	t161 = t189 * qJD(5) - t167 * t219 + t220 * t254;
	t238 = -t195 * t219 - t282 * t220;
	t301 = t238 * qJD(5);
	t162 = -t167 * t220 - t219 * t254 + t301;
	t182 = t238 ^ 2;
	t184 = 0.1e1 / t189 ^ 2;
	t173 = t182 * t184 + 0.1e1;
	t171 = 0.1e1 / t173;
	t183 = 0.1e1 / t189;
	t266 = t184 * t238;
	t241 = -t219 * t183 - t220 * t266;
	t275 = t162 * t183 * t184;
	t287 = -0.2e1 * t238;
	t251 = t275 * t287;
	t309 = (t241 * t168 - (((-t162 - t301) * t219 - t161 * t220) * t184 + (qJD(5) * t183 + t251) * t220) * t196) * t171;
	t308 = -t195 * t149 + t151 * t168;
	t260 = 0.2e1 * t279;
	t307 = -t167 * t163 - t260 * t311;
	t180 = t300 * t200;
	t305 = (t200 * t292 + t194) * t154 + t292 * t169 - t180;
	t170 = t195 * qJD(1) - t300 * t192;
	t304 = -(-t192 * t292 - t201) * t154 - t292 * t181 + t170;
	t294 = t181 * t198;
	t267 = t197 * t294;
	t302 = ((t169 * t201 - t180 * t192 - t181 * t194) * t198 - t170 * t197 - 0.2e1 * t263 * t267) * t176;
	t265 = t192 * t197;
	t290 = (-t175 * t265 + t174) * t176 - t174;
	t286 = -0.2e1 * t196;
	t278 = (-t161 * t266 - t182 * t275) / t173 ^ 2;
	t277 = (-t169 * t264 + t190 * t267) / t178 ^ 2;
	t270 = t171 * t184;
	t269 = t174 * t196;
	t268 = t175 * t196;
	t259 = -0.2e1 * t278;
	t258 = -0.2e1 * t277;
	t257 = t184 * t278;
	t256 = t197 * t277;
	t255 = t161 * t270;
	t253 = qJD(1) * t282;
	t187 = -t194 * t220 - t285 * t219;
	t239 = t194 * t219 - t285 * t220;
	t158 = 0.1e1 / t160;
	t153 = t290 * t196;
	t147 = t240 * t258 + t302;
	t146 = 0.2e1 * t240 * t277 - t302;
	t1 = [t256 * t286 + (t168 * t197 + t196 * t294) * t176, t146, 0, t147, 0, 0; -t192 * t163 * t261 + (-t169 * t163 + (-t149 * t192 - t153 * t168) * t164) * t158 - ((-t153 * t260 + t290 * t271) * t158 + (-t153 * t261 + ((t154 * t176 * t265 + t258) * t269 + (0.2e1 * t192 * t256 - t154 + (t154 - t242) * t176) * t268) * t158) * t164) * t196, t310 + (((t146 * t192 + t305) * t268 + (-t146 * t200 + t304) * t269 - t308) * t164 - t307) * t158, 0, -t310 + (((t147 * t192 - t305) * t268 + (-t147 * t200 - t304) * t269 + t308) * t164 + t307) * t158, 0, 0; (t257 * t287 - t255) * t187 - (-t162 * t270 + t183 * t259) * t239 + ((t187 * qJD(5) - t170 * t219 - t220 * t253) * t183 + (t239 * qJD(5) - t170 * t220 + t219 * t253) * t266 + t187 * t251) * t171, t241 * t278 * t286 + t309, 0, -t241 * t196 * t259 - t309, t259 + (t255 - (-t171 * t275 - t257) * t238) * t287, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:59:38
	% EndTime: 2019-10-10 10:59:40
	% DurationCPUTime: 1.72s
	% Computational Cost: add. (3949->95), mult. (10113->188), div. (729->12), fcn. (12551->11), ass. (0->97)
	t309 = sin(qJ(4));
	t310 = sin(qJ(2));
	t312 = cos(qJ(4));
	t313 = cos(qJ(2));
	t228 = -t313 * t309 + t310 * t312;
	t311 = sin(qJ(1));
	t219 = t228 * t311;
	t227 = t310 * t309 + t313 * t312;
	t210 = atan2(t219, t227);
	t201 = sin(t210);
	t202 = cos(t210);
	t221 = t227 * t311;
	t273 = -t201 * t227 + t202 * t219;
	t217 = t219 ^ 2;
	t225 = 0.1e1 / t227 ^ 2;
	t205 = t217 * t225 + 0.1e1;
	t203 = 0.1e1 / t205;
	t224 = 0.1e1 / t227;
	t292 = t219 * t228;
	t268 = -t221 * t224 - t225 * t292;
	t322 = t268 * t203;
	t178 = -t201 * t221 + t202 * t228 + t273 * t322;
	t314 = cos(qJ(1));
	t223 = t228 * t314;
	t339 = t178 * t223;
	t193 = t201 * t219 + t202 * t227;
	t190 = 0.1e1 / t193;
	t191 = 0.1e1 / t193 ^ 2;
	t222 = t227 * t314;
	t218 = t223 ^ 2;
	t187 = t191 * t218 + 0.1e1;
	t329 = qJD(2) - qJD(4);
	t198 = -t219 * qJD(1) + t329 * t222;
	t299 = t198 * t191;
	t199 = -t223 * qJD(1) - t329 * t221;
	t216 = t329 * t228;
	t293 = t219 * t225;
	t270 = -t199 * t224 + t216 * t293;
	t181 = t270 * t203;
	t176 = t273 * t181 - t199 * t201 - t202 * t216;
	t308 = t176 * t190 * t191;
	t290 = 0.2e1 * (-t218 * t308 + t223 * t299) / t187 ^ 2;
	t338 = (t190 * t222 + t191 * t339) * t290;
	t249 = qJ(5) + qJ(6);
	t246 = sin(t249);
	t247 = cos(t249);
	t248 = qJD(5) + qJD(6);
	t259 = qJD(1) * t314 + t222 * t248;
	t197 = t221 * qJD(1) + t329 * t223;
	t274 = -t311 * t248 - t197;
	t188 = t274 * t246 + t259 * t247;
	t189 = -t259 * t246 + t274 * t247;
	t213 = t222 * t246 + t311 * t247;
	t206 = t213 ^ 2;
	t214 = t222 * t247 - t311 * t246;
	t208 = 0.1e1 / t214 ^ 2;
	t196 = t206 * t208 + 0.1e1;
	t194 = 0.1e1 / t196;
	t207 = 0.1e1 / t214;
	t296 = t208 * t213;
	t269 = -t246 * t207 + t247 * t296;
	t304 = t189 * t207 * t208;
	t316 = 0.2e1 * t213;
	t283 = t304 * t316;
	t337 = (t269 * t198 - ((-t188 * t247 + (t213 * t248 - t189) * t246) * t208 + (t207 * t248 + t283) * t247) * t223) * t194;
	t336 = -t222 * t176 + t178 * t198;
	t289 = 0.2e1 * t308;
	t335 = -t197 * t190 - t289 * t339;
	t215 = t329 * t227;
	t333 = (t227 * t322 + t221) * t181 + t322 * t199 - t215;
	t200 = t222 * qJD(1) - t329 * t219;
	t332 = -(-t219 * t322 - t228) * t181 - t322 * t216 + t200;
	t323 = t216 * t225;
	t295 = t224 * t323;
	t330 = ((t199 * t228 - t215 * t219 - t216 * t221) * t225 - t200 * t224 - 0.2e1 * t292 * t295) * t203;
	t294 = t219 * t224;
	t319 = (-t202 * t294 + t201) * t203 - t201;
	t315 = -0.2e1 * t223;
	t307 = (t188 * t296 - t206 * t304) / t196 ^ 2;
	t306 = (-t199 * t293 + t217 * t295) / t205 ^ 2;
	t301 = t194 * t208;
	t298 = t201 * t223;
	t297 = t202 * t223;
	t288 = -0.2e1 * t307;
	t287 = -0.2e1 * t306;
	t286 = t208 * t307;
	t285 = t224 * t306;
	t284 = t188 * t301;
	t275 = -t314 * t248 - t200;
	t258 = qJD(1) * t311 + t221 * t248;
	t212 = -t221 * t247 - t314 * t246;
	t185 = 0.1e1 / t187;
	t180 = t319 * t223;
	t174 = t268 * t287 + t330;
	t173 = 0.2e1 * t268 * t306 - t330;
	t172 = t288 + (t284 + (-t194 * t304 - t286) * t213) * t316;
	t1 = [t285 * t315 + (t198 * t224 + t223 * t323) * t203, t173, 0, t174, 0, 0; -t219 * t190 * t290 + (-t199 * t190 + (-t176 * t219 - t180 * t198) * t191) * t185 - ((-t180 * t289 + t319 * t299) * t185 + (-t180 * t290 + ((t181 * t203 * t294 + t287) * t298 + (0.2e1 * t219 * t285 - t181 + (t181 - t270) * t203) * t297) * t185) * t191) * t223, t338 + (((t173 * t219 + t333) * t297 + (-t173 * t227 + t332) * t298 - t336) * t191 - t335) * t185, 0, -t338 + (((t174 * t219 - t333) * t297 + (-t174 * t227 - t332) * t298 + t336) * t191 + t335) * t185, 0, 0; (t286 * t316 - t284) * t212 + (-t189 * t301 + t207 * t288) * (-t221 * t246 + t314 * t247) + ((t275 * t246 - t258 * t247) * t207 - (t258 * t246 + t275 * t247) * t296 + t212 * t283) * t194, t269 * t307 * t315 + t337, 0, -t269 * t223 * t288 - t337, t172, t172;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end