% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR5
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
%   Wie in S6RRPRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(6));
	t93 = t99 ^ 2;
	t100 = cos(pkin(6));
	t95 = 0.1e1 / t100 ^ 2;
	t104 = cos(qJ(1));
	t98 = t104 ^ 2;
	t89 = t98 * t93 * t95 + 0.1e1;
	t102 = sin(qJ(1));
	t97 = t102 ^ 2;
	t126 = 0.1e1 / t89 ^ 2 * t97;
	t131 = t126 * t95;
	t122 = t104 * t99;
	t88 = atan2(t122, t100);
	t84 = sin(t88);
	t85 = cos(t88);
	t72 = t85 * t100 + t84 * t122;
	t67 = 0.1e1 / t72;
	t103 = cos(qJ(2));
	t118 = t104 * t103;
	t101 = sin(qJ(2));
	t121 = t102 * t101;
	t113 = t100 * t121 - t118;
	t77 = 0.1e1 / t113;
	t94 = 0.1e1 / t100;
	t68 = 0.1e1 / t72 ^ 2;
	t78 = 0.1e1 / t113 ^ 2;
	t119 = t104 * t101;
	t120 = t102 * t103;
	t81 = -t100 * t119 - t120;
	t82 = t100 * t120 + t119;
	t71 = t81 * qJD(1) - t82 * qJD(2);
	t128 = t71 * t77 * t78;
	t115 = t100 * t118;
	t70 = -qJD(1) * t115 - qJD(2) * t118 + (qJD(2) * t100 + qJD(1)) * t121;
	t129 = t70 * t78;
	t76 = t82 ^ 2;
	t75 = t76 * t78 + 0.1e1;
	t130 = (t76 * t128 - t82 * t129) / t75 ^ 2;
	t127 = t81 * t82;
	t125 = t93 * t94;
	t124 = t102 * t68;
	t123 = t104 * t68;
	t117 = qJD(1) * t104;
	t86 = 0.1e1 / t89;
	t116 = (t86 - 0.1e1) * t99;
	t114 = -0.2e1 * t94 * t131;
	t80 = t115 - t121;
	t63 = (-t104 * t85 * t86 * t125 + t84 * t116) * t102;
	t92 = t99 * t93;
	t73 = 0.1e1 / t75;
	t69 = t67 * t68;
	t66 = t97 * t93 * t68 + 0.1e1;
	t62 = qJD(1) * t63;
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (454->46), mult. (1493->122), div. (133->12), fcn. (1882->11), ass. (0->62)
	t133 = sin(qJ(1));
	t134 = cos(qJ(1));
	t131 = cos(pkin(6));
	t129 = sin(pkin(11));
	t132 = sin(qJ(2));
	t161 = cos(pkin(11));
	t166 = cos(qJ(2));
	t145 = -t132 * t129 + t166 * t161;
	t144 = t145 * t131;
	t146 = t129 * t166 + t132 * t161;
	t99 = -t133 * t144 - t134 * t146;
	t92 = t99 ^ 2;
	t109 = t146 * t131;
	t147 = t133 * t109 - t134 * t145;
	t94 = 0.1e1 / t147 ^ 2;
	t169 = t92 * t94;
	t125 = 0.1e1 / t131 ^ 2;
	t130 = sin(pkin(6));
	t123 = t130 ^ 2;
	t128 = t134 ^ 2;
	t119 = t128 * t123 * t125 + 0.1e1;
	t127 = t133 ^ 2;
	t159 = 0.1e1 / t119 ^ 2 * t127;
	t168 = t125 * t159;
	t167 = qJD(1) * t144 + t145 * qJD(2);
	t93 = 0.1e1 / t147;
	t155 = t134 * t130;
	t118 = atan2(t155, t131);
	t114 = sin(t118);
	t115 = cos(t118);
	t105 = t114 * t155 + t115 * t131;
	t102 = 0.1e1 / t105;
	t124 = 0.1e1 / t131;
	t103 = 0.1e1 / t105 ^ 2;
	t163 = t94 * t99;
	t111 = t146 * qJD(2);
	t108 = t131 * t111;
	t156 = t133 * t146;
	t83 = -qJD(1) * t156 - t133 * t108 + t167 * t134;
	t151 = t83 * t163;
	t107 = qJD(2) * t144;
	t97 = -t134 * t109 - t133 * t145;
	t84 = qJD(1) * t97 - t133 * t107 - t134 * t111;
	t95 = t93 * t94;
	t162 = t95 * t84;
	t88 = 0.1e1 + t169;
	t165 = (t162 * t92 - t151) / t88 ^ 2;
	t160 = t103 * t133;
	t158 = t123 * t124;
	t154 = qJD(1) * t134;
	t152 = -0.2e1 * t97 * t99;
	t116 = 0.1e1 / t119;
	t150 = (t116 - 0.1e1) * t130;
	t149 = -0.2e1 * t124 * t168;
	t85 = (-t115 * t116 * t134 * t158 + t114 * t150) * t133;
	t122 = t130 * t123;
	t104 = t102 * t103;
	t96 = t134 * t144 - t156;
	t91 = t127 * t123 * t103 + 0.1e1;
	t86 = 0.1e1 / t88;
	t82 = qJD(1) * t85;
	t1 = [(-t116 * t124 * t130 + t122 * t149) * t154, 0, 0, 0, 0, 0; (0.2e1 * (-t102 * t134 + t160 * t85) / t91 ^ 2 * (-t104 * t127 * t82 + t154 * t160) * t123 + ((0.2e1 * t104 * t133 * t85 - t103 * t134) * t82 + (-t133 * t102 + ((-t85 + (-t122 * t168 - t150) * t133 * t114) * t134 - (t128 * t123 ^ 2 * t149 + (-t159 + (0.2e1 * t127 - t128) * t116) * t158) * t133 * t115) * t103) * qJD(1)) / t91) * t130, 0, 0, 0, 0, 0; (t94 * t152 + 0.2e1 * t96 * t93) * t165 + ((-t97 * t83 - t96 * t84) * t94 - (-t134 * t108 - t167 * t133 - t146 * t154) * t93 + (qJD(1) * t147 - t134 * t107 + t133 * t111) * t163 - t152 * t162) * t86, 0.2e1 * (-t147 * t93 - t169) * t165 + (-0.2e1 * t151 + (t147 * t94 + 0.2e1 * t92 * t95 - t93) * t84) * t86, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:39
	% DurationCPUTime: 1.27s
	% Computational Cost: add. (3306->95), mult. (9869->202), div. (448->12), fcn. (12731->13), ass. (0->96)
	t216 = sin(qJ(4));
	t218 = cos(qJ(4));
	t271 = sin(pkin(11));
	t273 = cos(pkin(6));
	t241 = t273 * t271;
	t272 = cos(pkin(11));
	t242 = t273 * t272;
	t274 = sin(qJ(2));
	t275 = cos(qJ(2));
	t206 = t241 * t275 + t242 * t274;
	t208 = t271 * t274 - t272 * t275;
	t217 = sin(qJ(1));
	t219 = cos(qJ(1));
	t238 = t217 * t206 + t219 * t208;
	t215 = sin(pkin(6));
	t256 = t215 * t217;
	t233 = t216 * t238 + t218 * t256;
	t278 = qJD(4) * t233;
	t230 = t275 * t271 + t274 * t272;
	t205 = t230 * t215;
	t196 = qJD(2) * t205;
	t204 = t208 * t215;
	t201 = 0.1e1 / t204 ^ 2;
	t258 = t196 * t201;
	t277 = t274 * t241 - t275 * t242;
	t229 = t208 * qJD(2);
	t188 = -t217 * t230 - t219 * t277;
	t176 = atan2(t188, t204);
	t171 = sin(t176);
	t172 = cos(t176);
	t185 = t188 ^ 2;
	t175 = t185 * t201 + 0.1e1;
	t173 = 0.1e1 / t175;
	t200 = 0.1e1 / t204;
	t260 = t188 * t200;
	t276 = t173 * (t172 * t260 - t171) + t171;
	t160 = t171 * t188 + t172 * t204;
	t157 = 0.1e1 / t160;
	t184 = t216 * t256 - t218 * t238;
	t178 = 0.1e1 / t184;
	t158 = 0.1e1 / t160 ^ 2;
	t179 = 0.1e1 / t184 ^ 2;
	t226 = t217 * t277;
	t191 = -t219 * t230 + t226;
	t186 = t191 ^ 2;
	t156 = t186 * t158 + 0.1e1;
	t199 = t206 * qJD(2);
	t166 = qJD(1) * t188 - t217 * t199 - t219 * t229;
	t264 = t166 * t158;
	t254 = qJD(1) * t219;
	t169 = qJD(1) * t226 - t219 * t199 + t217 * t229 - t230 * t254;
	t237 = t169 * t200 - t188 * t258;
	t151 = t237 * t173;
	t240 = -t171 * t204 + t172 * t188;
	t147 = t151 * t240 + t171 * t169 + t172 * t196;
	t269 = t147 * t157 * t158;
	t270 = (-t186 * t269 - t191 * t264) / t156 ^ 2;
	t239 = -t219 * t206 + t217 * t208;
	t259 = t188 * t205;
	t235 = -t200 * t239 + t201 * t259;
	t152 = t235 * t173;
	t148 = -t152 * t240 + t171 * t239 + t172 * t205;
	t268 = t148 * t191;
	t198 = t277 * qJD(2);
	t207 = t230 * qJD(2);
	t167 = qJD(1) * t239 + t217 * t198 - t219 * t207;
	t248 = t215 * t254;
	t161 = qJD(4) * t184 + t167 * t216 - t218 * t248;
	t177 = t233 ^ 2;
	t165 = t177 * t179 + 0.1e1;
	t261 = t179 * t233;
	t162 = t167 * t218 + t216 * t248 + t278;
	t265 = t162 * t178 * t179;
	t267 = (-t161 * t261 - t177 * t265) / t165 ^ 2;
	t257 = t200 * t258;
	t266 = (t188 * t201 * t169 - t185 * t257) / t175 ^ 2;
	t263 = t171 * t191;
	t262 = t172 * t191;
	t255 = t215 * t219;
	t253 = -0.2e1 * t270;
	t252 = -0.2e1 * t269;
	t251 = 0.2e1 * t267;
	t250 = 0.2e1 * t266;
	t249 = qJD(1) * t256;
	t247 = -0.2e1 * t200 * t266;
	t246 = -0.2e1 * t233 * t265;
	t236 = -t216 * t178 - t218 * t261;
	t234 = -t216 * t239 + t218 * t255;
	t182 = t216 * t255 + t218 * t239;
	t228 = qJD(1) * t238 + t219 * t198 + t217 * t207;
	t197 = t215 * t229;
	t163 = 0.1e1 / t165;
	t154 = 0.1e1 / t156;
	t150 = t276 * t191;
	t146 = t235 * t250 + (0.2e1 * t257 * t259 + t228 * t200 + (-t169 * t205 + t188 * t197 - t196 * t239) * t201) * t173;
	t1 = [t191 * t247 + (-t166 * t200 - t191 * t258) * t173, t146, 0, 0, 0, 0; t188 * t157 * t253 + (t169 * t157 + (-t147 * t188 - t150 * t166) * t158) * t154 + ((t150 * t252 - t276 * t264) * t154 + (t150 * t253 + ((-t151 * t173 * t260 + t250) * t263 + (t188 * t247 + t151 + (-t151 + t237) * t173) * t262) * t154) * t158) * t191, 0.2e1 * (t157 * t238 - t158 * t268) * t270 + (t167 * t157 + t252 * t268 + (t238 * t147 - t148 * t166 + (t146 * t188 - t152 * t169 - t197 + (t152 * t204 + t239) * t151) * t262 + (-t146 * t204 + t152 * t196 + t228 + (t152 * t188 - t205) * t151) * t263) * t158) * t154, 0, 0, 0, 0; (t178 * t234 - t182 * t261) * t251 + ((qJD(4) * t182 + t216 * t228 + t218 * t249) * t178 + t182 * t246 + (t234 * t162 + (qJD(4) * t234 - t216 * t249 + t218 * t228) * t233 - t182 * t161) * t179) * t163, t236 * t191 * t251 + (t236 * t166 + ((qJD(4) * t178 + t246) * t218 + (-t161 * t218 + (-t162 - t278) * t216) * t179) * t191) * t163, 0, -0.2e1 * t267 - 0.2e1 * (t161 * t179 * t163 - (-t163 * t265 - t179 * t267) * t233) * t233, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:41
	% DurationCPUTime: 2.48s
	% Computational Cost: add. (7685->143), mult. (22101->293), div. (702->12), fcn. (28611->15), ass. (0->125)
	t276 = sin(pkin(11));
	t278 = cos(pkin(11));
	t281 = sin(qJ(2));
	t284 = cos(qJ(2));
	t265 = t276 * t281 - t284 * t278;
	t279 = cos(pkin(6));
	t262 = t265 * t279;
	t257 = qJD(2) * t262;
	t302 = t276 * t284 + t278 * t281;
	t264 = t302 * qJD(2);
	t282 = sin(qJ(1));
	t285 = cos(qJ(1));
	t345 = sin(pkin(6));
	t312 = t285 * t345;
	t263 = t302 * t279;
	t327 = t282 * t263;
	t350 = -t282 * t264 - qJD(1) * t327 + (-qJD(1) * t265 - t257) * t285 - qJD(4) * t312;
	t280 = sin(qJ(4));
	t283 = cos(qJ(4));
	t303 = -t265 * t285 - t327;
	t314 = t282 * t345;
	t242 = t280 * t314 + t283 * t303;
	t246 = t282 * t262 - t285 * t302;
	t275 = sin(pkin(12));
	t277 = cos(pkin(12));
	t217 = t242 * t275 + t246 * t277;
	t244 = t263 * t285 - t282 * t265;
	t292 = -t244 * qJD(1) + t282 * t257 - t264 * t285;
	t296 = -t280 * t303 + t283 * t314;
	t311 = qJD(1) * t345;
	t307 = t285 * t311;
	t207 = t296 * qJD(4) + t280 * t307 + t283 * t292;
	t243 = -t262 * t285 - t282 * t302;
	t258 = t279 * t264;
	t297 = t265 * qJD(2);
	t222 = t243 * qJD(1) - t282 * t258 - t285 * t297;
	t202 = t207 * t277 + t222 * t275;
	t218 = t242 * t277 - t246 * t275;
	t212 = 0.1e1 / t218;
	t213 = 0.1e1 / t218 ^ 2;
	t340 = t202 * t212 * t213;
	t310 = 0.2e1 * t217 * t340;
	t236 = t244 * t280 + t283 * t312;
	t313 = t284 * t345;
	t315 = t281 * t345;
	t293 = -t276 * t313 - t278 * t315;
	t251 = -t279 * t283 - t280 * t293;
	t231 = atan2(-t236, t251);
	t226 = sin(t231);
	t227 = cos(t231);
	t200 = -t226 * t236 + t227 * t251;
	t198 = 0.1e1 / t200 ^ 2;
	t235 = t296 ^ 2;
	t196 = t198 * t235 + 0.1e1;
	t206 = t242 * qJD(4) + t280 * t292 - t283 * t307;
	t339 = t206 * t198;
	t234 = t236 ^ 2;
	t249 = 0.1e1 / t251 ^ 2;
	t230 = t234 * t249 + 0.1e1;
	t228 = 0.1e1 / t230;
	t305 = t282 * t311;
	t324 = qJD(4) * t283;
	t208 = t244 * t324 + t350 * t280 - t283 * t305;
	t252 = t279 * t280 - t283 * t293;
	t260 = -t276 * t315 + t278 * t313;
	t256 = t260 * qJD(2);
	t232 = t252 * qJD(4) + t256 * t280;
	t248 = 0.1e1 / t251;
	t331 = t236 * t249;
	t301 = -t208 * t248 + t232 * t331;
	t190 = t301 * t228;
	t304 = -t226 * t251 - t227 * t236;
	t185 = t304 * t190 - t208 * t226 + t227 * t232;
	t197 = 0.1e1 / t200;
	t199 = t197 * t198;
	t343 = t185 * t199;
	t323 = 0.2e1 * (-t235 * t343 - t296 * t339) / t196 ^ 2;
	t349 = t232 * t249;
	t299 = -t243 * t248 + t260 * t331;
	t348 = t280 * t299;
	t209 = (-qJD(4) * t244 + t305) * t280 + t350 * t283;
	t347 = -0.2e1 * t236;
	t346 = -0.2e1 * t296;
	t333 = t248 * t349;
	t342 = (t208 * t331 - t234 * t333) / t230 ^ 2;
	t341 = t198 * t296;
	t338 = t212 * t275;
	t337 = t213 * t217;
	t336 = t217 * t277;
	t335 = t226 * t296;
	t334 = t227 * t296;
	t332 = t236 * t248;
	t330 = t246 * t280;
	t329 = t246 * t283;
	t201 = t207 * t275 - t222 * t277;
	t211 = t217 ^ 2;
	t205 = t211 * t213 + 0.1e1;
	t322 = 0.2e1 * (t201 * t337 - t211 * t340) / t205 ^ 2;
	t321 = -0.2e1 * t342;
	t320 = t199 * t346;
	t319 = t248 * t342;
	t318 = t198 * t335;
	t317 = t198 * t334;
	t309 = t333 * t347;
	t238 = t244 * t283 - t280 * t312;
	t300 = -t238 * t248 + t252 * t331;
	t298 = -qJD(4) * t330 - t222 * t283;
	t295 = -t226 + (t227 * t332 + t226) * t228;
	t255 = t293 * qJD(2);
	t233 = -t251 * qJD(4) + t256 * t283;
	t224 = t246 * qJD(1) - t285 * t258 + t282 * t297;
	t220 = t275 * t303 + t277 * t329;
	t219 = t275 * t329 - t277 * t303;
	t216 = -t238 * t277 + t243 * t275;
	t215 = -t238 * t275 - t243 * t277;
	t203 = 0.1e1 / t205;
	t194 = 0.1e1 / t196;
	t193 = t228 * t348;
	t192 = t300 * t228;
	t189 = t295 * t296;
	t187 = (-t226 * t243 + t227 * t260) * t280 + t304 * t193;
	t186 = t304 * t192 - t226 * t238 + t227 * t252;
	t184 = t300 * t321 + (t252 * t309 - t209 * t248 + (t208 * t252 + t232 * t238 + t233 * t236) * t249) * t228;
	t182 = t321 * t348 + (t299 * t324 + (t260 * t309 - t224 * t248 + (t208 * t260 + t232 * t243 + t236 * t255) * t249) * t280) * t228;
	t1 = [t319 * t346 + (-t206 * t248 - t296 * t349) * t228, t182, 0, t184, 0, 0; t236 * t197 * t323 + (-t208 * t197 + (t185 * t236 + t189 * t206) * t198) * t194 - (-t189 * t198 * t323 + (-0.2e1 * t189 * t343 + (-t190 * t228 * t332 + t321) * t318 + (t319 * t347 - t190 + (t190 - t301) * t228) * t317 - t295 * t339) * t194) * t296, (-t187 * t341 - t197 * t330) * t323 + (-t187 * t339 + (-t222 * t280 + t246 * t324) * t197 + (t187 * t320 - t198 * t330) * t185 + (t260 * t324 - t182 * t236 - t193 * t208 + t255 * t280 + (-t193 * t251 - t243 * t280) * t190) * t317 + (-t243 * t324 - t182 * t251 - t193 * t232 - t224 * t280 + (t193 * t236 - t260 * t280) * t190) * t318) * t194, 0, (-t186 * t341 - t197 * t242) * t323 + (t186 * t185 * t320 + t207 * t197 + (-t242 * t185 - t186 * t206 + (-t184 * t236 - t192 * t208 + t233 + (-t192 * t251 - t238) * t190) * t334 + (-t184 * t251 - t192 * t232 - t209 + (t192 * t236 - t252) * t190) * t335) * t198) * t194, 0, 0; (-t212 * t215 + t216 * t337) * t322 + ((-t209 * t275 - t224 * t277) * t212 + t216 * t310 + (-t215 * t202 - (-t209 * t277 + t224 * t275) * t217 - t216 * t201) * t213) * t203, (-t212 * t219 + t220 * t337) * t322 + ((t298 * t275 - t277 * t292) * t212 + t220 * t310 + (-t219 * t202 - (t275 * t292 + t298 * t277) * t217 - t220 * t201) * t213) * t203, 0, -(-t213 * t336 + t338) * t296 * t322 + (t296 * t277 * t310 - t206 * t338 + (t206 * t336 - (t201 * t277 + t202 * t275) * t296) * t213) * t203, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:39
	% EndTime: 2019-10-10 10:11:41
	% DurationCPUTime: 2.68s
	% Computational Cost: add. (8842->154), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->131)
	t320 = sin(pkin(11));
	t321 = cos(pkin(11));
	t324 = sin(qJ(2));
	t326 = cos(qJ(2));
	t307 = t324 * t320 - t326 * t321;
	t322 = cos(pkin(6));
	t304 = t307 * t322;
	t299 = qJD(2) * t304;
	t348 = t326 * t320 + t324 * t321;
	t306 = t348 * qJD(2);
	t327 = cos(qJ(1));
	t392 = sin(pkin(6));
	t357 = t327 * t392;
	t305 = t348 * t322;
	t393 = sin(qJ(1));
	t360 = t393 * t305;
	t398 = -t393 * t306 - qJD(1) * t360 + (-qJD(1) * t307 - t299) * t327 - qJD(4) * t357;
	t323 = sin(qJ(4));
	t325 = cos(qJ(4));
	t341 = -t327 * t305 + t393 * t307;
	t277 = -t323 * t341 + t325 * t357;
	t358 = t326 * t392;
	t359 = t324 * t392;
	t335 = -t320 * t358 - t321 * t359;
	t293 = -t322 * t325 - t323 * t335;
	t272 = atan2(-t277, t293);
	t267 = sin(t272);
	t268 = cos(t272);
	t243 = -t267 * t277 + t268 * t293;
	t241 = 0.1e1 / t243 ^ 2;
	t342 = -t327 * t307 - t360;
	t353 = t393 * t392;
	t338 = -t323 * t342 + t325 * t353;
	t276 = t338 ^ 2;
	t239 = t276 * t241 + 0.1e1;
	t283 = t323 * t353 + t325 * t342;
	t334 = t341 * qJD(1) + t393 * t299 - t327 * t306;
	t351 = qJD(1) * t357;
	t247 = t283 * qJD(4) + t323 * t334 - t325 * t351;
	t385 = t247 * t241;
	t275 = t277 ^ 2;
	t291 = 0.1e1 / t293 ^ 2;
	t271 = t275 * t291 + 0.1e1;
	t269 = 0.1e1 / t271;
	t347 = qJD(1) * t353;
	t371 = qJD(4) * t325;
	t249 = t398 * t323 - t325 * t347 - t341 * t371;
	t294 = t322 * t323 - t325 * t335;
	t302 = -t320 * t359 + t321 * t358;
	t298 = t302 * qJD(2);
	t273 = t294 * qJD(4) + t298 * t323;
	t290 = 0.1e1 / t293;
	t379 = t277 * t291;
	t346 = -t249 * t290 + t273 * t379;
	t231 = t346 * t269;
	t349 = -t267 * t293 - t268 * t277;
	t226 = t349 * t231 - t267 * t249 + t268 * t273;
	t240 = 0.1e1 / t243;
	t242 = t240 * t241;
	t390 = t226 * t242;
	t369 = 0.2e1 * (-t276 * t390 - t338 * t385) / t239 ^ 2;
	t397 = t273 * t291;
	t285 = -t327 * t304 - t348 * t393;
	t343 = -t285 * t290 + t302 * t379;
	t396 = t323 * t343;
	t250 = (qJD(4) * t341 + t347) * t323 + t398 * t325;
	t288 = t393 * t304 - t327 * t348;
	t319 = pkin(12) + qJ(6);
	t317 = sin(t319);
	t318 = cos(t319);
	t259 = t283 * t318 - t288 * t317;
	t253 = 0.1e1 / t259;
	t254 = 0.1e1 / t259 ^ 2;
	t395 = -0.2e1 * t277;
	t394 = -0.2e1 * t338;
	t248 = t338 * qJD(4) + t323 * t351 + t325 * t334;
	t300 = t322 * t306;
	t340 = t307 * qJD(2);
	t263 = t285 * qJD(1) - t393 * t300 - t327 * t340;
	t234 = t259 * qJD(6) + t248 * t317 - t263 * t318;
	t258 = t283 * t317 + t288 * t318;
	t252 = t258 ^ 2;
	t246 = t252 * t254 + 0.1e1;
	t384 = t254 * t258;
	t370 = qJD(6) * t258;
	t235 = t248 * t318 + t263 * t317 - t370;
	t387 = t235 * t253 * t254;
	t389 = (t234 * t384 - t252 * t387) / t246 ^ 2;
	t381 = t290 * t397;
	t388 = (t249 * t379 - t275 * t381) / t271 ^ 2;
	t386 = t241 * t338;
	t383 = t267 * t338;
	t382 = t268 * t338;
	t380 = t277 * t290;
	t378 = t288 * t323;
	t377 = t288 * t325;
	t376 = t317 * t253;
	t375 = t318 * t258;
	t368 = -0.2e1 * t389;
	t367 = 0.2e1 * t389;
	t366 = -0.2e1 * t388;
	t365 = t242 * t394;
	t364 = t290 * t388;
	t363 = t241 * t383;
	t362 = t241 * t382;
	t361 = t258 * t387;
	t356 = 0.2e1 * t361;
	t355 = t381 * t395;
	t279 = -t323 * t357 - t325 * t341;
	t350 = qJD(6) * t377 - t334;
	t257 = -t279 * t318 + t285 * t317;
	t256 = -t279 * t317 - t285 * t318;
	t345 = t254 * t375 - t376;
	t344 = -t279 * t290 + t294 * t379;
	t339 = -t267 + (t268 * t380 + t267) * t269;
	t336 = -qJD(4) * t378 + qJD(6) * t342 - t263 * t325;
	t297 = t335 * qJD(2);
	t274 = -t293 * qJD(4) + t298 * t325;
	t265 = t288 * qJD(1) - t327 * t300 + t393 * t340;
	t261 = t317 * t342 + t318 * t377;
	t260 = t317 * t377 - t318 * t342;
	t244 = 0.1e1 / t246;
	t237 = 0.1e1 / t239;
	t236 = t269 * t396;
	t233 = t344 * t269;
	t230 = t339 * t338;
	t228 = (-t267 * t285 + t268 * t302) * t323 + t349 * t236;
	t227 = t349 * t233 - t267 * t279 + t268 * t294;
	t225 = t344 * t366 + (t294 * t355 - t250 * t290 + (t249 * t294 + t273 * t279 + t274 * t277) * t291) * t269;
	t223 = t366 * t396 + (t343 * t371 + (t302 * t355 - t265 * t290 + (t249 * t302 + t273 * t285 + t277 * t297) * t291) * t323) * t269;
	t1 = [t364 * t394 + (-t247 * t290 - t338 * t397) * t269, t223, 0, t225, 0, 0; t277 * t240 * t369 + (-t249 * t240 + (t226 * t277 + t230 * t247) * t241) * t237 - (-t230 * t241 * t369 + (-0.2e1 * t230 * t390 + (-t231 * t269 * t380 + t366) * t363 + (t364 * t395 - t231 + (t231 - t346) * t269) * t362 - t339 * t385) * t237) * t338, (-t228 * t386 - t240 * t378) * t369 + (-t228 * t385 + (-t263 * t323 + t288 * t371) * t240 + (t228 * t365 - t241 * t378) * t226 + (t302 * t371 - t223 * t277 - t236 * t249 + t297 * t323 + (-t236 * t293 - t285 * t323) * t231) * t362 + (-t285 * t371 - t223 * t293 - t236 * t273 - t265 * t323 + (t236 * t277 - t302 * t323) * t231) * t363) * t237, 0, (-t227 * t386 - t240 * t283) * t369 + (t227 * t226 * t365 + t248 * t240 + (-t283 * t226 - t227 * t247 + (-t225 * t277 - t233 * t249 + t274 + (-t233 * t293 - t279) * t231) * t382 + (-t225 * t293 - t233 * t273 - t250 + (t233 * t277 - t294) * t231) * t383) * t241) * t237, 0, 0; (-t253 * t256 + t257 * t384) * t367 + ((t257 * qJD(6) - t250 * t317 - t265 * t318) * t253 + t257 * t356 + (-t256 * t235 - (-t256 * qJD(6) - t250 * t318 + t265 * t317) * t258 - t257 * t234) * t254) * t244, (-t253 * t260 + t261 * t384) * t367 + (t261 * t356 + t350 * t253 * t318 + t336 * t376 + (t350 * t258 * t317 - t261 * t234 - t260 * t235 - t336 * t375) * t254) * t244, 0, -t345 * t338 * t368 + (t345 * t247 - ((-qJD(6) * t253 - 0.2e1 * t361) * t318 + (t234 * t318 + (t235 - t370) * t317) * t254) * t338) * t244, 0, t368 + 0.2e1 * (t234 * t254 * t244 + (-t244 * t387 - t254 * t389) * t258) * t258;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end