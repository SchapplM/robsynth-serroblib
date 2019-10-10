% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR4
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
%   Wie in S6RRPRRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:49
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
	% StartTime: 2019-10-10 10:55:48
	% EndTime: 2019-10-10 10:55:49
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (454->46), mult. (1493->122), div. (133->12), fcn. (1882->11), ass. (0->62)
	t133 = sin(qJ(1));
	t134 = cos(qJ(1));
	t131 = cos(pkin(6));
	t129 = sin(pkin(12));
	t132 = sin(qJ(2));
	t161 = cos(pkin(12));
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
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:50
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (3306->95), mult. (9869->202), div. (448->12), fcn. (12731->13), ass. (0->96)
	t216 = sin(qJ(4));
	t218 = cos(qJ(4));
	t271 = sin(pkin(12));
	t273 = cos(pkin(6));
	t241 = t273 * t271;
	t272 = cos(pkin(12));
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
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:50
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (3975->97), mult. (10490->203), div. (466->12), fcn. (13487->13), ass. (0->100)
	t245 = sin(pkin(6));
	t299 = sin(pkin(12));
	t300 = cos(pkin(12));
	t302 = sin(qJ(2));
	t303 = cos(qJ(2));
	t258 = t303 * t299 + t302 * t300;
	t231 = t258 * t245;
	t222 = qJD(2) * t231;
	t234 = t302 * t299 - t303 * t300;
	t230 = t234 * t245;
	t227 = 0.1e1 / t230 ^ 2;
	t286 = t222 * t227;
	t301 = cos(pkin(6));
	t269 = t301 * t299;
	t270 = t301 * t300;
	t305 = t302 * t269 - t303 * t270;
	t257 = t234 * qJD(2);
	t246 = sin(qJ(1));
	t247 = cos(qJ(1));
	t214 = -t246 * t258 - t247 * t305;
	t202 = atan2(t214, t230);
	t197 = sin(t202);
	t198 = cos(t202);
	t211 = t214 ^ 2;
	t201 = t211 * t227 + 0.1e1;
	t199 = 0.1e1 / t201;
	t226 = 0.1e1 / t230;
	t288 = t214 * t226;
	t304 = (t198 * t288 - t197) * t199 + t197;
	t186 = t197 * t214 + t198 * t230;
	t183 = 0.1e1 / t186;
	t244 = qJ(4) + qJ(5);
	t241 = sin(t244);
	t242 = cos(t244);
	t232 = t303 * t269 + t302 * t270;
	t266 = t246 * t232 + t247 * t234;
	t284 = t245 * t246;
	t210 = t241 * t284 - t242 * t266;
	t204 = 0.1e1 / t210;
	t184 = 0.1e1 / t186 ^ 2;
	t205 = 0.1e1 / t210 ^ 2;
	t254 = t246 * t305;
	t217 = -t247 * t258 + t254;
	t212 = t217 ^ 2;
	t182 = t212 * t184 + 0.1e1;
	t225 = t232 * qJD(2);
	t192 = t214 * qJD(1) - t246 * t225 - t247 * t257;
	t292 = t192 * t184;
	t282 = qJD(1) * t247;
	t195 = qJD(1) * t254 - t247 * t225 + t246 * t257 - t258 * t282;
	t265 = t195 * t226 - t214 * t286;
	t177 = t265 * t199;
	t268 = -t197 * t230 + t198 * t214;
	t173 = t268 * t177 + t197 * t195 + t198 * t222;
	t297 = t173 * t183 * t184;
	t298 = (-t212 * t297 - t217 * t292) / t182 ^ 2;
	t267 = -t247 * t232 + t246 * t234;
	t287 = t214 * t231;
	t263 = -t226 * t267 + t227 * t287;
	t178 = t263 * t199;
	t174 = -t268 * t178 + t197 * t267 + t198 * t231;
	t296 = t174 * t217;
	t243 = qJD(4) + qJD(5);
	t261 = t243 * t266 + t245 * t282;
	t224 = t305 * qJD(2);
	t233 = t258 * qJD(2);
	t193 = t267 * qJD(1) + t246 * t224 - t247 * t233;
	t273 = t243 * t284 + t193;
	t187 = t273 * t241 - t261 * t242;
	t209 = -t241 * t266 - t242 * t284;
	t203 = t209 ^ 2;
	t191 = t203 * t205 + 0.1e1;
	t289 = t205 * t209;
	t188 = t261 * t241 + t273 * t242;
	t293 = t188 * t204 * t205;
	t295 = (t187 * t289 - t203 * t293) / t191 ^ 2;
	t285 = t226 * t286;
	t294 = (t214 * t227 * t195 - t211 * t285) / t201 ^ 2;
	t291 = t197 * t217;
	t290 = t198 * t217;
	t283 = t245 * t247;
	t281 = -0.2e1 * t298;
	t280 = -0.2e1 * t297;
	t279 = 0.2e1 * t295;
	t278 = 0.2e1 * t294;
	t277 = -0.2e1 * t226 * t294;
	t276 = 0.2e1 * t209 * t293;
	t256 = t266 * qJD(1) + t247 * t224 + t246 * t233;
	t272 = t243 * t283 + t256;
	t264 = -t241 * t204 + t242 * t289;
	t262 = qJD(1) * t284 + t243 * t267;
	t223 = t245 * t257;
	t208 = t241 * t283 + t242 * t267;
	t207 = t241 * t267 - t242 * t283;
	t189 = 0.1e1 / t191;
	t180 = 0.1e1 / t182;
	t176 = t304 * t217;
	t172 = t263 * t278 + (0.2e1 * t285 * t287 + t256 * t226 + (-t195 * t231 + t214 * t223 - t222 * t267) * t227) * t199;
	t170 = -0.2e1 * t295 + 0.2e1 * (t187 * t205 * t189 + (-t189 * t293 - t205 * t295) * t209) * t209;
	t1 = [t217 * t277 + (-t192 * t226 - t217 * t286) * t199, t172, 0, 0, 0, 0; t214 * t183 * t281 + (t195 * t183 + (-t173 * t214 - t176 * t192) * t184) * t180 + ((t176 * t280 - t304 * t292) * t180 + (t176 * t281 + ((-t177 * t199 * t288 + t278) * t291 + (t214 * t277 + t177 + (-t177 + t265) * t199) * t290) * t180) * t184) * t217, 0.2e1 * (t183 * t266 - t184 * t296) * t298 + (t193 * t183 + t280 * t296 + (t266 * t173 - t174 * t192 + (t172 * t214 - t178 * t195 - t223 + (t178 * t230 + t267) * t177) * t290 + (-t172 * t230 + t178 * t222 + t256 + (t178 * t214 - t231) * t177) * t291) * t184) * t180, 0, 0, 0, 0; (-t204 * t207 + t208 * t289) * t279 + ((t272 * t241 + t262 * t242) * t204 + t208 * t276 + (-t207 * t188 - (-t262 * t241 + t272 * t242) * t209 - t208 * t187) * t205) * t189, t264 * t217 * t279 + (t264 * t192 + ((t204 * t243 + t276) * t242 + (-t187 * t242 + (t209 * t243 - t188) * t241) * t205) * t217) * t189, 0, t170, t170, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:55:49
	% EndTime: 2019-10-10 10:55:52
	% DurationCPUTime: 3.31s
	% Computational Cost: add. (17152->158), mult. (32251->311), div. (989->12), fcn. (41704->15), ass. (0->139)
	t359 = qJ(4) + qJ(5);
	t356 = sin(t359);
	t362 = cos(pkin(6));
	t360 = sin(pkin(12));
	t361 = cos(pkin(12));
	t364 = sin(qJ(2));
	t366 = cos(qJ(2));
	t388 = t360 * t366 + t361 * t364;
	t344 = t388 * t362;
	t348 = t360 * t364 - t366 * t361;
	t367 = cos(qJ(1));
	t436 = sin(qJ(1));
	t381 = -t367 * t344 + t436 * t348;
	t357 = cos(t359);
	t435 = sin(pkin(6));
	t398 = t367 * t435;
	t393 = t357 * t398;
	t316 = -t356 * t381 + t393;
	t399 = t366 * t435;
	t400 = t364 * t435;
	t375 = -t360 * t399 - t361 * t400;
	t332 = -t356 * t375 - t357 * t362;
	t304 = atan2(-t316, t332);
	t299 = sin(t304);
	t300 = cos(t304);
	t282 = -t299 * t316 + t300 * t332;
	t280 = 0.1e1 / t282 ^ 2;
	t401 = t436 * t344;
	t382 = -t367 * t348 - t401;
	t392 = t436 * t435;
	t321 = t356 * t382 - t357 * t392;
	t315 = t321 ^ 2;
	t278 = t280 * t315 + 0.1e1;
	t358 = qJD(4) + qJD(5);
	t343 = t348 * t362;
	t338 = qJD(2) * t343;
	t345 = t388 * qJD(2);
	t374 = t381 * qJD(1) + t436 * t338 - t367 * t345;
	t379 = t358 * t392 + t374;
	t391 = qJD(1) * t398;
	t415 = t357 * t358;
	t286 = t379 * t356 - t357 * t391 + t382 * t415;
	t428 = t280 * t321;
	t314 = t316 ^ 2;
	t330 = 0.1e1 / t332 ^ 2;
	t303 = t314 * t330 + 0.1e1;
	t301 = 0.1e1 / t303;
	t311 = -t436 * t345 - qJD(1) * t401 + (-qJD(1) * t348 - t338) * t367;
	t353 = t356 * t398;
	t387 = qJD(1) * t392;
	t288 = t311 * t356 - t358 * t353 - t357 * t387 - t381 * t415;
	t341 = -t360 * t400 + t361 * t399;
	t394 = t341 * qJD(2) + t358 * t362;
	t312 = t394 * t356 - t375 * t415;
	t329 = 0.1e1 / t332;
	t420 = t316 * t330;
	t386 = -t288 * t329 + t312 * t420;
	t270 = t386 * t301;
	t389 = -t299 * t332 - t300 * t316;
	t265 = t389 * t270 - t288 * t299 + t300 * t312;
	t279 = 0.1e1 / t282;
	t281 = t279 * t280;
	t433 = t265 * t281;
	t410 = 0.2e1 * (t286 * t428 - t315 * t433) / t278 ^ 2;
	t440 = t312 * t330;
	t324 = -t367 * t343 - t388 * t436;
	t383 = -t324 * t329 + t341 * t420;
	t439 = t356 * t383;
	t289 = (t358 * t381 + t387) * t356 + t311 * t357 - t358 * t393;
	t322 = t356 * t392 + t357 * t382;
	t365 = cos(qJ(6));
	t327 = t436 * t343 - t367 * t388;
	t363 = sin(qJ(6));
	t418 = t327 * t363;
	t298 = t322 * t365 - t418;
	t292 = 0.1e1 / t298;
	t293 = 0.1e1 / t298 ^ 2;
	t438 = -0.2e1 * t316;
	t437 = 0.2e1 * t321;
	t287 = t379 * t357 + (-t358 * t382 + t391) * t356;
	t339 = t362 * t345;
	t380 = t348 * qJD(2);
	t308 = t324 * qJD(1) - t436 * t339 - t367 * t380;
	t274 = t298 * qJD(6) + t287 * t363 - t308 * t365;
	t417 = t327 * t365;
	t297 = t322 * t363 + t417;
	t291 = t297 ^ 2;
	t285 = t291 * t293 + 0.1e1;
	t426 = t293 * t297;
	t411 = qJD(6) * t297;
	t275 = t287 * t365 + t308 * t363 - t411;
	t430 = t275 * t292 * t293;
	t432 = (t274 * t426 - t291 * t430) / t285 ^ 2;
	t422 = t329 * t440;
	t431 = (t288 * t420 - t314 * t422) / t303 ^ 2;
	t429 = t280 * t286;
	t427 = t292 * t363;
	t425 = t297 * t365;
	t424 = t299 * t321;
	t423 = t300 * t321;
	t421 = t316 * t329;
	t419 = t327 * t356;
	t416 = t356 * t358;
	t409 = -0.2e1 * t432;
	t408 = 0.2e1 * t432;
	t407 = -0.2e1 * t431;
	t406 = t281 * t437;
	t405 = t329 * t431;
	t404 = t297 * t430;
	t403 = t280 * t424;
	t402 = t280 * t423;
	t397 = 0.2e1 * t404;
	t396 = t422 * t438;
	t318 = -t357 * t381 - t353;
	t390 = qJD(6) * t327 * t357 - t374;
	t296 = -t318 * t365 + t324 * t363;
	t295 = -t318 * t363 - t324 * t365;
	t385 = t293 * t425 - t427;
	t333 = t356 * t362 - t357 * t375;
	t384 = -t318 * t329 + t333 * t420;
	t378 = -t299 + (t300 * t421 + t299) * t301;
	t376 = qJD(6) * t382 - t308 * t357 - t327 * t416;
	t336 = t375 * qJD(2);
	t313 = t394 * t357 + t375 * t416;
	t310 = t327 * qJD(1) - t367 * t339 + t436 * t380;
	t306 = t357 * t417 + t363 * t382;
	t305 = t357 * t418 - t365 * t382;
	t283 = 0.1e1 / t285;
	t276 = 0.1e1 / t278;
	t273 = t301 * t439;
	t272 = t384 * t301;
	t269 = t378 * t321;
	t267 = (-t299 * t324 + t300 * t341) * t356 + t389 * t273;
	t266 = t389 * t272 - t299 * t318 + t300 * t333;
	t263 = t384 * t407 + (t333 * t396 - t289 * t329 + (t288 * t333 + t312 * t318 + t313 * t316) * t330) * t301;
	t262 = t407 * t439 + (t383 * t415 + (t341 * t396 - t310 * t329 + (t288 * t341 + t312 * t324 + t316 * t336) * t330) * t356) * t301;
	t261 = t385 * t321 * t409 + (t385 * t286 + ((-qJD(6) * t292 - 0.2e1 * t404) * t365 + (t274 * t365 + (t275 - t411) * t363) * t293) * t321) * t283;
	t260 = (t266 * t428 - t279 * t322) * t410 + (t266 * t265 * t406 + t287 * t279 + (-t322 * t265 - t266 * t286 - (-t263 * t316 - t272 * t288 + t313 + (-t272 * t332 - t318) * t270) * t423 - (-t263 * t332 - t272 * t312 - t289 + (t272 * t316 - t333) * t270) * t424) * t280) * t276;
	t1 = [t405 * t437 + (-t286 * t329 + t321 * t440) * t301, t262, 0, t263, t263, 0; t316 * t279 * t410 + (-t288 * t279 + (t265 * t316 - t269 * t286) * t280) * t276 + (t269 * t280 * t410 + (0.2e1 * t269 * t433 - (-t270 * t301 * t421 + t407) * t403 - (t405 * t438 - t270 + (t270 - t386) * t301) * t402 - t378 * t429) * t276) * t321, (t267 * t428 - t279 * t419) * t410 + (-t267 * t429 + (-t308 * t356 + t327 * t415) * t279 + (t267 * t406 - t280 * t419) * t265 - (t341 * t415 - t262 * t316 - t273 * t288 + t336 * t356 + (-t273 * t332 - t324 * t356) * t270) * t402 - (-t324 * t415 - t262 * t332 - t273 * t312 - t310 * t356 + (t273 * t316 - t341 * t356) * t270) * t403) * t276, 0, t260, t260, 0; (-t292 * t295 + t296 * t426) * t408 + ((t296 * qJD(6) - t289 * t363 - t310 * t365) * t292 + t296 * t397 + (-t295 * t275 - (-t295 * qJD(6) - t289 * t365 + t310 * t363) * t297 - t296 * t274) * t293) * t283, (-t292 * t305 + t306 * t426) * t408 + (t306 * t397 + t390 * t292 * t365 + t376 * t427 + (t390 * t297 * t363 - t306 * t274 - t305 * t275 - t376 * t425) * t293) * t283, 0, t261, t261, t409 + 0.2e1 * (t274 * t283 * t293 + (-t283 * t430 - t293 * t432) * t297) * t297;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end