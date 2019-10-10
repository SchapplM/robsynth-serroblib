% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR6
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
%   Wie in S6RRPRPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:33
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
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:33
	% DurationCPUTime: 0.46s
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
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:34
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
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:35
	% DurationCPUTime: 2.16s
	% Computational Cost: add. (6793->127), mult. (19653->266), div. (691->12), fcn. (25511->13), ass. (0->115)
	t241 = sin(pkin(11));
	t242 = cos(pkin(11));
	t245 = sin(qJ(2));
	t248 = cos(qJ(2));
	t231 = t245 * t241 - t248 * t242;
	t243 = cos(pkin(6));
	t228 = t231 * t243;
	t223 = qJD(2) * t228;
	t264 = t248 * t241 + t245 * t242;
	t230 = t264 * qJD(2);
	t246 = sin(qJ(1));
	t249 = cos(qJ(1));
	t306 = sin(pkin(6));
	t273 = t249 * t306;
	t229 = t264 * t243;
	t289 = t246 * t229;
	t311 = -t246 * t230 - qJD(1) * t289 + (-qJD(1) * t231 - t223) * t249 - qJD(4) * t273;
	t209 = t249 * t229 - t246 * t231;
	t244 = sin(qJ(4));
	t247 = cos(qJ(4));
	t198 = t209 * t244 + t247 * t273;
	t274 = t248 * t306;
	t276 = t245 * t306;
	t256 = -t241 * t274 - t242 * t276;
	t217 = -t243 * t247 - t244 * t256;
	t192 = atan2(-t198, t217);
	t187 = sin(t192);
	t188 = cos(t192);
	t175 = -t187 * t198 + t188 * t217;
	t173 = 0.1e1 / t175 ^ 2;
	t265 = -t249 * t231 - t289;
	t275 = t246 * t306;
	t259 = -t244 * t265 + t247 * t275;
	t196 = t259 ^ 2;
	t171 = t196 * t173 + 0.1e1;
	t184 = -t209 * qJD(1) + t246 * t223 - t249 * t230;
	t204 = t244 * t275 + t247 * t265;
	t272 = qJD(1) * t306;
	t267 = t249 * t272;
	t176 = t204 * qJD(4) + t184 * t244 - t247 * t267;
	t300 = t176 * t173;
	t195 = t198 ^ 2;
	t215 = 0.1e1 / t217 ^ 2;
	t191 = t195 * t215 + 0.1e1;
	t189 = 0.1e1 / t191;
	t268 = t246 * t272;
	t286 = qJD(4) * t247;
	t178 = t209 * t286 + t311 * t244 - t247 * t268;
	t218 = t243 * t244 - t247 * t256;
	t226 = -t241 * t276 + t242 * t274;
	t222 = t226 * qJD(2);
	t193 = t218 * qJD(4) + t222 * t244;
	t214 = 0.1e1 / t217;
	t294 = t198 * t215;
	t263 = -t178 * t214 + t193 * t294;
	t164 = t263 * t189;
	t266 = -t187 * t217 - t188 * t198;
	t160 = t266 * t164 - t187 * t178 + t188 * t193;
	t172 = 0.1e1 / t175;
	t174 = t172 * t173;
	t304 = t160 * t174;
	t285 = 0.2e1 * (-t196 * t304 - t259 * t300) / t171 ^ 2;
	t310 = t193 * t215;
	t208 = -t249 * t228 - t246 * t264;
	t261 = -t208 * t214 + t226 * t294;
	t309 = t244 * t261;
	t179 = (-qJD(4) * t209 + t268) * t244 + t311 * t247;
	t211 = t246 * t228 - t249 * t264;
	t205 = 0.1e1 / t211;
	t206 = 0.1e1 / t211 ^ 2;
	t308 = -0.2e1 * t198;
	t307 = -0.2e1 * t259;
	t177 = t259 * qJD(4) + t184 * t247 + t244 * t267;
	t197 = t204 ^ 2;
	t182 = t197 * t206 + 0.1e1;
	t224 = t243 * t230;
	t260 = t231 * qJD(2);
	t183 = t208 * qJD(1) - t246 * t224 - t249 * t260;
	t207 = t205 * t206;
	t293 = t204 * t206;
	t303 = (t197 * t207 * t183 + t177 * t293) / t182 ^ 2;
	t296 = t214 * t310;
	t302 = (t178 * t294 - t195 * t296) / t191 ^ 2;
	t301 = t173 * t259;
	t299 = t183 * t206;
	t298 = t187 * t259;
	t297 = t188 * t259;
	t295 = t198 * t214;
	t292 = t206 * t208;
	t291 = t211 * t244;
	t284 = 0.2e1 * t303;
	t283 = -0.2e1 * t302;
	t282 = t174 * t307;
	t281 = -0.2e1 * t204 * t207;
	t280 = t205 * t303;
	t279 = t214 * t302;
	t278 = t173 * t298;
	t277 = t173 * t297;
	t271 = t296 * t308;
	t200 = t209 * t247 - t244 * t273;
	t262 = -t200 * t214 + t218 * t294;
	t258 = -t187 + (t188 * t295 + t187) * t189;
	t221 = t256 * qJD(2);
	t194 = -t217 * qJD(4) + t222 * t247;
	t185 = t211 * qJD(1) - t249 * t224 + t246 * t260;
	t180 = 0.1e1 / t182;
	t169 = 0.1e1 / t171;
	t168 = t189 * t309;
	t167 = t262 * t189;
	t163 = t258 * t259;
	t162 = (-t187 * t208 + t188 * t226) * t244 + t266 * t168;
	t161 = t266 * t167 - t187 * t200 + t188 * t218;
	t159 = t262 * t283 + (t218 * t271 - t179 * t214 + (t178 * t218 + t193 * t200 + t194 * t198) * t215) * t189;
	t157 = t283 * t309 + (t261 * t286 + (t226 * t271 - t185 * t214 + (t178 * t226 + t193 * t208 + t198 * t221) * t215) * t244) * t189;
	t1 = [t279 * t307 + (-t176 * t214 - t259 * t310) * t189, t157, 0, t159, 0, 0; t198 * t172 * t285 + (-t178 * t172 + (t160 * t198 + t163 * t176) * t173) * t169 - (-t163 * t173 * t285 + (-0.2e1 * t163 * t304 + (-t164 * t189 * t295 + t283) * t278 + (t279 * t308 - t164 + (t164 - t263) * t189) * t277 - t258 * t300) * t169) * t259, (-t162 * t301 - t172 * t291) * t285 + (-t162 * t300 + (-t183 * t244 + t211 * t286) * t172 + (t162 * t282 - t173 * t291) * t160 + (t226 * t286 - t157 * t198 - t168 * t178 + t221 * t244 + (-t168 * t217 - t208 * t244) * t164) * t277 + (-t208 * t286 - t157 * t217 - t168 * t193 - t185 * t244 + (t168 * t198 - t226 * t244) * t164) * t278) * t169, 0, (-t161 * t301 - t172 * t204) * t285 + (t161 * t160 * t282 + t177 * t172 + (-t204 * t160 - t161 * t176 + (-t159 * t198 - t167 * t178 + t194 + (-t167 * t217 - t200) * t164) * t297 + (-t159 * t217 - t167 * t193 - t179 + (t167 * t198 - t218) * t164) * t298) * t173) * t169, 0, 0; t204 * t284 * t292 - 0.2e1 * t200 * t280 + (t208 * t183 * t281 - t177 * t292 + t179 * t205 - t185 * t293 + t200 * t299) * t180, (t205 * t211 * t247 + t265 * t293) * t284 + (qJD(4) * t205 * t291 + (-t177 * t265 - t184 * t204) * t206 + (t265 * t281 + (-t206 * t211 + t205) * t247) * t183) * t180, 0, 0.2e1 * t259 * t280 + (t176 * t205 - t259 * t299) * t180, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:36
	% DurationCPUTime: 2.68s
	% Computational Cost: add. (8421->154), mult. (24081->309), div. (726->12), fcn. (31139->15), ass. (0->130)
	t310 = sin(pkin(11));
	t312 = cos(pkin(11));
	t316 = sin(qJ(2));
	t320 = cos(qJ(2));
	t301 = t316 * t310 - t320 * t312;
	t313 = cos(pkin(6));
	t332 = t301 * t313;
	t294 = qJD(2) * t332;
	t337 = t320 * t310 + t316 * t312;
	t299 = t337 * t313;
	t300 = t337 * qJD(2);
	t317 = sin(qJ(1));
	t321 = cos(qJ(1));
	t359 = qJD(1) * t317;
	t311 = sin(pkin(6));
	t365 = t311 * t321;
	t388 = -t317 * t300 - t299 * t359 + (-qJD(1) * t301 - t294) * t321 - qJD(4) * t365;
	t315 = sin(qJ(4));
	t338 = -t317 * t299 - t321 * t301;
	t319 = cos(qJ(4));
	t366 = t311 * t319;
	t276 = t315 * t338 - t317 * t366;
	t283 = t317 * t332 - t321 * t337;
	t314 = sin(qJ(6));
	t318 = cos(qJ(6));
	t339 = t276 * t318 + t283 * t314;
	t387 = t339 * qJD(6);
	t281 = t321 * t299 - t317 * t301;
	t363 = t315 * t321;
	t272 = t281 * t319 - t311 * t363;
	t298 = t337 * t311;
	t289 = t298 * t319 + t313 * t315;
	t265 = atan2(-t272, t289);
	t260 = sin(t265);
	t261 = cos(t265);
	t236 = -t260 * t272 + t261 * t289;
	t234 = 0.1e1 / t236 ^ 2;
	t277 = t317 * t311 * t315 + t319 * t338;
	t270 = t277 ^ 2;
	t232 = t270 * t234 + 0.1e1;
	t328 = -t281 * qJD(1) + t317 * t294 - t321 * t300;
	t356 = qJD(4) * t319;
	t357 = qJD(4) * t315;
	t241 = t328 * t319 - t338 * t357 + (qJD(1) * t363 + t317 * t356) * t311;
	t376 = t241 * t234;
	t269 = t272 ^ 2;
	t286 = 0.1e1 / t289 ^ 2;
	t264 = t269 * t286 + 0.1e1;
	t262 = 0.1e1 / t264;
	t343 = t388 * t319;
	t348 = t311 * t359;
	t243 = -t281 * t357 + t315 * t348 + t343;
	t288 = -t298 * t315 + t313 * t319;
	t297 = t301 * t311;
	t293 = qJD(2) * t297;
	t267 = t288 * qJD(4) - t293 * t319;
	t285 = 0.1e1 / t289;
	t370 = t272 * t286;
	t336 = -t243 * t285 + t267 * t370;
	t224 = t336 * t262;
	t341 = -t260 * t289 - t261 * t272;
	t219 = t341 * t224 - t260 * t243 + t261 * t267;
	t233 = 0.1e1 / t236;
	t235 = t233 * t234;
	t381 = t219 * t235;
	t355 = 0.2e1 * (-t270 * t381 + t277 * t376) / t232 ^ 2;
	t386 = t267 * t286;
	t280 = -t317 * t337 - t321 * t332;
	t333 = -t280 * t285 - t297 * t370;
	t385 = t319 * t333;
	t368 = t283 * t318;
	t252 = t276 * t314 - t368;
	t246 = 0.1e1 / t252;
	t247 = 0.1e1 / t252 ^ 2;
	t384 = -0.2e1 * t272;
	t383 = 0.2e1 * t277;
	t347 = qJD(1) * t366;
	t240 = t277 * qJD(4) + t315 * t328 - t321 * t347;
	t295 = t313 * t300;
	t331 = t301 * qJD(2);
	t256 = t280 * qJD(1) - t317 * t295 - t321 * t331;
	t227 = t252 * qJD(6) - t240 * t318 + t256 * t314;
	t245 = t339 ^ 2;
	t239 = t245 * t247 + 0.1e1;
	t375 = t247 * t339;
	t228 = t240 * t314 + t256 * t318 + t387;
	t378 = t228 * t246 * t247;
	t380 = (-t227 * t375 - t245 * t378) / t239 ^ 2;
	t372 = t285 * t386;
	t379 = (t243 * t370 - t269 * t372) / t264 ^ 2;
	t377 = t234 * t277;
	t374 = t260 * t277;
	t373 = t261 * t277;
	t371 = t272 * t285;
	t369 = t283 * t315;
	t367 = t283 * t319;
	t364 = t314 * t339;
	t361 = t318 * t246;
	t354 = 0.2e1 * t380;
	t353 = -0.2e1 * t379;
	t352 = t235 * t383;
	t351 = t285 * t379;
	t350 = t234 * t374;
	t349 = t234 * t373;
	t345 = -0.2e1 * t339 * t378;
	t344 = t372 * t384;
	t342 = qJD(6) * t369 + t328;
	t271 = t281 * t315 + t319 * t365;
	t340 = -t271 * t318 - t280 * t314;
	t250 = -t271 * t314 + t280 * t318;
	t335 = -t247 * t364 + t361;
	t334 = t271 * t285 + t288 * t370;
	t330 = -t260 + (t261 * t371 + t260) * t262;
	t242 = t281 * t356 + t388 * t315 - t317 * t347;
	t329 = -qJD(6) * t338 - t256 * t315 + t283 * t356;
	t292 = t311 * t300;
	t266 = -t289 * qJD(4) + t293 * t315;
	t258 = t283 * qJD(1) - t321 * t295 + t317 * t331;
	t254 = t314 * t369 + t318 * t338;
	t253 = t314 * t338 - t315 * t368;
	t237 = 0.1e1 / t239;
	t230 = 0.1e1 / t232;
	t229 = t262 * t385;
	t226 = t334 * t262;
	t223 = t330 * t277;
	t221 = (-t260 * t280 - t261 * t297) * t319 + t341 * t229;
	t220 = t341 * t226 + t260 * t271 + t261 * t288;
	t218 = t334 * t353 + (t288 * t344 + t242 * t285 + (t243 * t288 + t266 * t272 - t267 * t271) * t286) * t262;
	t216 = t353 * t385 + (-t333 * t357 + (-t297 * t344 - t258 * t285 + (-t243 * t297 + t267 * t280 - t272 * t292) * t286) * t319) * t262;
	t1 = [t351 * t383 + (-t241 * t285 + t277 * t386) * t262, t216, 0, t218, 0, 0; t272 * t233 * t355 + (((qJD(4) * t281 - t348) * t315 - t343) * t233 + (t219 * t272 - t223 * t241) * t234) * t230 + (t223 * t234 * t355 + (0.2e1 * t223 * t381 - (-t224 * t262 * t371 + t353) * t350 - (t351 * t384 - t224 + (t224 - t336) * t262) * t349 - t330 * t376) * t230) * t277, (t221 * t377 - t233 * t367) * t355 + (-t221 * t376 + (-t256 * t319 - t283 * t357) * t233 + (t221 * t352 - t234 * t367) * t219 - (t297 * t357 - t216 * t272 - t229 * t243 - t292 * t319 + (-t229 * t289 - t280 * t319) * t224) * t349 - (t280 * t357 - t216 * t289 - t229 * t267 - t258 * t319 + (t229 * t272 + t297 * t319) * t224) * t350) * t230, 0, (t220 * t377 + t233 * t276) * t355 + (t220 * t219 * t352 - t240 * t233 + (t276 * t219 - t220 * t241 - (-t218 * t272 - t226 * t243 + t266 + (-t226 * t289 + t271) * t224) * t373 - (-t218 * t289 - t226 * t267 + t242 + (t226 * t272 - t288) * t224) * t374) * t234) * t230, 0, 0; (t246 * t340 - t250 * t375) * t354 + ((t250 * qJD(6) + t242 * t318 + t258 * t314) * t246 + t250 * t345 + (t340 * t228 + (t340 * qJD(6) - t242 * t314 + t258 * t318) * t339 - t250 * t227) * t247) * t237, (-t246 * t253 - t254 * t375) * t354 + (t254 * t345 + t342 * t246 * t314 - t329 * t361 + (t318 * t339 * t342 - t254 * t227 - t253 * t228 + t329 * t364) * t247) * t237, 0, t335 * t277 * t354 + (-t335 * t241 + ((qJD(6) * t246 + t345) * t314 + (-t227 * t314 + (t228 + t387) * t318) * t247) * t277) * t237, 0, -0.2e1 * t380 - 0.2e1 * (t227 * t247 * t237 - (-t237 * t378 - t247 * t380) * t339) * t339;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end