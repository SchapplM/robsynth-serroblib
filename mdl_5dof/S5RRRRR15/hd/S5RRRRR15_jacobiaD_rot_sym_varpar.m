% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRRR15_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRR15_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR15_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(5));
	t93 = t99 ^ 2;
	t100 = cos(pkin(5));
	t95 = 0.1e1 / t100 ^ 2;
	t104 = cos(qJ(1));
	t98 = t104 ^ 2;
	t89 = t93 * t95 * t98 + 0.1e1;
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
	t118 = t103 * t104;
	t101 = sin(qJ(2));
	t120 = t102 * t101;
	t113 = t100 * t120 - t118;
	t77 = 0.1e1 / t113;
	t94 = 0.1e1 / t100;
	t68 = 0.1e1 / t72 ^ 2;
	t78 = 0.1e1 / t113 ^ 2;
	t119 = t102 * t103;
	t121 = t101 * t104;
	t81 = -t100 * t121 - t119;
	t82 = t100 * t119 + t121;
	t71 = t81 * qJD(1) - t82 * qJD(2);
	t128 = t71 * t77 * t78;
	t115 = t100 * t118;
	t70 = -qJD(1) * t115 - qJD(2) * t118 + (qJD(2) * t100 + qJD(1)) * t120;
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
	t80 = t115 - t120;
	t63 = (-t104 * t85 * t86 * t125 + t84 * t116) * t102;
	t92 = t99 * t93;
	t73 = 0.1e1 / t75;
	t69 = t67 * t68;
	t66 = t68 * t93 * t97 + 0.1e1;
	t62 = qJD(1) * t63;
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:08
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (1384->49), mult. (1182->120), div. (158->12), fcn. (1387->13), ass. (0->69)
	t156 = sin(qJ(1));
	t200 = t156 / 0.2e1;
	t155 = cos(pkin(5));
	t148 = 0.1e1 / t155 ^ 2;
	t154 = sin(pkin(5));
	t146 = t154 ^ 2;
	t157 = cos(qJ(1));
	t152 = t157 ^ 2;
	t139 = t152 * t146 * t148 + 0.1e1;
	t151 = t156 ^ 2;
	t186 = 0.1e1 / t139 ^ 2 * t151;
	t199 = t148 * t186;
	t153 = qJ(2) + qJ(3);
	t141 = pkin(5) + t153;
	t142 = pkin(5) - t153;
	t169 = cos(t142) + cos(t141);
	t143 = sin(t153);
	t181 = t157 * t143;
	t126 = t169 * t200 + t181;
	t119 = t126 ^ 2;
	t144 = cos(t153);
	t140 = t157 * t144;
	t168 = sin(t142) - sin(t141);
	t170 = t168 * t200 + t140;
	t121 = 0.1e1 / t170 ^ 2;
	t111 = t119 * t121 + 0.1e1;
	t109 = 0.1e1 / t111;
	t195 = -t157 / 0.2e1;
	t125 = -t156 * t144 - t168 * t195;
	t150 = qJD(2) + qJD(3);
	t129 = t169 * t150;
	t196 = -t156 / 0.2e1;
	t108 = t125 * qJD(1) + t129 * t196 - t150 * t181;
	t120 = 0.1e1 / t170;
	t187 = t120 * t121 * t108;
	t166 = -qJD(1) * t169 / 0.2e1;
	t167 = t168 * t150;
	t183 = t156 * t143;
	t107 = qJD(1) * t183 - t150 * t140 + t157 * t166 + t167 * t196;
	t188 = t121 * t126;
	t174 = t107 * t188;
	t192 = (-t119 * t187 - t174) / t111 ^ 2;
	t198 = t109 * t187 + t121 * t192;
	t190 = t109 * t121;
	t197 = -t108 * t190 - 0.2e1 * t120 * t192;
	t180 = t157 * t154;
	t138 = atan2(t180, t155);
	t133 = sin(t138);
	t134 = cos(t138);
	t118 = t133 * t180 + t134 * t155;
	t115 = 0.1e1 / t118;
	t147 = 0.1e1 / t155;
	t116 = 0.1e1 / t118 ^ 2;
	t194 = t157 / 0.2e1;
	t191 = t109 * t120;
	t189 = t116 * t156;
	t185 = t146 * t147;
	t182 = t156 * t150;
	t179 = qJD(1) * t157;
	t136 = 0.1e1 / t139;
	t173 = (t136 - 0.1e1) * t154;
	t171 = -0.2e1 * t147 * t199;
	t106 = (-t134 * t136 * t157 * t185 + t133 * t173) * t156;
	t145 = t154 * t146;
	t117 = t115 * t116;
	t114 = t151 * t146 * t116 + 0.1e1;
	t105 = qJD(1) * t106;
	t102 = t108 * t191 - 0.2e1 * t109 * t174 - 0.2e1 * t198 * t119 + t197 * t170;
	t1 = [(-t136 * t147 * t154 + t145 * t171) * t179, 0, 0, 0, 0; (0.2e1 * (t106 * t189 - t115 * t157) / t114 ^ 2 * (-t105 * t117 * t151 + t179 * t189) * t146 + ((0.2e1 * t106 * t117 * t156 - t116 * t157) * t105 + (-t156 * t115 + ((-t106 + (-t145 * t199 - t173) * t156 * t133) * t157 - (t152 * t146 ^ 2 * t171 + (-t186 + (0.2e1 * t151 - t152) * t136) * t185) * t156 * t134) * t116) * qJD(1)) / t114) * t154, 0, 0, 0, 0; (-t143 * t179 - t144 * t182 + t156 * t166 + t167 * t194) * t191 - (-t170 * qJD(1) + t129 * t195 + t143 * t182) * t109 * t188 + t197 * (t169 * t194 - t183) + (t107 * t190 + 0.2e1 * t198 * t126) * t125, t102, t102, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:08
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (2588->51), mult. (1659->122), div. (183->12), fcn. (1662->13), ass. (0->65)
	t156 = qJ(2) + qJ(3) + qJ(4);
	t177 = pkin(5) + t156;
	t178 = pkin(5) - t156;
	t140 = cos(t178) / 0.2e1 + cos(t177) / 0.2e1;
	t166 = sin(qJ(1));
	t153 = sin(t156);
	t167 = cos(qJ(1));
	t187 = t167 * t153;
	t134 = t166 * t140 + t187;
	t128 = t134 ^ 2;
	t139 = sin(t177) / 0.2e1 - sin(t178) / 0.2e1;
	t154 = cos(t156);
	t148 = t167 * t154;
	t179 = -t166 * t139 + t148;
	t130 = 0.1e1 / t179 ^ 2;
	t199 = t128 * t130;
	t165 = cos(pkin(5));
	t160 = 0.1e1 / t165 ^ 2;
	t164 = sin(pkin(5));
	t158 = t164 ^ 2;
	t163 = t167 ^ 2;
	t147 = t163 * t158 * t160 + 0.1e1;
	t162 = t166 ^ 2;
	t192 = 0.1e1 / t147 ^ 2 * t162;
	t198 = t160 * t192;
	t186 = t167 * t164;
	t146 = atan2(t186, t165);
	t142 = sin(t146);
	t143 = cos(t146);
	t126 = t142 * t186 + t143 * t165;
	t123 = 0.1e1 / t126;
	t129 = 0.1e1 / t179;
	t159 = 0.1e1 / t165;
	t124 = 0.1e1 / t126 ^ 2;
	t119 = 0.1e1 + t199;
	t155 = qJD(2) + qJD(3) + qJD(4);
	t176 = t139 * t155;
	t184 = qJD(1) * t167;
	t185 = qJD(1) * t166;
	t115 = -t140 * t184 - t155 * t148 + t153 * t185 + t166 * t176;
	t193 = t130 * t134;
	t182 = t115 * t193;
	t188 = t166 * t154;
	t133 = -t167 * t139 - t188;
	t137 = t140 * t155;
	t116 = t133 * qJD(1) - t166 * t137 - t155 * t187;
	t131 = t129 * t130;
	t196 = t116 * t131;
	t197 = (-t128 * t196 - t182) / t119 ^ 2;
	t117 = 0.1e1 / t119;
	t195 = t117 * t130;
	t194 = t124 * t166;
	t191 = t158 * t159;
	t189 = t166 * t153;
	t183 = 0.2e1 * t133 * t134;
	t144 = 0.1e1 / t147;
	t181 = (t144 - 0.1e1) * t164;
	t180 = -0.2e1 * t159 * t198;
	t114 = (-t143 * t144 * t167 * t191 + t142 * t181) * t166;
	t157 = t164 * t158;
	t125 = t123 * t124;
	t122 = t162 * t158 * t124 + 0.1e1;
	t113 = qJD(1) * t114;
	t110 = 0.2e1 * (-t129 * t179 - t199) * t197 + (-0.2e1 * t182 + (-0.2e1 * t128 * t131 - t130 * t179 + t129) * t116) * t117;
	t1 = [(-t144 * t159 * t164 + t157 * t180) * t184, 0, 0, 0, 0; (0.2e1 * (t114 * t194 - t123 * t167) / t122 ^ 2 * (-t113 * t125 * t162 + t184 * t194) * t158 + ((0.2e1 * t114 * t125 * t166 - t124 * t167) * t113 + (-t166 * t123 + ((-t114 + (-t157 * t198 - t181) * t166 * t142) * t167 - (t163 * t158 ^ 2 * t180 + (-t192 + (0.2e1 * t162 - t163) * t144) * t191) * t166 * t143) * t124) * qJD(1)) / t122) * t164, 0, 0, 0, 0; t133 * t115 * t195 + t130 * t183 * t197 + (-t116 * t195 - 0.2e1 * t129 * t197) * (t167 * t140 - t189) + ((-t140 * t185 - t153 * t184 - t155 * t188 - t167 * t176) * t129 - (-t179 * qJD(1) - t167 * t137 + t155 * t189) * t193 + t183 * t196) * t117, t110, t110, t110, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:08
	% EndTime: 2024-09-27 22:28:09
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (14119->121), mult. (17828->257), div. (963->12), fcn. (22205->13), ass. (0->127)
	t304 = qJ(2) + qJ(3) + qJ(4);
	t301 = sin(t304);
	t302 = cos(t304);
	t311 = cos(qJ(5));
	t307 = cos(pkin(6));
	t308 = cos(pkin(5));
	t310 = sin(qJ(1));
	t368 = t311 * t310;
	t348 = t308 * t368;
	t309 = sin(qJ(5));
	t312 = cos(qJ(1));
	t366 = t312 * t309;
	t324 = t307 * t348 + t366;
	t369 = t310 * t309;
	t349 = t308 * t369;
	t365 = t312 * t311;
	t325 = t307 * t365 - t349;
	t305 = sin(pkin(6));
	t306 = sin(pkin(5));
	t371 = t306 * t310;
	t351 = t305 * t371;
	t259 = -t325 * t301 - t324 * t302 + t311 * t351;
	t346 = t308 * t365;
	t290 = t307 * t369 - t346;
	t352 = t308 * t366;
	t323 = t307 * t352 + t368;
	t372 = t305 * t306;
	t347 = t312 * t372;
	t403 = t290 * t301 - t323 * t302 + t309 * t347;
	t287 = t307 * t349 - t365;
	t292 = t307 * t366 + t348;
	t326 = -t287 * t302 - t292 * t301 + t309 * t351;
	t254 = 0.1e1 / t326 ^ 2;
	t256 = t259 ^ 2;
	t402 = t256 * t254;
	t291 = t307 * t368 + t352;
	t303 = qJD(2) + qJD(3) + qJD(4);
	t401 = -t291 * qJD(1) - t292 * qJD(5) - t324 * t303;
	t335 = qJD(5) * t351;
	t336 = qJD(1) * t347;
	t341 = t290 * qJD(1) - t325 * qJD(5) + t287 * t303;
	t342 = t323 * qJD(1) + t324 * qJD(5) + t292 * t303;
	t242 = t341 * t301 - t342 * t302 + t309 * t336 + t311 * t335;
	t253 = 0.1e1 / t326;
	t249 = 0.1e1 + t402;
	t289 = -t307 * t346 + t369;
	t343 = t289 * qJD(1) + t287 * qJD(5) - t303 * t325;
	t241 = -t401 * t301 + t343 * t302 - t309 * t335 + t311 * t336;
	t387 = t254 * t259;
	t357 = t241 * t387;
	t255 = t253 * t254;
	t386 = t255 * t242;
	t393 = (-t256 * t386 + t357) / t249 ^ 2;
	t362 = -0.2e1 * t393;
	t247 = 0.1e1 / t249;
	t389 = t247 * t254;
	t399 = -t242 * t389 + t253 * t362;
	t398 = t289 * t302 + t291 * t301 + t311 * t347;
	t397 = -0.2e1 * t259 * t247 * t386 + t241 * t389 + t362 * t387;
	t327 = t302 * t305 * t308 + t306 * t307;
	t376 = t301 * t310;
	t354 = t305 * t376;
	t278 = t327 * t312 - t354;
	t286 = -t302 * t372 + t308 * t307;
	t268 = atan2(t278, t286);
	t263 = sin(t268);
	t264 = cos(t268);
	t246 = t263 * t278 + t264 * t286;
	t243 = 0.1e1 / t246;
	t282 = 0.1e1 / t286;
	t244 = 0.1e1 / t246 ^ 2;
	t283 = 0.1e1 / t286 ^ 2;
	t370 = t310 * t302;
	t375 = t301 * t312;
	t330 = t308 * t370 + t375;
	t350 = t307 * t371;
	t277 = t330 * t305 + t350;
	t275 = t277 ^ 2;
	t239 = t275 * t244 + 0.1e1;
	t353 = t308 * t376;
	t373 = t303 * t305;
	t250 = -t353 * t373 - qJD(1) * t354 + (t327 * qJD(1) + t302 * t373) * t312;
	t388 = t250 * t244;
	t276 = t278 ^ 2;
	t267 = t276 * t283 + 0.1e1;
	t265 = 0.1e1 / t267;
	t329 = t308 * t375 + t370;
	t251 = qJD(1) * t350 + (t330 * qJD(1) + t329 * t303) * t305;
	t355 = t303 * t372;
	t340 = t301 * t355;
	t332 = t283 * t340;
	t319 = -t251 * t282 - t278 * t332;
	t235 = t319 * t265;
	t331 = -t263 * t286 + t264 * t278;
	t231 = t331 * t235 - t263 * t251 + t264 * t340;
	t394 = t231 * t243 * t244;
	t395 = (-t275 * t394 + t277 * t388) / t239 ^ 2;
	t284 = t282 * t283;
	t382 = t278 * t283;
	t392 = (-t276 * t284 * t340 - t251 * t382) / t267 ^ 2;
	t391 = t244 * t277;
	t390 = t247 * t253;
	t385 = t263 * t277;
	t384 = t264 * t277;
	t383 = t278 * t282;
	t374 = t302 * t303;
	t367 = t312 * t302;
	t364 = 0.2e1 * t395;
	t363 = 0.2e1 * t394;
	t361 = 0.2e1 * t392;
	t359 = t247 * t387;
	t356 = t301 * t372;
	t344 = t282 * t361;
	t328 = t353 - t367;
	t322 = -t263 + (-t264 * t383 + t263) * t265;
	t281 = t329 * t305;
	t321 = t281 * t282 + t356 * t382;
	t280 = t328 * t305;
	t252 = ((-t308 * t367 + t376) * t303 + t328 * qJD(1)) * t305;
	t237 = 0.1e1 / t239;
	t236 = t321 * t265;
	t234 = t322 * t277;
	t232 = -t331 * t236 - t263 * t281 + t264 * t356;
	t230 = t321 * t361 + (t252 * t282 + (0.2e1 * t278 * t284 * t301 ^ 2 * t355 + (-t278 * t374 + (t281 * t303 + t251) * t301) * t283) * t372) * t265;
	t228 = (t343 * t301 + t401 * t302) * t390 + (t342 * t301 + t341 * t302) * t359 + t399 * (-t324 * t301 + t302 * t325) + t397 * (t287 * t301 - t292 * t302);
	t227 = (t232 * t391 + t243 * t280) * t364 + (t232 * t277 * t363 + (t280 * t231 - t232 * t250 - (t235 * t236 * t278 - t230 * t286 + t252) * t385 - (t230 * t278 + t236 * t251 + (t236 * t286 - t281) * t235) * t384) * t244 + ((-t329 * qJD(1) - t330 * t303) * t243 - (t264 * t374 + (t236 * t303 - t235) * t301 * t263) * t306 * t391) * t305) * t237;
	t1 = [t277 * t344 + (-t250 * t282 + t277 * t332) * t265, t230, t230, t230, 0; -0.2e1 * t278 * t243 * t395 + (-t251 * t243 + (-t231 * t278 - t234 * t250) * t244) * t237 + ((t234 * t363 - t322 * t388) * t237 + (t234 * t364 + (-(t235 * t265 * t383 - 0.2e1 * t392) * t385 - (t278 * t344 - t235 + (t235 - t319) * t265) * t384) * t237) * t244) * t277, t227, t227, t227, 0; ((t289 * t301 - t291 * t302) * t303 + t403 * qJD(5) + t259 * qJD(1)) * t390 + ((t290 * t302 + t301 * t323) * t303 + t398 * qJD(5) - t326 * qJD(1)) * t359 - t399 * t398 + t397 * t403, t228, t228, t228, 0.2e1 * (-t253 * t326 - t402) * t393 + (0.2e1 * t357 + (-t254 * t326 - 0.2e1 * t256 * t255 + t253) * t242) * t247;];
	JaD_rot = t1;
end