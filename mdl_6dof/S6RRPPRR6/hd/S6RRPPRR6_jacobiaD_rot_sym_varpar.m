% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR6
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
%   Wie in S6RRPPRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:55
	% DurationCPUTime: 0.76s
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
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:56
	% DurationCPUTime: 1.78s
	% Computational Cost: add. (7491->94), mult. (9702->187), div. (711->12), fcn. (12089->11), ass. (0->95)
	t272 = pkin(10) + qJ(5);
	t261 = sin(t272);
	t262 = cos(t272);
	t291 = sin(qJ(2));
	t293 = cos(qJ(2));
	t209 = -t293 * t261 + t291 * t262;
	t292 = sin(qJ(1));
	t200 = t209 * t292;
	t208 = t291 * t261 + t293 * t262;
	t187 = atan2(t200, t208);
	t182 = sin(t187);
	t183 = cos(t187);
	t202 = t208 * t292;
	t256 = -t182 * t208 + t183 * t200;
	t198 = t200 ^ 2;
	t206 = 0.1e1 / t208 ^ 2;
	t186 = t198 * t206 + 0.1e1;
	t184 = 0.1e1 / t186;
	t205 = 0.1e1 / t208;
	t274 = t200 * t209;
	t248 = -t202 * t205 - t206 * t274;
	t301 = t248 * t184;
	t159 = -t182 * t202 + t183 * t209 + t256 * t301;
	t294 = cos(qJ(1));
	t204 = t209 * t294;
	t320 = t159 * t204;
	t172 = t182 * t200 + t183 * t208;
	t169 = 0.1e1 / t172;
	t170 = 0.1e1 / t172 ^ 2;
	t203 = t208 * t294;
	t199 = t204 ^ 2;
	t168 = t199 * t170 + 0.1e1;
	t309 = qJD(2) - qJD(5);
	t176 = -t200 * qJD(1) + t309 * t203;
	t282 = t176 * t170;
	t177 = -t204 * qJD(1) - t309 * t202;
	t189 = t309 * t209;
	t275 = t200 * t206;
	t250 = -t177 * t205 + t189 * t275;
	t162 = t250 * t184;
	t157 = t256 * t162 - t182 * t177 - t183 * t189;
	t290 = t157 * t169 * t170;
	t271 = 0.2e1 * (-t199 * t290 + t204 * t282) / t168 ^ 2;
	t319 = (t169 * t203 + t170 * t320) * t271;
	t175 = qJD(1) * t202 + t309 * t204;
	t227 = sin(qJ(6));
	t228 = cos(qJ(6));
	t197 = t203 * t228 - t292 * t227;
	t264 = qJD(1) * t294;
	t173 = t197 * qJD(6) - t175 * t227 + t228 * t264;
	t246 = -t203 * t227 - t292 * t228;
	t310 = t246 * qJD(6);
	t174 = -t175 * t228 - t227 * t264 + t310;
	t190 = t246 ^ 2;
	t192 = 0.1e1 / t197 ^ 2;
	t181 = t190 * t192 + 0.1e1;
	t179 = 0.1e1 / t181;
	t191 = 0.1e1 / t197;
	t277 = t192 * t246;
	t249 = -t227 * t191 - t228 * t277;
	t284 = t174 * t191 * t192;
	t296 = -0.2e1 * t246;
	t260 = t284 * t296;
	t318 = (t249 * t176 - (((-t174 - t310) * t227 - t173 * t228) * t192 + (qJD(6) * t191 + t260) * t228) * t204) * t179;
	t317 = -t203 * t157 + t159 * t176;
	t270 = 0.2e1 * t290;
	t316 = -t175 * t169 - t270 * t320;
	t188 = t309 * t208;
	t314 = (t208 * t301 + t202) * t162 + t301 * t177 - t188;
	t178 = qJD(1) * t203 - t309 * t200;
	t313 = -(-t200 * t301 - t209) * t162 - t301 * t189 + t178;
	t303 = t189 * t206;
	t278 = t205 * t303;
	t311 = ((t177 * t209 - t188 * t200 - t189 * t202) * t206 - t178 * t205 - 0.2e1 * t274 * t278) * t184;
	t276 = t200 * t205;
	t299 = (-t183 * t276 + t182) * t184 - t182;
	t295 = -0.2e1 * t204;
	t289 = (-t173 * t277 - t190 * t284) / t181 ^ 2;
	t288 = (-t177 * t275 + t198 * t278) / t186 ^ 2;
	t281 = t179 * t192;
	t280 = t182 * t204;
	t279 = t183 * t204;
	t269 = -0.2e1 * t289;
	t268 = -0.2e1 * t288;
	t267 = t192 * t289;
	t266 = t205 * t288;
	t265 = t173 * t281;
	t263 = qJD(1) * t292;
	t195 = -t202 * t228 - t294 * t227;
	t247 = t202 * t227 - t294 * t228;
	t166 = 0.1e1 / t168;
	t161 = t299 * t204;
	t155 = t248 * t268 + t311;
	t154 = 0.2e1 * t248 * t288 - t311;
	t1 = [t266 * t295 + (t176 * t205 + t204 * t303) * t184, t154, 0, 0, t155, 0; -t200 * t169 * t271 + (-t177 * t169 + (-t157 * t200 - t161 * t176) * t170) * t166 - ((-t161 * t270 + t299 * t282) * t166 + (-t161 * t271 + ((t162 * t184 * t276 + t268) * t280 + (0.2e1 * t200 * t266 - t162 + (t162 - t250) * t184) * t279) * t166) * t170) * t204, t319 + (((t154 * t200 + t314) * t279 + (-t154 * t208 + t313) * t280 - t317) * t170 - t316) * t166, 0, 0, -t319 + (((t155 * t200 - t314) * t279 + (-t155 * t208 - t313) * t280 + t317) * t170 + t316) * t166, 0; (t267 * t296 - t265) * t195 - (-t174 * t281 + t191 * t269) * t247 + ((t195 * qJD(6) - t178 * t227 - t228 * t263) * t191 + (t247 * qJD(6) - t178 * t228 + t227 * t263) * t277 + t195 * t260) * t179, t249 * t289 * t295 + t318, 0, 0, -t249 * t204 * t269 - t318, t269 + (t265 - (-t179 * t284 - t267) * t246) * t296;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end