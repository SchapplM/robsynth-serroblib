% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR7
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
%   Wie in S6RRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:25
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
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:26
	% DurationCPUTime: 1.77s
	% Computational Cost: add. (7491->94), mult. (9702->187), div. (711->12), fcn. (12089->11), ass. (0->95)
	t274 = qJ(4) + pkin(10);
	t263 = sin(t274);
	t264 = cos(t274);
	t293 = sin(qJ(2));
	t295 = cos(qJ(2));
	t211 = -t295 * t263 + t293 * t264;
	t294 = sin(qJ(1));
	t202 = t211 * t294;
	t210 = t293 * t263 + t295 * t264;
	t189 = atan2(t202, t210);
	t184 = sin(t189);
	t185 = cos(t189);
	t204 = t210 * t294;
	t258 = -t184 * t210 + t185 * t202;
	t200 = t202 ^ 2;
	t208 = 0.1e1 / t210 ^ 2;
	t188 = t200 * t208 + 0.1e1;
	t186 = 0.1e1 / t188;
	t207 = 0.1e1 / t210;
	t276 = t202 * t211;
	t250 = -t204 * t207 - t208 * t276;
	t303 = t250 * t186;
	t161 = -t184 * t204 + t185 * t211 + t258 * t303;
	t296 = cos(qJ(1));
	t206 = t211 * t296;
	t322 = t161 * t206;
	t174 = t184 * t202 + t185 * t210;
	t171 = 0.1e1 / t174;
	t172 = 0.1e1 / t174 ^ 2;
	t205 = t210 * t296;
	t201 = t206 ^ 2;
	t170 = t172 * t201 + 0.1e1;
	t311 = qJD(2) - qJD(4);
	t178 = -t202 * qJD(1) + t311 * t205;
	t284 = t178 * t172;
	t179 = -t206 * qJD(1) - t311 * t204;
	t191 = t311 * t211;
	t277 = t202 * t208;
	t252 = -t179 * t207 + t191 * t277;
	t164 = t252 * t186;
	t159 = t258 * t164 - t179 * t184 - t185 * t191;
	t292 = t159 * t171 * t172;
	t273 = 0.2e1 * (-t201 * t292 + t206 * t284) / t170 ^ 2;
	t321 = (t171 * t205 + t172 * t322) * t273;
	t177 = t204 * qJD(1) + t311 * t206;
	t229 = sin(qJ(6));
	t230 = cos(qJ(6));
	t199 = t205 * t230 - t294 * t229;
	t266 = qJD(1) * t296;
	t175 = t199 * qJD(6) - t177 * t229 + t230 * t266;
	t248 = -t205 * t229 - t294 * t230;
	t312 = t248 * qJD(6);
	t176 = -t177 * t230 - t229 * t266 + t312;
	t192 = t248 ^ 2;
	t194 = 0.1e1 / t199 ^ 2;
	t183 = t192 * t194 + 0.1e1;
	t181 = 0.1e1 / t183;
	t193 = 0.1e1 / t199;
	t279 = t194 * t248;
	t251 = -t229 * t193 - t230 * t279;
	t286 = t176 * t193 * t194;
	t298 = -0.2e1 * t248;
	t262 = t286 * t298;
	t320 = (t251 * t178 - (((-t176 - t312) * t229 - t175 * t230) * t194 + (qJD(6) * t193 + t262) * t230) * t206) * t181;
	t319 = -t205 * t159 + t161 * t178;
	t272 = 0.2e1 * t292;
	t318 = -t177 * t171 - t272 * t322;
	t190 = t311 * t210;
	t316 = (t210 * t303 + t204) * t164 + t303 * t179 - t190;
	t180 = t205 * qJD(1) - t311 * t202;
	t315 = -(-t202 * t303 - t211) * t164 - t303 * t191 + t180;
	t305 = t191 * t208;
	t280 = t207 * t305;
	t313 = ((t179 * t211 - t190 * t202 - t191 * t204) * t208 - t180 * t207 - 0.2e1 * t276 * t280) * t186;
	t278 = t202 * t207;
	t301 = (-t185 * t278 + t184) * t186 - t184;
	t297 = -0.2e1 * t206;
	t291 = (-t175 * t279 - t192 * t286) / t183 ^ 2;
	t290 = (-t179 * t277 + t200 * t280) / t188 ^ 2;
	t283 = t181 * t194;
	t282 = t184 * t206;
	t281 = t185 * t206;
	t271 = -0.2e1 * t291;
	t270 = -0.2e1 * t290;
	t269 = t194 * t291;
	t268 = t207 * t290;
	t267 = t175 * t283;
	t265 = qJD(1) * t294;
	t197 = -t204 * t230 - t296 * t229;
	t249 = t204 * t229 - t296 * t230;
	t168 = 0.1e1 / t170;
	t163 = t301 * t206;
	t157 = t250 * t270 + t313;
	t156 = 0.2e1 * t250 * t290 - t313;
	t1 = [t268 * t297 + (t178 * t207 + t206 * t305) * t186, t156, 0, t157, 0, 0; -t202 * t171 * t273 + (-t179 * t171 + (-t159 * t202 - t163 * t178) * t172) * t168 - ((-t163 * t272 + t301 * t284) * t168 + (-t163 * t273 + ((t164 * t186 * t278 + t270) * t282 + (0.2e1 * t202 * t268 - t164 + (t164 - t252) * t186) * t281) * t168) * t172) * t206, t321 + (((t156 * t202 + t316) * t281 + (-t156 * t210 + t315) * t282 - t319) * t172 - t318) * t168, 0, -t321 + (((t157 * t202 - t316) * t281 + (-t157 * t210 - t315) * t282 + t319) * t172 + t318) * t168, 0, 0; (t269 * t298 - t267) * t197 - (-t176 * t283 + t193 * t271) * t249 + ((t197 * qJD(6) - t180 * t229 - t230 * t265) * t193 + (t249 * qJD(6) - t180 * t230 + t229 * t265) * t279 + t197 * t262) * t181, t251 * t291 * t297 + t320, 0, -t251 * t206 * t271 - t320, 0, t271 + (t267 - (-t181 * t286 - t269) * t248) * t298;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end