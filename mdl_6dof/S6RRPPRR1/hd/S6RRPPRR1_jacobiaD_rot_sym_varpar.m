% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR1
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
%   Wie in S6RRPPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:46
	% DurationCPUTime: 0.81s
	% Computational Cost: add. (1893->71), mult. (1835->159), div. (470->13), fcn. (2177->7), ass. (0->71)
	t92 = sin(qJ(1));
	t122 = qJD(1) * t92;
	t141 = 0.2e1 * t92;
	t119 = qJD(2) * t92;
	t93 = cos(qJ(1));
	t121 = qJD(1) * t93;
	t85 = qJ(2) + pkin(10);
	t83 = sin(t85);
	t112 = t83 * t121;
	t79 = t83 ^ 2;
	t84 = cos(t85);
	t81 = 0.1e1 / t84 ^ 2;
	t127 = t79 * t81;
	t87 = t92 ^ 2;
	t74 = t87 * t127 + 0.1e1;
	t72 = 0.1e1 / t74;
	t80 = 0.1e1 / t84;
	t58 = (-(-t84 * t119 - t112) * t80 + t119 * t127) * t72;
	t139 = t58 - t119;
	t88 = t93 ^ 2;
	t89 = 0.1e1 / t93;
	t138 = (t87 / t88 + 0.1e1) * t89 * t122;
	t123 = t92 * t83;
	t71 = atan2(-t123, -t84);
	t69 = sin(t71);
	t114 = t69 * t123;
	t70 = cos(t71);
	t65 = -t70 * t84 - t114;
	t62 = 0.1e1 / t65;
	t63 = 0.1e1 / t65 ^ 2;
	t137 = -0.2e1 * t83;
	t136 = t72 - 0.1e1;
	t129 = t70 * t83;
	t54 = (-t58 * t92 + qJD(2)) * t129 + (t139 * t84 - t112) * t69;
	t135 = t54 * t62 * t63;
	t134 = t58 * t83;
	t133 = t63 * t83;
	t132 = t63 * t93;
	t125 = t80 * t83;
	t78 = t83 * t79;
	t82 = t80 * t81;
	t101 = qJD(2) * (t78 * t82 + t125);
	t105 = t79 * t92 * t121;
	t131 = (t87 * t101 + t81 * t105) / t74 ^ 2;
	t130 = t69 * t92;
	t128 = t79 * t80;
	t126 = t79 * t88;
	t124 = t87 / t93 ^ 2;
	t120 = qJD(2) * t84;
	t61 = t63 * t126 + 0.1e1;
	t118 = 0.2e1 * (-t126 * t135 + (t83 * t88 * t120 - t105) * t63) / t61 ^ 2;
	t117 = 0.2e1 * t135;
	t77 = t81 * t124 + 0.1e1;
	t116 = 0.2e1 * (t82 * qJD(2) * t83 * t124 + t81 * t138) / t77 ^ 2;
	t115 = t83 * t132;
	t113 = t72 * t128;
	t111 = 0.1e1 + t127;
	t110 = 0.1e1 + t124;
	t109 = t83 * t118;
	t108 = t131 * t137;
	t107 = t131 * t141;
	t106 = t92 * t113;
	t104 = t111 * t93;
	t102 = t110 * t83 * t81;
	t75 = 0.1e1 / t77;
	t67 = t111 * t92 * t72;
	t59 = 0.1e1 / t61;
	t57 = (t136 * t83 * t69 - t70 * t106) * t93;
	t56 = -t84 * t130 + t129 + (-t70 * t123 + t69 * t84) * t67;
	t55 = -t111 * t107 + (qJD(1) * t104 + t101 * t141) * t72;
	t1 = [t93 * t80 * t108 + (qJD(2) * t104 - t122 * t125) * t72, t55, 0, 0, 0, 0; (t62 * t109 + (-t62 * t120 + (qJD(1) * t57 + t54) * t133) * t59) * t92 + (t63 * t109 * t57 + (-((t58 * t106 + t136 * t120 + t108) * t69 + (t107 * t128 - t134 + (t134 + (-t78 * t81 + t137) * t119) * t72) * t70) * t115 + (t83 * t117 - t63 * t120) * t57 + (-t62 + ((-t87 + t88) * t70 * t113 + t136 * t114) * t63) * t83 * qJD(1)) * t59) * t93, (t56 * t133 - t62 * t84) * t93 * t118 + ((-t62 * t122 + (-qJD(2) * t56 - t54) * t132) * t84 + (-t93 * qJD(2) * t62 - (-t55 * t70 * t92 - t139 * t69 + (-qJD(2) * t69 - t121 * t70 + t130 * t58) * t67) * t115 + (t93 * t117 + t63 * t122) * t56 - ((t55 - t121) * t69 + ((-t67 * t92 + 0.1e1) * qJD(2) + (t67 - t92) * t58) * t70) * t84 * t132) * t83) * t59, 0, 0, 0, 0; t110 * t80 * t116 + (-qJD(2) * t102 - 0.2e1 * t80 * t138) * t75, t89 * t81 * t116 * t123 + ((-0.2e1 * t79 * t82 - t80) * t89 * t119 - qJD(1) * t102) * t75, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:47
	% DurationCPUTime: 1.81s
	% Computational Cost: add. (7491->94), mult. (9702->187), div. (711->12), fcn. (12089->11), ass. (0->95)
	t274 = qJ(2) + pkin(10);
	t263 = sin(t274);
	t264 = cos(t274);
	t293 = sin(qJ(5));
	t295 = cos(qJ(5));
	t211 = t263 * t295 - t264 * t293;
	t294 = sin(qJ(1));
	t202 = t211 * t294;
	t210 = t263 * t293 + t264 * t295;
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
	t170 = t201 * t172 + 0.1e1;
	t311 = qJD(2) - qJD(5);
	t178 = -t202 * qJD(1) + t311 * t205;
	t284 = t178 * t172;
	t179 = -t206 * qJD(1) - t311 * t204;
	t191 = t311 * t211;
	t277 = t202 * t208;
	t252 = -t179 * t207 + t191 * t277;
	t164 = t252 * t186;
	t159 = t258 * t164 - t184 * t179 - t185 * t191;
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
	t1 = [t268 * t297 + (t178 * t207 + t206 * t305) * t186, t156, 0, 0, t157, 0; -t202 * t171 * t273 + (-t179 * t171 + (-t159 * t202 - t163 * t178) * t172) * t168 - ((-t163 * t272 + t301 * t284) * t168 + (-t163 * t273 + ((t164 * t186 * t278 + t270) * t282 + (0.2e1 * t202 * t268 - t164 + (t164 - t252) * t186) * t281) * t168) * t172) * t206, t321 + (((t156 * t202 + t316) * t281 + (-t156 * t210 + t315) * t282 - t319) * t172 - t318) * t168, 0, 0, -t321 + (((t157 * t202 - t316) * t281 + (-t157 * t210 - t315) * t282 + t319) * t172 + t318) * t168, 0; (t269 * t298 - t267) * t197 - (-t176 * t283 + t193 * t271) * t249 + ((t197 * qJD(6) - t180 * t229 - t230 * t265) * t193 + (t249 * qJD(6) - t180 * t230 + t229 * t265) * t279 + t197 * t262) * t181, t251 * t291 * t297 + t320, 0, 0, -t251 * t206 * t271 - t320, t271 + (t267 - (-t181 * t286 - t269) * t248) * t298;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end