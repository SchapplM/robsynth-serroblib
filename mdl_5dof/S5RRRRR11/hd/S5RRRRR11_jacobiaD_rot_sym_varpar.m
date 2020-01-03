% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR11
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
%   Wie in S5RRRRR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRR11_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR11_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:45:34
	% EndTime: 2019-12-31 22:45:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:45:35
	% EndTime: 2019-12-31 22:45:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:45:35
	% EndTime: 2019-12-31 22:45:35
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(5));
	t93 = t99 ^ 2;
	t100 = cos(pkin(5));
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
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:45:35
	% EndTime: 2019-12-31 22:45:35
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (1479->91), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->91)
	t171 = sin(qJ(2));
	t172 = sin(qJ(1));
	t225 = cos(pkin(5));
	t195 = t172 * t225;
	t193 = t171 * t195;
	t174 = cos(qJ(2));
	t175 = cos(qJ(1));
	t209 = t175 * t174;
	t157 = -t193 + t209;
	t170 = sin(qJ(3));
	t173 = cos(qJ(3));
	t169 = sin(pkin(5));
	t213 = t169 * t172;
	t185 = -t157 * t170 + t173 * t213;
	t227 = t185 * qJD(3);
	t194 = t175 * t225;
	t192 = t174 * t194;
	t210 = t172 * t171;
	t153 = -t192 + t210;
	t212 = t169 * t174;
	t147 = atan2(-t153, -t212);
	t145 = sin(t147);
	t146 = cos(t147);
	t151 = t153 ^ 2;
	t165 = 0.1e1 / t169 ^ 2;
	t167 = 0.1e1 / t174 ^ 2;
	t150 = t151 * t165 * t167 + 0.1e1;
	t148 = 0.1e1 / t150;
	t164 = 0.1e1 / t169;
	t166 = 0.1e1 / t174;
	t199 = t153 * t164 * t166;
	t226 = (t146 * t199 - t145) * t148 + t145;
	t129 = -t145 * t153 - t146 * t212;
	t126 = 0.1e1 / t129;
	t144 = t157 * t173 + t170 * t213;
	t138 = 0.1e1 / t144;
	t127 = 0.1e1 / t129 ^ 2;
	t139 = 0.1e1 / t144 ^ 2;
	t182 = -t171 * t194 - t172 * t174;
	t183 = -t175 * t171 - t174 * t195;
	t135 = -t183 * qJD(1) - t182 * qJD(2);
	t207 = qJD(2) * t171;
	t196 = t167 * t207;
	t184 = t135 * t166 + t153 * t196;
	t215 = t148 * t164;
	t118 = t184 * t215;
	t188 = t145 * t212 - t146 * t153;
	t200 = t146 * t169 * t171;
	t114 = qJD(2) * t200 + t188 * t118 - t145 * t135;
	t224 = t114 * t126 * t127;
	t214 = t167 * t171;
	t187 = t153 * t214 - t166 * t182;
	t119 = t187 * t215;
	t115 = t188 * t119 + t145 * t182 + t200;
	t223 = t115 * t183;
	t134 = t182 * qJD(1) + t183 * qJD(2);
	t208 = qJD(1) * t169;
	t197 = t175 * t208;
	t124 = t144 * qJD(3) + t134 * t170 - t173 * t197;
	t137 = t185 ^ 2;
	t132 = t137 * t139 + 0.1e1;
	t218 = t139 * t185;
	t125 = t134 * t173 + t170 * t197 + t227;
	t220 = t125 * t138 * t139;
	t222 = (-t124 * t218 - t137 * t220) / t132 ^ 2;
	t168 = t166 * t167;
	t221 = (t135 * t153 * t167 + t151 * t168 * t207) * t165 / t150 ^ 2;
	t191 = qJD(2) * t225 + qJD(1);
	t206 = qJD(2) * t174;
	t133 = -qJD(1) * t192 - t175 * t206 + t191 * t210;
	t219 = t133 * t127;
	t217 = t145 * t183;
	t216 = t146 * t183;
	t211 = t169 * t175;
	t152 = t183 ^ 2;
	t122 = t127 * t152 + 0.1e1;
	t205 = 0.2e1 * (-t152 * t224 + t183 * t219) / t122 ^ 2;
	t204 = 0.2e1 * t224;
	t203 = 0.2e1 * t222;
	t202 = -0.2e1 * t221;
	t201 = t185 * t220;
	t198 = t172 * t208;
	t189 = t170 * t138 + t173 * t218;
	t186 = -t170 * t182 + t173 * t211;
	t142 = t170 * t211 + t173 * t182;
	t136 = -qJD(1) * t193 - t172 * t207 + t191 * t209;
	t130 = 0.1e1 / t132;
	t120 = 0.1e1 / t122;
	t117 = t226 * t183;
	t113 = (t187 * t202 + (t135 * t214 + t136 * t166 + (-t182 * t214 + (0.2e1 * t168 * t171 ^ 2 + t166) * t153) * qJD(2)) * t148) * t164;
	t1 = [(-t183 * t166 * t202 + (-t133 * t166 - t183 * t196) * t148) * t164, t113, 0, 0, 0; t153 * t126 * t205 + (-t135 * t126 + (t114 * t153 + t117 * t133) * t127) * t120 - ((t117 * t204 - t226 * t219) * t120 + (t117 * t205 + ((t118 * t148 * t199 + t202) * t217 + (0.2e1 * t199 * t221 - t118 + (-t184 * t164 + t118) * t148) * t216) * t120) * t127) * t183, (-t126 * t157 - t127 * t223) * t205 + (-t204 * t223 + t134 * t126 + (-t157 * t114 + t115 * t133 + (t169 * t206 - t113 * t153 - t119 * t135 + (t119 * t212 + t182) * t118) * t216 + (t118 * t119 * t153 - t136 + (t113 * t174 + (-qJD(2) * t119 - t118) * t171) * t169) * t217) * t127) * t120, 0, 0, 0; (t138 * t186 - t142 * t218) * t203 + ((t142 * qJD(3) - t136 * t170 + t173 * t198) * t138 - 0.2e1 * t142 * t201 + (t186 * t125 + (t186 * qJD(3) - t136 * t173 - t170 * t198) * t185 - t142 * t124) * t139) * t130, -t189 * t183 * t203 + (t189 * t133 - ((-qJD(3) * t138 + 0.2e1 * t201) * t173 + (t124 * t173 + (t125 + t227) * t170) * t139) * t183) * t130, -0.2e1 * t222 - 0.2e1 * (t124 * t139 * t130 - (-t130 * t220 - t139 * t222) * t185) * t185, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:45:35
	% EndTime: 2019-12-31 22:45:37
	% DurationCPUTime: 1.41s
	% Computational Cost: add. (4522->148), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->126)
	t256 = cos(pkin(5));
	t259 = sin(qJ(2));
	t331 = sin(qJ(1));
	t293 = t331 * t259;
	t283 = t256 * t293;
	t287 = t331 * qJD(2);
	t262 = cos(qJ(2));
	t263 = cos(qJ(1));
	t309 = t263 * t262;
	t255 = sin(pkin(5));
	t313 = t255 * t263;
	t336 = -qJD(1) * t283 - t259 * t287 + (qJD(2) * t256 + qJD(1)) * t309 - qJD(3) * t313;
	t258 = sin(qJ(3));
	t261 = cos(qJ(3));
	t292 = t331 * t262;
	t310 = t263 * t259;
	t275 = -t256 * t310 - t292;
	t228 = -t258 * t275 + t261 * t313;
	t315 = t255 * t259;
	t239 = -t256 * t261 + t258 * t315;
	t217 = atan2(-t228, t239);
	t212 = sin(t217);
	t213 = cos(t217);
	t195 = -t212 * t228 + t213 * t239;
	t193 = 0.1e1 / t195 ^ 2;
	t244 = -t283 + t309;
	t294 = t255 * t331;
	t274 = -t244 * t258 + t261 * t294;
	t225 = t274 ^ 2;
	t191 = t225 * t193 + 0.1e1;
	t273 = -t256 * t292 - t310;
	t221 = t275 * qJD(1) + t273 * qJD(2);
	t234 = t244 * t261 + t258 * t294;
	t291 = qJD(1) * t313;
	t199 = t234 * qJD(3) + t221 * t258 - t261 * t291;
	t324 = t199 * t193;
	t224 = t228 ^ 2;
	t237 = 0.1e1 / t239 ^ 2;
	t216 = t224 * t237 + 0.1e1;
	t214 = 0.1e1 / t216;
	t288 = t331 * qJD(1);
	t282 = t255 * t288;
	t306 = qJD(3) * t261;
	t201 = t336 * t258 - t261 * t282 - t275 * t306;
	t240 = t256 * t258 + t261 * t315;
	t307 = qJD(2) * t262;
	t290 = t255 * t307;
	t226 = t240 * qJD(3) + t258 * t290;
	t236 = 0.1e1 / t239;
	t318 = t228 * t237;
	t279 = -t201 * t236 + t226 * t318;
	t183 = t279 * t214;
	t280 = -t212 * t239 - t213 * t228;
	t178 = t280 * t183 - t212 * t201 + t213 * t226;
	t192 = 0.1e1 / t195;
	t194 = t192 * t193;
	t329 = t178 * t194;
	t304 = 0.2e1 * (-t225 * t329 - t274 * t324) / t191 ^ 2;
	t335 = t226 * t237;
	t295 = t256 * t309;
	t241 = -t293 + t295;
	t314 = t255 * t262;
	t276 = -t236 * t241 + t314 * t318;
	t334 = t258 * t276;
	t202 = (qJD(3) * t275 + t282) * t258 + t336 * t261;
	t257 = sin(qJ(4));
	t260 = cos(qJ(4));
	t211 = t234 * t260 - t257 * t273;
	t205 = 0.1e1 / t211;
	t206 = 0.1e1 / t211 ^ 2;
	t333 = -0.2e1 * t228;
	t332 = -0.2e1 * t274;
	t200 = t274 * qJD(3) + t221 * t261 + t258 * t291;
	t220 = -qJD(1) * t295 - t263 * t307 + (t256 * t287 + t288) * t259;
	t187 = t211 * qJD(4) + t200 * t257 + t220 * t260;
	t210 = t234 * t257 + t260 * t273;
	t204 = t210 ^ 2;
	t198 = t204 * t206 + 0.1e1;
	t323 = t206 * t210;
	t305 = qJD(4) * t210;
	t188 = t200 * t260 - t220 * t257 - t305;
	t326 = t188 * t205 * t206;
	t328 = (t187 * t323 - t204 * t326) / t198 ^ 2;
	t320 = t236 * t335;
	t327 = (t201 * t318 - t224 * t320) / t216 ^ 2;
	t325 = t193 * t274;
	t322 = t212 * t274;
	t321 = t213 * t274;
	t319 = t228 * t236;
	t317 = t273 * t258;
	t316 = t273 * t261;
	t312 = t257 * t205;
	t311 = t260 * t210;
	t308 = qJD(2) * t259;
	t303 = -0.2e1 * t328;
	t302 = 0.2e1 * t328;
	t301 = -0.2e1 * t327;
	t300 = t194 * t332;
	t299 = t236 * t327;
	t298 = t193 * t322;
	t297 = t193 * t321;
	t296 = t210 * t326;
	t286 = 0.2e1 * t296;
	t285 = t320 * t333;
	t230 = -t258 * t313 - t261 * t275;
	t281 = -qJD(4) * t316 + t221;
	t209 = -t230 * t260 + t241 * t257;
	t208 = -t230 * t257 - t241 * t260;
	t278 = t206 * t311 - t312;
	t277 = -t230 * t236 + t240 * t318;
	t271 = -t212 + (t213 * t319 + t212) * t214;
	t270 = -qJD(3) * t317 + qJD(4) * t244 + t220 * t261;
	t227 = -t239 * qJD(3) + t261 * t290;
	t222 = t273 * qJD(1) + t275 * qJD(2);
	t219 = t244 * t257 + t260 * t316;
	t218 = -t244 * t260 + t257 * t316;
	t196 = 0.1e1 / t198;
	t189 = 0.1e1 / t191;
	t186 = t214 * t334;
	t185 = t277 * t214;
	t182 = t271 * t274;
	t180 = (-t212 * t241 + t213 * t314) * t258 + t280 * t186;
	t179 = t280 * t185 - t212 * t230 + t213 * t240;
	t177 = t277 * t301 + (t240 * t285 - t202 * t236 + (t201 * t240 + t226 * t230 + t227 * t228) * t237) * t214;
	t175 = t301 * t334 + (t276 * t306 + (t285 * t314 - t222 * t236 + (t226 * t241 + (t201 * t262 - t228 * t308) * t255) * t237) * t258) * t214;
	t1 = [t299 * t332 + (-t199 * t236 - t274 * t335) * t214, t175, t177, 0, 0; t228 * t192 * t304 + (-t201 * t192 + (t178 * t228 + t182 * t199) * t193) * t189 - (-t182 * t193 * t304 + (-0.2e1 * t182 * t329 + (-t183 * t214 * t319 + t301) * t298 + (t299 * t333 - t183 + (t183 - t279) * t214) * t297 - t271 * t324) * t189) * t274, (-t180 * t325 - t192 * t317) * t304 + (-t180 * t324 + (t220 * t258 + t273 * t306) * t192 + (t180 * t300 - t193 * t317) * t178 + (-t175 * t228 - t186 * t201 + (-t258 * t308 + t262 * t306) * t255 + (-t186 * t239 - t241 * t258) * t183) * t297 + (-t241 * t306 - t175 * t239 - t186 * t226 - t222 * t258 + (t186 * t228 - t258 * t314) * t183) * t298) * t189, (-t179 * t325 - t192 * t234) * t304 + (t179 * t178 * t300 + t200 * t192 + (-t234 * t178 - t179 * t199 + (-t177 * t228 - t185 * t201 + t227 + (-t185 * t239 - t230) * t183) * t321 + (-t177 * t239 - t185 * t226 - t202 + (t185 * t228 - t240) * t183) * t322) * t193) * t189, 0, 0; (-t205 * t208 + t209 * t323) * t302 + ((t209 * qJD(4) - t202 * t257 - t222 * t260) * t205 + t209 * t286 + (-t208 * t188 - (-t208 * qJD(4) - t202 * t260 + t222 * t257) * t210 - t209 * t187) * t206) * t196, (-t205 * t218 + t219 * t323) * t302 + (t219 * t286 - t281 * t205 * t260 + t270 * t312 + (-t281 * t210 * t257 - t219 * t187 - t218 * t188 - t270 * t311) * t206) * t196, -t278 * t274 * t303 + (t278 * t199 - ((-qJD(4) * t205 - 0.2e1 * t296) * t260 + (t187 * t260 + (t188 - t305) * t257) * t206) * t274) * t196, t303 + 0.2e1 * (t187 * t206 * t196 + (-t196 * t326 - t206 * t328) * t210) * t210, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:45:35
	% EndTime: 2019-12-31 22:45:37
	% DurationCPUTime: 1.48s
	% Computational Cost: add. (5342->150), mult. (14155->298), div. (744->12), fcn. (17850->13), ass. (0->132)
	t283 = cos(pkin(5));
	t285 = sin(qJ(2));
	t359 = sin(qJ(1));
	t322 = t359 * t285;
	t308 = t283 * t322;
	t316 = t359 * qJD(2);
	t287 = cos(qJ(2));
	t288 = cos(qJ(1));
	t337 = t288 * t287;
	t282 = sin(pkin(5));
	t339 = t282 * t288;
	t364 = -qJD(1) * t308 - t285 * t316 + (qJD(2) * t283 + qJD(1)) * t337 - qJD(3) * t339;
	t284 = sin(qJ(3));
	t286 = cos(qJ(3));
	t321 = t359 * t287;
	t338 = t288 * t285;
	t300 = -t283 * t338 - t321;
	t251 = -t284 * t300 + t286 * t339;
	t341 = t282 * t285;
	t262 = -t283 * t286 + t284 * t341;
	t240 = atan2(-t251, t262);
	t235 = sin(t240);
	t236 = cos(t240);
	t218 = -t235 * t251 + t236 * t262;
	t216 = 0.1e1 / t218 ^ 2;
	t267 = -t308 + t337;
	t323 = t282 * t359;
	t299 = -t267 * t284 + t286 * t323;
	t248 = t299 ^ 2;
	t214 = t248 * t216 + 0.1e1;
	t298 = -t283 * t321 - t338;
	t244 = t300 * qJD(1) + t298 * qJD(2);
	t257 = t267 * t286 + t284 * t323;
	t320 = qJD(1) * t339;
	t222 = t257 * qJD(3) + t244 * t284 - t286 * t320;
	t352 = t222 * t216;
	t247 = t251 ^ 2;
	t260 = 0.1e1 / t262 ^ 2;
	t239 = t247 * t260 + 0.1e1;
	t237 = 0.1e1 / t239;
	t317 = t359 * qJD(1);
	t307 = t282 * t317;
	t334 = qJD(3) * t286;
	t224 = t364 * t284 - t286 * t307 - t300 * t334;
	t263 = t283 * t284 + t286 * t341;
	t335 = qJD(2) * t287;
	t319 = t282 * t335;
	t249 = t263 * qJD(3) + t284 * t319;
	t259 = 0.1e1 / t262;
	t346 = t251 * t260;
	t304 = -t224 * t259 + t249 * t346;
	t206 = t304 * t237;
	t305 = -t235 * t262 - t236 * t251;
	t201 = t305 * t206 - t235 * t224 + t236 * t249;
	t215 = 0.1e1 / t218;
	t217 = t215 * t216;
	t357 = t201 * t217;
	t333 = 0.2e1 * (-t248 * t357 - t299 * t352) / t214 ^ 2;
	t363 = t249 * t260;
	t324 = t283 * t337;
	t264 = -t322 + t324;
	t340 = t282 * t287;
	t301 = -t259 * t264 + t340 * t346;
	t362 = t284 * t301;
	t225 = (qJD(3) * t300 + t307) * t284 + t364 * t286;
	t281 = qJ(4) + qJ(5);
	t278 = sin(t281);
	t279 = cos(t281);
	t234 = t257 * t279 - t278 * t298;
	t228 = 0.1e1 / t234;
	t229 = 0.1e1 / t234 ^ 2;
	t361 = -0.2e1 * t251;
	t360 = -0.2e1 * t299;
	t243 = -qJD(1) * t324 - t288 * t335 + (t283 * t316 + t317) * t285;
	t280 = qJD(4) + qJD(5);
	t311 = t257 * t280 + t243;
	t223 = t299 * qJD(3) + t244 * t286 + t284 * t320;
	t313 = -t280 * t298 + t223;
	t209 = t313 * t278 + t311 * t279;
	t233 = t257 * t278 + t279 * t298;
	t227 = t233 ^ 2;
	t221 = t227 * t229 + 0.1e1;
	t351 = t229 * t233;
	t210 = -t311 * t278 + t313 * t279;
	t354 = t210 * t228 * t229;
	t356 = (t209 * t351 - t227 * t354) / t221 ^ 2;
	t348 = t259 * t363;
	t355 = (t224 * t346 - t247 * t348) / t239 ^ 2;
	t353 = t216 * t299;
	t350 = t235 * t299;
	t349 = t236 * t299;
	t347 = t251 * t259;
	t345 = t298 * t284;
	t344 = t298 * t286;
	t343 = t278 * t228;
	t342 = t279 * t233;
	t336 = qJD(2) * t285;
	t332 = -0.2e1 * t356;
	t331 = 0.2e1 * t356;
	t330 = -0.2e1 * t355;
	t329 = t217 * t360;
	t328 = t259 * t355;
	t327 = t216 * t350;
	t326 = t216 * t349;
	t325 = t233 * t354;
	t315 = 0.2e1 * t325;
	t314 = t348 * t361;
	t312 = t264 * t280 - t225;
	t245 = t298 * qJD(1) + t300 * qJD(2);
	t253 = -t284 * t339 - t286 * t300;
	t310 = -t253 * t280 - t245;
	t306 = -t280 * t344 + t244;
	t303 = t229 * t342 - t343;
	t302 = -t253 * t259 + t263 * t346;
	t296 = -t235 + (t236 * t347 + t235) * t237;
	t295 = -qJD(3) * t345 + t243 * t286 + t267 * t280;
	t250 = -t262 * qJD(3) + t286 * t319;
	t242 = t267 * t278 + t279 * t344;
	t241 = -t267 * t279 + t278 * t344;
	t232 = -t253 * t279 + t264 * t278;
	t231 = -t253 * t278 - t264 * t279;
	t219 = 0.1e1 / t221;
	t212 = 0.1e1 / t214;
	t211 = t237 * t362;
	t208 = t302 * t237;
	t205 = t296 * t299;
	t203 = (-t235 * t264 + t236 * t340) * t284 + t305 * t211;
	t202 = t305 * t208 - t235 * t253 + t236 * t263;
	t200 = t302 * t330 + (t263 * t314 - t225 * t259 + (t224 * t263 + t249 * t253 + t250 * t251) * t260) * t237;
	t198 = t330 * t362 + (t301 * t334 + (t314 * t340 - t245 * t259 + (t249 * t264 + (t224 * t287 - t251 * t336) * t282) * t260) * t284) * t237;
	t197 = t332 + 0.2e1 * (t209 * t229 * t219 + (-t219 * t354 - t229 * t356) * t233) * t233;
	t1 = [t328 * t360 + (-t222 * t259 - t299 * t363) * t237, t198, t200, 0, 0; t251 * t215 * t333 + (-t224 * t215 + (t201 * t251 + t205 * t222) * t216) * t212 - (-t205 * t216 * t333 + (-0.2e1 * t205 * t357 + (-t206 * t237 * t347 + t330) * t327 + (t328 * t361 - t206 + (t206 - t304) * t237) * t326 - t296 * t352) * t212) * t299, (-t203 * t353 - t215 * t345) * t333 + (-t203 * t352 + (t243 * t284 + t298 * t334) * t215 + (t203 * t329 - t216 * t345) * t201 + (-t198 * t251 - t211 * t224 + (-t284 * t336 + t287 * t334) * t282 + (-t211 * t262 - t264 * t284) * t206) * t326 + (-t264 * t334 - t198 * t262 - t211 * t249 - t245 * t284 + (t211 * t251 - t284 * t340) * t206) * t327) * t212, (-t202 * t353 - t215 * t257) * t333 + (t202 * t201 * t329 + t223 * t215 + (-t257 * t201 - t202 * t222 + (-t200 * t251 - t208 * t224 + t250 + (-t208 * t262 - t253) * t206) * t349 + (-t200 * t262 - t208 * t249 - t225 + (t208 * t251 - t263) * t206) * t350) * t216) * t212, 0, 0; (-t228 * t231 + t232 * t351) * t331 + ((t312 * t278 + t310 * t279) * t228 + t232 * t315 + (-t231 * t210 - (-t310 * t278 + t312 * t279) * t233 - t232 * t209) * t229) * t219, (-t228 * t241 + t242 * t351) * t331 + (t242 * t315 - t306 * t228 * t279 + t295 * t343 + (-t306 * t233 * t278 - t242 * t209 - t241 * t210 - t295 * t342) * t229) * t219, -t303 * t299 * t332 + (t303 * t222 - ((-t228 * t280 - 0.2e1 * t325) * t279 + (t209 * t279 + (-t233 * t280 + t210) * t278) * t229) * t299) * t219, t197, t197;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end