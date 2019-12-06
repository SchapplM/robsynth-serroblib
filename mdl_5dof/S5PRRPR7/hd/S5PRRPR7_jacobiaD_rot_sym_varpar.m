% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPR7
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
%   Wie in S5PRRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(9));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(9)) * cos(pkin(5));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:32
	% EndTime: 2019-12-05 16:38:32
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (756->55), mult. (2271->133), div. (423->14), fcn. (2956->11), ass. (0->65)
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t136 = sin(pkin(9));
	t166 = cos(pkin(5));
	t155 = t136 * t166;
	t165 = cos(pkin(9));
	t127 = -t139 * t155 + t165 * t141;
	t138 = sin(qJ(3));
	t140 = cos(qJ(3));
	t137 = sin(pkin(5));
	t160 = t136 * t137;
	t149 = -t127 * t138 + t140 * t160;
	t170 = t149 * qJD(3);
	t152 = t166 * t165;
	t123 = t136 * t139 - t141 * t152;
	t159 = t137 * t141;
	t113 = atan2(-t123, -t159);
	t111 = sin(t113);
	t112 = cos(t113);
	t98 = -t111 * t123 - t112 * t159;
	t95 = 0.1e1 / t98;
	t110 = t127 * t140 + t138 * t160;
	t106 = 0.1e1 / t110;
	t133 = 0.1e1 / t141;
	t107 = 0.1e1 / t110 ^ 2;
	t134 = 0.1e1 / t141 ^ 2;
	t96 = 0.1e1 / t98 ^ 2;
	t105 = t149 ^ 2;
	t102 = t105 * t107 + 0.1e1;
	t148 = -t165 * t139 - t141 * t155;
	t119 = t148 * qJD(2);
	t103 = t110 * qJD(3) + t119 * t138;
	t163 = t107 * t149;
	t104 = t119 * t140 + t170;
	t164 = t104 * t106 * t107;
	t169 = 0.1e1 / t102 ^ 2 * (-t103 * t163 - t105 * t164);
	t125 = t136 * t141 + t139 * t152;
	t161 = t134 * t139;
	t156 = t123 * t161;
	t150 = t125 * t133 + t156;
	t121 = t123 ^ 2;
	t132 = 0.1e1 / t137 ^ 2;
	t116 = t121 * t132 * t134 + 0.1e1;
	t114 = 0.1e1 / t116;
	t131 = 0.1e1 / t137;
	t162 = t114 * t131;
	t91 = t150 * t162;
	t168 = t123 * t91;
	t167 = t148 * t96;
	t158 = qJD(2) * t139;
	t157 = -0.2e1 * t169;
	t151 = -t106 * t138 - t140 * t163;
	t135 = t133 * t134;
	t122 = t148 ^ 2;
	t120 = t127 * qJD(2);
	t118 = t125 * qJD(2);
	t117 = t123 * qJD(2);
	t100 = 0.1e1 / t102;
	t97 = t95 * t96;
	t94 = t122 * t96 + 0.1e1;
	t90 = (qJD(2) * t156 + t118 * t133) * t162;
	t88 = (t137 * t139 - t168) * t112 + (t91 * t159 - t125) * t111;
	t87 = (-t123 * t90 + t137 * t158) * t112 + (t90 * t159 - t118) * t111;
	t86 = (-0.2e1 * t150 * (t118 * t123 * t134 + t121 * t135 * t158) * t132 / t116 ^ 2 + (t118 * t161 - t117 * t133 + (t125 * t161 + (0.2e1 * t135 * t139 ^ 2 + t133) * t123) * qJD(2)) * t114) * t131;
	t1 = [0, t86, 0, 0, 0; 0, 0.2e1 * (-t127 * t95 - t88 * t167) / t94 ^ 2 * (-t122 * t87 * t97 - t120 * t167) + (-t88 * t120 * t96 + t119 * t95 + (-0.2e1 * t148 * t88 * t97 - t127 * t96) * t87 - (-(-t118 * t91 - t123 * t86 - t125 * t90 + (t90 * t91 + qJD(2)) * t159) * t112 - (t90 * t168 + t117 + (t141 * t86 + (-qJD(2) * t91 - t90) * t139) * t137) * t111) * t167) / t94, 0, 0, 0; 0, -t151 * t148 * t157 + (t151 * t120 - ((-qJD(3) * t106 + 0.2e1 * t149 * t164) * t140 + (t103 * t140 + (t104 + t170) * t138) * t107) * t148) * t100, t157 - 0.2e1 * (t100 * t103 * t107 - (-t100 * t164 - t107 * t169) * t149) * t149, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:32
	% EndTime: 2019-12-05 16:38:33
	% DurationCPUTime: 0.93s
	% Computational Cost: add. (2688->101), mult. (8196->213), div. (535->12), fcn. (10572->13), ass. (0->98)
	t187 = sin(pkin(9));
	t190 = cos(pkin(9));
	t193 = sin(qJ(2));
	t191 = cos(pkin(5));
	t195 = cos(qJ(2));
	t218 = t191 * t195;
	t176 = -t187 * t193 + t190 * t218;
	t169 = t176 * qJD(2);
	t219 = t191 * t193;
	t177 = t187 * t195 + t190 * t219;
	t192 = sin(qJ(3));
	t188 = sin(pkin(5));
	t222 = t188 * t192;
	t209 = t190 * t222;
	t194 = cos(qJ(3));
	t215 = qJD(3) * t194;
	t142 = -qJD(3) * t209 + t169 * t192 + t177 * t215;
	t221 = t188 * t194;
	t162 = t177 * t192 + t190 * t221;
	t160 = t162 ^ 2;
	t180 = -t191 * t194 + t193 * t222;
	t174 = 0.1e1 / t180 ^ 2;
	t156 = t160 * t174 + 0.1e1;
	t154 = 0.1e1 / t156;
	t181 = t191 * t192 + t193 * t221;
	t216 = qJD(2) * t195;
	t208 = t188 * t216;
	t167 = t181 * qJD(3) + t192 * t208;
	t173 = 0.1e1 / t180;
	t226 = t162 * t174;
	t126 = (-t142 * t173 + t167 * t226) * t154;
	t157 = atan2(-t162, t180);
	t152 = sin(t157);
	t153 = cos(t157);
	t206 = -t152 * t180 - t153 * t162;
	t122 = t206 * t126 - t142 * t152 + t153 * t167;
	t136 = -t152 * t162 + t153 * t180;
	t133 = 0.1e1 / t136;
	t134 = 0.1e1 / t136 ^ 2;
	t238 = t122 * t133 * t134;
	t210 = t187 * t219;
	t179 = t190 * t195 - t210;
	t166 = t179 * t194 + t187 * t222;
	t178 = t187 * t218 + t190 * t193;
	t186 = sin(pkin(10));
	t189 = cos(pkin(10));
	t150 = t166 * t186 - t178 * t189;
	t171 = t178 * qJD(2);
	t204 = -t179 * t192 + t187 * t221;
	t145 = t204 * qJD(3) - t171 * t194;
	t172 = -qJD(2) * t210 + t190 * t216;
	t141 = t145 * t189 + t172 * t186;
	t151 = t166 * t189 + t178 * t186;
	t147 = 0.1e1 / t151;
	t148 = 0.1e1 / t151 ^ 2;
	t232 = t141 * t147 * t148;
	t237 = 0.2e1 * t150 * t232;
	t236 = -0.2e1 * t204 * t238;
	t220 = t188 * t195;
	t202 = -t173 * t176 + t220 * t226;
	t235 = t192 * t202;
	t225 = t167 * t173 * t174;
	t234 = -0.2e1 * (t142 * t226 - t160 * t225) / t156 ^ 2;
	t233 = t134 * t204;
	t231 = t147 * t186;
	t230 = t148 * t150;
	t229 = t150 * t189;
	t228 = t152 * t204;
	t227 = t153 * t204;
	t224 = t178 * t192;
	t223 = t178 * t194;
	t217 = qJD(2) * t193;
	t161 = t204 ^ 2;
	t132 = t134 * t161 + 0.1e1;
	t144 = t166 * qJD(3) - t171 * t192;
	t214 = 0.2e1 * (-t144 * t233 - t161 * t238) / t132 ^ 2;
	t146 = t150 ^ 2;
	t139 = t146 * t148 + 0.1e1;
	t140 = t145 * t186 - t172 * t189;
	t213 = 0.2e1 * (t140 * t230 - t146 * t232) / t139 ^ 2;
	t207 = -0.2e1 * t162 * t225;
	t164 = t177 * t194 - t209;
	t205 = -t164 * t173 + t181 * t226;
	t203 = qJD(3) * t224 - t172 * t194;
	t170 = t177 * qJD(2);
	t168 = -t180 * qJD(3) + t194 * t208;
	t159 = t179 * t186 - t189 * t223;
	t158 = -t179 * t189 - t186 * t223;
	t143 = -t162 * qJD(3) + t169 * t194;
	t137 = 0.1e1 / t139;
	t129 = 0.1e1 / t132;
	t128 = t154 * t235;
	t127 = t205 * t154;
	t124 = (-t152 * t176 + t153 * t220) * t192 + t206 * t128;
	t123 = t206 * t127 - t152 * t164 + t153 * t181;
	t121 = t205 * t234 + (t181 * t207 - t143 * t173 + (t142 * t181 + t162 * t168 + t164 * t167) * t174) * t154;
	t119 = t234 * t235 + (t202 * t215 + (t207 * t220 + t170 * t173 + (t167 * t176 + (t142 * t195 - t162 * t217) * t188) * t174) * t192) * t154;
	t1 = [0, t119, t121, 0, 0; 0, (-t124 * t233 + t133 * t224) * t214 + ((-t172 * t192 - t178 * t215) * t133 + t124 * t236 + (-t124 * t144 + t224 * t122 + (-t119 * t162 - t128 * t142 + (-t192 * t217 + t195 * t215) * t188 + (-t128 * t180 - t176 * t192) * t126) * t227 + (-t176 * t215 - t119 * t180 - t128 * t167 + t170 * t192 + (t128 * t162 - t192 * t220) * t126) * t228) * t134) * t129, (-t123 * t233 - t133 * t166) * t214 + (t123 * t236 + t145 * t133 + (-t166 * t122 - t123 * t144 + (-t121 * t162 - t127 * t142 + t168 + (-t127 * t180 - t164) * t126) * t227 + (-t121 * t180 - t127 * t167 - t143 + (t127 * t162 - t181) * t126) * t228) * t134) * t129, 0, 0; 0, (-t147 * t158 + t159 * t230) * t213 + ((t171 * t189 + t203 * t186) * t147 + t159 * t237 + (-t158 * t141 - (-t171 * t186 + t203 * t189) * t150 - t159 * t140) * t148) * t137, -(-t148 * t229 + t231) * t204 * t213 + (t204 * t189 * t237 - t144 * t231 + (t144 * t229 - (t140 * t189 + t141 * t186) * t204) * t148) * t137, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:32
	% EndTime: 2019-12-05 16:38:33
	% DurationCPUTime: 1.24s
	% Computational Cost: add. (5037->126), mult. (14827->267), div. (532->12), fcn. (19143->15), ass. (0->118)
	t279 = cos(qJ(2));
	t335 = cos(pkin(9));
	t336 = cos(pkin(5));
	t302 = t336 * t335;
	t334 = sin(pkin(9));
	t337 = sin(qJ(2));
	t298 = -t279 * t302 + t334 * t337;
	t301 = t336 * t334;
	t267 = t335 * t279 - t337 * t301;
	t287 = -t334 * t279 - t337 * t302;
	t339 = qJD(2) * t287;
	t272 = sin(pkin(10));
	t274 = cos(pkin(10));
	t276 = sin(qJ(3));
	t278 = cos(qJ(3));
	t273 = sin(pkin(5));
	t311 = t273 * t335;
	t292 = t276 * t311 + t278 * t287;
	t240 = -t292 * t272 - t298 * t274;
	t238 = t240 ^ 2;
	t313 = t278 * t337;
	t289 = t273 * t313 + t336 * t276;
	t319 = t279 * t274;
	t253 = t289 * t272 + t273 * t319;
	t251 = 0.1e1 / t253 ^ 2;
	t224 = t238 * t251 + 0.1e1;
	t254 = t276 * t287 - t278 * t311;
	t262 = t298 * qJD(2);
	t226 = (t254 * qJD(3) - t262 * t278) * t272 + t274 * t339;
	t309 = t336 * t278;
	t312 = t337 * t274;
	t314 = t276 * t337;
	t318 = qJD(2) * t279;
	t244 = qJD(3) * t272 * t309 + ((-qJD(3) * t314 + t278 * t318) * t272 - qJD(2) * t312) * t273;
	t250 = 0.1e1 / t253;
	t326 = t244 * t250 * t251;
	t327 = t240 * t251;
	t338 = -0.2e1 * (t226 * t327 - t238 * t326) / t224 ^ 2;
	t225 = atan2(-t240, t253);
	t220 = sin(t225);
	t221 = cos(t225);
	t299 = t220 * t240 - t221 * t253;
	t205 = 0.1e1 / t299;
	t310 = t273 * t334;
	t256 = t267 * t278 + t276 * t310;
	t286 = -t279 * t301 - t335 * t337;
	t243 = t256 * t274 - t272 * t286;
	t275 = sin(qJ(5));
	t277 = cos(qJ(5));
	t291 = -t267 * t276 + t278 * t310;
	t219 = t243 * t277 - t275 * t291;
	t215 = 0.1e1 / t219;
	t206 = 0.1e1 / t299 ^ 2;
	t216 = 0.1e1 / t219 ^ 2;
	t222 = 0.1e1 / t224;
	t198 = (-t226 * t250 + t244 * t327) * t222;
	t300 = -t220 * t253 - t221 * t240;
	t194 = t300 * t198 - t220 * t226 + t221 * t244;
	t333 = t194 * t205 * t206;
	t263 = t286 * qJD(2);
	t237 = t291 * qJD(3) + t263 * t278;
	t264 = t267 * qJD(2);
	t228 = t237 * t274 + t264 * t272;
	t236 = t256 * qJD(3) + t263 * t276;
	t209 = t219 * qJD(5) + t228 * t275 - t236 * t277;
	t218 = t243 * t275 + t277 * t291;
	t214 = t218 ^ 2;
	t213 = t214 * t216 + 0.1e1;
	t329 = t216 * t218;
	t210 = -t218 * qJD(5) + t228 * t277 + t236 * t275;
	t330 = t210 * t215 * t216;
	t332 = (t209 * t329 - t214 * t330) / t213 ^ 2;
	t242 = t256 * t272 + t274 * t286;
	t331 = t206 * t242;
	t227 = t237 * t272 - t264 * t274;
	t328 = t227 * t206;
	t325 = t286 * t276;
	t324 = t286 * t278;
	t323 = t272 * t278;
	t322 = t272 * t279;
	t321 = t274 * t275;
	t320 = t274 * t277;
	t317 = qJD(3) * t276;
	t239 = t242 ^ 2;
	t203 = t239 * t206 + 0.1e1;
	t316 = 0.2e1 * (t239 * t333 + t242 * t328) / t203 ^ 2;
	t315 = 0.2e1 * t332;
	t307 = 0.2e1 * t218 * t330;
	t306 = -0.2e1 * t240 * t326;
	t305 = -0.2e1 * t242 * t333;
	t294 = -t264 * t278 - t286 * t317;
	t303 = qJD(5) * t325 + t263 * t272 + t294 * t274;
	t290 = t272 * t298;
	t245 = t274 * t287 - t278 * t290;
	t259 = (t278 * t322 - t312) * t273;
	t297 = -t245 * t250 + t259 * t327;
	t268 = -t273 * t314 + t309;
	t296 = -t250 * t254 + t268 * t327;
	t232 = t256 * t275 + t291 * t320;
	t295 = t256 * t277 - t291 * t321;
	t247 = t267 * t272 + t274 * t324;
	t288 = -qJD(3) * t324 + qJD(5) * t247 + t264 * t276;
	t257 = -t273 * t276 * t318 - t289 * qJD(3);
	t248 = (-t317 * t322 + (-t272 * t313 - t319) * qJD(2)) * t273;
	t246 = -t267 * t274 + t286 * t323;
	t235 = t292 * qJD(3) + t262 * t276;
	t234 = t247 * t277 + t275 * t325;
	t233 = t247 * t275 - t277 * t325;
	t229 = t262 * t274 + t290 * t317 + t323 * t339;
	t211 = 0.1e1 / t213;
	t201 = 0.1e1 / t203;
	t200 = t296 * t272 * t222;
	t199 = t297 * t222;
	t196 = (-t220 * t254 + t221 * t268) * t272 + t300 * t200;
	t195 = t300 * t199 - t220 * t245 + t221 * t259;
	t193 = (t296 * t338 + (t268 * t306 - t235 * t250 + (t226 * t268 + t240 * t257 + t244 * t254) * t251) * t222) * t272;
	t192 = t297 * t338 + (t259 * t306 - t229 * t250 + (t226 * t259 + t240 * t248 + t244 * t245) * t251) * t222;
	t1 = [0, t192, t193, 0, 0; 0, (t195 * t331 + t205 * t246) * t316 + (-(-t263 * t274 + t294 * t272) * t205 + t195 * t305 + (-t246 * t194 - t195 * t227 + (-(-t192 * t240 - t199 * t226 + t248 + (-t199 * t253 - t245) * t198) * t221 - (-t192 * t253 - t199 * t244 - t229 + (t199 * t240 - t259) * t198) * t220) * t242) * t206) * t201, (t205 * t272 * t291 + t196 * t331) * t316 + (-(t300 * t193 + (t299 * t198 - t220 * t244 - t221 * t226) * t200) * t331 + (t305 - t328) * t196 + (t236 * t205 + (-t291 * t194 - (-t220 * t235 + t221 * t257 + (-t220 * t268 - t221 * t254) * t198) * t242) * t206) * t272) * t201, 0, 0; 0, (-t215 * t233 + t234 * t329) * t315 + (t234 * t307 + (-t234 * t209 - t233 * t210 + (t275 * t288 - t277 * t303) * t218) * t216 + (t303 * t275 + t288 * t277) * t215) * t211, (t215 * t295 + t232 * t329) * t315 + ((t232 * qJD(5) - t236 * t321 - t237 * t277) * t215 + t232 * t307 + (t295 * t210 - (t295 * qJD(5) - t236 * t320 + t237 * t275) * t218 - t232 * t209) * t216) * t211, 0, -0.2e1 * t332 + 0.2e1 * (t209 * t216 * t211 + (-t211 * t330 - t216 * t332) * t218) * t218;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end