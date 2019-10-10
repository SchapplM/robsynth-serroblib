% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPP3
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
%   Wie in S6PRRRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(10));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(10)) * cos(pkin(6));
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
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:51
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (756->55), mult. (2271->133), div. (423->14), fcn. (2956->11), ass. (0->65)
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t136 = sin(pkin(10));
	t166 = cos(pkin(6));
	t155 = t136 * t166;
	t165 = cos(pkin(10));
	t127 = -t139 * t155 + t165 * t141;
	t138 = sin(qJ(3));
	t140 = cos(qJ(3));
	t137 = sin(pkin(6));
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
	t117 = qJD(2) * t123;
	t100 = 0.1e1 / t102;
	t97 = t95 * t96;
	t94 = t122 * t96 + 0.1e1;
	t90 = (qJD(2) * t156 + t118 * t133) * t162;
	t88 = (t137 * t139 - t168) * t112 + (t91 * t159 - t125) * t111;
	t87 = (-t123 * t90 + t137 * t158) * t112 + (t90 * t159 - t118) * t111;
	t86 = (-0.2e1 * t150 * (t118 * t123 * t134 + t121 * t135 * t158) * t132 / t116 ^ 2 + (t118 * t161 - t117 * t133 + (t125 * t161 + (0.2e1 * t135 * t139 ^ 2 + t133) * t123) * qJD(2)) * t114) * t131;
	t1 = [0, t86, 0, 0, 0, 0; 0, 0.2e1 * (-t127 * t95 - t88 * t167) / t94 ^ 2 * (-t122 * t97 * t87 - t120 * t167) + (-t88 * t120 * t96 + t119 * t95 + (-0.2e1 * t148 * t88 * t97 - t127 * t96) * t87 - (-(-t118 * t91 - t123 * t86 - t125 * t90 + (t90 * t91 + qJD(2)) * t159) * t112 - (t90 * t168 + t117 + (t141 * t86 + (-qJD(2) * t91 - t90) * t139) * t137) * t111) * t167) / t94, 0, 0, 0, 0; 0, -t151 * t148 * t157 + (t151 * t120 - ((-qJD(3) * t106 + 0.2e1 * t149 * t164) * t140 + (t103 * t140 + (t104 + t170) * t138) * t107) * t148) * t100, t157 - 0.2e1 * (t100 * t103 * t107 - (-t100 * t164 - t107 * t169) * t149) * t149, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:52
	% DurationCPUTime: 1.34s
	% Computational Cost: add. (3002->109), mult. (9085->226), div. (559->12), fcn. (11668->13), ass. (0->103)
	t207 = sin(pkin(10));
	t209 = cos(pkin(10));
	t213 = sin(qJ(2));
	t210 = cos(pkin(6));
	t216 = cos(qJ(2));
	t242 = t210 * t216;
	t197 = -t207 * t213 + t209 * t242;
	t190 = t197 * qJD(2);
	t243 = t210 * t213;
	t198 = t207 * t216 + t209 * t243;
	t212 = sin(qJ(3));
	t208 = sin(pkin(6));
	t246 = t208 * t212;
	t232 = t209 * t246;
	t215 = cos(qJ(3));
	t239 = qJD(3) * t215;
	t162 = -qJD(3) * t232 + t190 * t212 + t198 * t239;
	t245 = t208 * t215;
	t182 = t198 * t212 + t209 * t245;
	t180 = t182 ^ 2;
	t201 = -t210 * t215 + t213 * t246;
	t195 = 0.1e1 / t201 ^ 2;
	t176 = t180 * t195 + 0.1e1;
	t174 = 0.1e1 / t176;
	t202 = t210 * t212 + t213 * t245;
	t240 = qJD(2) * t216;
	t231 = t208 * t240;
	t187 = t202 * qJD(3) + t212 * t231;
	t194 = 0.1e1 / t201;
	t250 = t182 * t195;
	t146 = (-t162 * t194 + t187 * t250) * t174;
	t177 = atan2(-t182, t201);
	t172 = sin(t177);
	t173 = cos(t177);
	t228 = -t172 * t201 - t173 * t182;
	t142 = t228 * t146 - t172 * t162 + t173 * t187;
	t156 = -t172 * t182 + t173 * t201;
	t153 = 0.1e1 / t156;
	t154 = 0.1e1 / t156 ^ 2;
	t263 = t142 * t153 * t154;
	t233 = t207 * t243;
	t200 = t209 * t216 - t233;
	t225 = -t200 * t212 + t207 * t245;
	t262 = -0.2e1 * t225 * t263;
	t244 = t208 * t216;
	t224 = -t194 * t197 + t244 * t250;
	t261 = t212 * t224;
	t249 = t187 * t194 * t195;
	t260 = -0.2e1 * (t162 * t250 - t180 * t249) / t176 ^ 2;
	t186 = t200 * t215 + t207 * t246;
	t199 = t207 * t242 + t209 * t213;
	t211 = sin(qJ(4));
	t214 = cos(qJ(4));
	t171 = t186 * t214 + t199 * t211;
	t167 = 0.1e1 / t171;
	t168 = 0.1e1 / t171 ^ 2;
	t192 = t199 * qJD(2);
	t165 = t225 * qJD(3) - t192 * t215;
	t193 = -qJD(2) * t233 + t209 * t240;
	t157 = t171 * qJD(4) + t165 * t211 - t193 * t214;
	t170 = t186 * t211 - t199 * t214;
	t166 = t170 ^ 2;
	t161 = t166 * t168 + 0.1e1;
	t254 = t168 * t170;
	t238 = qJD(4) * t170;
	t158 = t165 * t214 + t193 * t211 - t238;
	t257 = t158 * t167 * t168;
	t259 = (t157 * t254 - t166 * t257) / t161 ^ 2;
	t258 = t154 * t225;
	t164 = t186 * qJD(3) - t192 * t212;
	t256 = t164 * t154;
	t255 = t167 * t211;
	t253 = t170 * t214;
	t252 = t172 * t225;
	t251 = t173 * t225;
	t248 = t199 * t212;
	t247 = t199 * t215;
	t241 = qJD(2) * t213;
	t181 = t225 ^ 2;
	t152 = t181 * t154 + 0.1e1;
	t237 = 0.2e1 * (-t181 * t263 - t225 * t256) / t152 ^ 2;
	t236 = -0.2e1 * t259;
	t234 = t170 * t257;
	t230 = -0.2e1 * t182 * t249;
	t229 = qJD(4) * t247 - t192;
	t227 = t168 * t253 - t255;
	t184 = t198 * t215 - t232;
	t226 = -t184 * t194 + t202 * t250;
	t223 = qJD(3) * t248 + qJD(4) * t200 - t193 * t215;
	t191 = t198 * qJD(2);
	t188 = -t201 * qJD(3) + t215 * t231;
	t179 = t200 * t211 - t214 * t247;
	t178 = -t200 * t214 - t211 * t247;
	t163 = -t182 * qJD(3) + t190 * t215;
	t159 = 0.1e1 / t161;
	t149 = 0.1e1 / t152;
	t148 = t174 * t261;
	t147 = t226 * t174;
	t144 = (-t172 * t197 + t173 * t244) * t212 + t228 * t148;
	t143 = t228 * t147 - t172 * t184 + t173 * t202;
	t141 = t226 * t260 + (t202 * t230 - t163 * t194 + (t162 * t202 + t182 * t188 + t184 * t187) * t195) * t174;
	t139 = t260 * t261 + (t224 * t239 + (t230 * t244 + t191 * t194 + (t187 * t197 + (t162 * t216 - t182 * t241) * t208) * t195) * t212) * t174;
	t1 = [0, t139, t141, 0, 0, 0; 0, (-t144 * t258 + t153 * t248) * t237 + ((-t193 * t212 - t199 * t239) * t153 + (-t256 + t262) * t144 + (t248 * t142 + (-t139 * t182 - t148 * t162 + (-t212 * t241 + t216 * t239) * t208 + (-t148 * t201 - t197 * t212) * t146) * t251 + (-t197 * t239 - t139 * t201 - t148 * t187 + t191 * t212 + (t148 * t182 - t212 * t244) * t146) * t252) * t154) * t149, (-t143 * t258 - t153 * t186) * t237 + (t143 * t262 + t165 * t153 + (-t186 * t142 - t143 * t164 + (-t141 * t182 - t147 * t162 + t188 + (-t147 * t201 - t184) * t146) * t251 + (-t141 * t201 - t147 * t187 - t163 + (t147 * t182 - t202) * t146) * t252) * t154) * t149, 0, 0, 0; 0, 0.2e1 * (-t167 * t178 + t179 * t254) * t259 + (0.2e1 * t179 * t234 - t229 * t167 * t214 + t223 * t255 + (-t229 * t170 * t211 - t179 * t157 - t178 * t158 - t223 * t253) * t168) * t159, -t227 * t225 * t236 + (t227 * t164 - ((-qJD(4) * t167 - 0.2e1 * t234) * t214 + (t157 * t214 + (t158 - t238) * t211) * t168) * t225) * t159, t236 + 0.2e1 * (t157 * t168 * t159 + (-t159 * t257 - t168 * t259) * t170) * t170, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:53
	% DurationCPUTime: 2.47s
	% Computational Cost: add. (7193->155), mult. (20987->306), div. (784->12), fcn. (27058->13), ass. (0->122)
	t249 = sin(qJ(3));
	t252 = cos(qJ(3));
	t250 = sin(qJ(2));
	t253 = cos(qJ(2));
	t309 = cos(pkin(10));
	t310 = cos(pkin(6));
	t275 = t310 * t309;
	t308 = sin(pkin(10));
	t262 = -t250 * t275 - t308 * t253;
	t247 = sin(pkin(6));
	t280 = t247 * t309;
	t220 = t249 * t262 - t252 * t280;
	t273 = t308 * t250 - t253 * t275;
	t236 = t273 * qJD(2);
	t203 = t220 * qJD(3) - t236 * t252;
	t248 = sin(qJ(4));
	t251 = cos(qJ(4));
	t264 = t249 * t280 + t252 * t262;
	t211 = t273 * t248 - t251 * t264;
	t260 = t262 * qJD(2);
	t185 = t211 * qJD(4) + t203 * t248 + t251 * t260;
	t265 = t273 * t251;
	t209 = -t248 * t264 - t265;
	t206 = t209 ^ 2;
	t294 = t247 * t250;
	t243 = t310 * t249 + t252 * t294;
	t292 = t251 * t253;
	t229 = t243 * t248 + t247 * t292;
	t225 = 0.1e1 / t229 ^ 2;
	t197 = t206 * t225 + 0.1e1;
	t195 = 0.1e1 / t197;
	t242 = -t249 * t294 + t310 * t252;
	t291 = qJD(2) * t247;
	t281 = t253 * t291;
	t228 = t242 * qJD(3) + t252 * t281;
	t293 = t248 * t253;
	t230 = t243 * t251 - t247 * t293;
	t282 = t250 * t291;
	t199 = t230 * qJD(4) + t228 * t248 - t251 * t282;
	t224 = 0.1e1 / t229;
	t301 = t209 * t225;
	t172 = (-t185 * t224 + t199 * t301) * t195;
	t198 = atan2(-t209, t229);
	t192 = sin(t198);
	t193 = cos(t198);
	t272 = -t192 * t229 - t193 * t209;
	t168 = t272 * t172 - t192 * t185 + t193 * t199;
	t184 = -t192 * t209 + t193 * t229;
	t181 = 0.1e1 / t184;
	t182 = 0.1e1 / t184 ^ 2;
	t313 = t168 * t181 * t182;
	t274 = t310 * t308;
	t241 = -t250 * t274 + t309 * t253;
	t279 = t247 * t308;
	t223 = t241 * t252 + t249 * t279;
	t261 = -t309 * t250 - t253 * t274;
	t212 = t223 * t248 + t251 * t261;
	t277 = 0.2e1 * t212 * t313;
	t267 = -t220 * t224 + t242 * t301;
	t312 = t248 * t267;
	t303 = t199 * t224 * t225;
	t311 = -0.2e1 * (t185 * t301 - t206 * t303) / t197 ^ 2;
	t263 = -t241 * t249 + t252 * t279;
	t217 = 0.1e1 / t263;
	t218 = 0.1e1 / t263 ^ 2;
	t307 = t182 * t212;
	t237 = t261 * qJD(2);
	t205 = t263 * qJD(3) + t237 * t252;
	t213 = t223 * t251 - t248 * t261;
	t238 = t241 * qJD(2);
	t187 = t213 * qJD(4) + t205 * t248 - t238 * t251;
	t306 = t187 * t182;
	t305 = t192 * t212;
	t304 = t193 * t212;
	t204 = t223 * qJD(3) + t237 * t249;
	t219 = t217 * t218;
	t302 = t204 * t219;
	t300 = t213 * t218;
	t299 = t213 * t223;
	t298 = t217 * t263;
	t297 = t263 * t248;
	t296 = t261 * t249;
	t295 = t261 * t252;
	t290 = qJD(2) * t252;
	t289 = qJD(3) * t249;
	t288 = qJD(4) * t248;
	t287 = qJD(4) * t251;
	t286 = t252 * qJD(4);
	t207 = t212 ^ 2;
	t180 = t207 * t182 + 0.1e1;
	t285 = 0.2e1 * (-t207 * t313 + t212 * t306) / t180 ^ 2;
	t188 = -t212 * qJD(4) + t205 * t251 + t238 * t248;
	t208 = t213 ^ 2;
	t194 = t208 * t218 + 0.1e1;
	t284 = 0.2e1 * (t188 * t300 + t208 * t302) / t194 ^ 2;
	t276 = -0.2e1 * t209 * t303;
	t271 = qJD(4) * t241 - t238 * t252;
	t269 = -t211 * t224 + t230 * t301;
	t266 = t252 * t273;
	t214 = -t248 * t266 + t251 * t262;
	t232 = (-t250 * t251 + t252 * t293) * t247;
	t268 = -t214 * t224 + t232 * t301;
	t227 = -t243 * qJD(3) - t249 * t281;
	t216 = t241 * t248 + t251 * t295;
	t215 = -t241 * t251 + t248 * t295;
	t202 = t264 * qJD(3) + t236 * t249;
	t201 = ((-qJD(2) + t286) * t292 + (-t253 * t289 + (qJD(4) - t290) * t250) * t248) * t247;
	t200 = -t229 * qJD(4) + t228 * t251 + t248 * t282;
	t190 = 0.1e1 / t194;
	t189 = t236 * t251 - t262 * t288 - t266 * t287 + (t262 * t290 + t273 * t289) * t248;
	t186 = qJD(4) * t265 + t203 * t251 - t248 * t260 + t264 * t288;
	t178 = 0.1e1 / t180;
	t177 = t195 * t312;
	t176 = t268 * t195;
	t174 = t269 * t195;
	t171 = (-t192 * t220 + t193 * t242) * t248 + t272 * t177;
	t170 = t272 * t176 - t192 * t214 + t193 * t232;
	t169 = t272 * t174 - t192 * t211 + t193 * t230;
	t167 = t268 * t311 + (t232 * t276 - t189 * t224 + (t185 * t232 + t199 * t214 + t201 * t209) * t225) * t195;
	t165 = t269 * t311 + (t230 * t276 - t186 * t224 + (t185 * t230 + t199 * t211 + t200 * t209) * t225) * t195;
	t164 = t311 * t312 + (t267 * t287 + (t242 * t276 - t202 * t224 + (t185 * t242 + t199 * t220 + t209 * t227) * t225) * t248) * t195;
	t1 = [0, t167, t164, t165, 0, 0; 0, (t170 * t307 - t181 * t215) * t285 + (t170 * t277 + (-t215 * t168 - t170 * t187 - (-t167 * t209 - t176 * t185 + t201 + (-t176 * t229 - t214) * t172) * t304 - (-t167 * t229 - t176 * t199 - t189 + (t176 * t209 - t232) * t172) * t305) * t182 + ((t261 * t286 - t237) * t251 + (-t261 * t289 + t271) * t248) * t181) * t178, (t171 * t307 - t181 * t297) * t285 + ((-t204 * t248 + t263 * t287) * t181 + (-t306 + t277) * t171 + (-t297 * t168 - (t242 * t287 - t164 * t209 - t177 * t185 + t227 * t248 + (-t177 * t229 - t220 * t248) * t172) * t304 - (-t220 * t287 - t164 * t229 - t177 * t199 - t202 * t248 + (t177 * t209 - t242 * t248) * t172) * t305) * t182) * t178, (t169 * t307 - t181 * t213) * t285 + (t169 * t277 + t188 * t181 + (-t213 * t168 - t169 * t187 - (-t165 * t209 - t174 * t185 + t200 + (-t174 * t229 - t211) * t172) * t304 - (-t165 * t229 - t174 * t199 - t186 + (t174 * t209 - t230) * t172) * t305) * t182) * t178, 0, 0; 0, (t216 * t217 + t296 * t300) * t284 + (-(t237 * t248 + t271 * t251) * t217 - (-(t248 * t286 + t251 * t289) * t217 + 0.2e1 * t249 * t213 * t302) * t261 + (-t188 * t296 - t216 * t204 + (-qJD(3) * t295 + t238 * t249) * t213) * t218) * t190, (t218 * t299 + t251 * t298) * t284 + (t288 * t298 + (-t188 * t223 - t205 * t213) * t218 + (-0.2e1 * t219 * t299 + (-t218 * t263 + t217) * t251) * t204) * t190, -t212 * t217 * t284 + (t204 * t212 * t218 + t187 * t217) * t190, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:51
	% EndTime: 2019-10-09 22:44:53
	% DurationCPUTime: 2.42s
	% Computational Cost: add. (7193->157), mult. (20987->307), div. (784->12), fcn. (27058->13), ass. (0->119)
	t247 = sin(pkin(10));
	t249 = cos(pkin(10));
	t256 = cos(qJ(2));
	t250 = cos(pkin(6));
	t253 = sin(qJ(2));
	t289 = t250 * t253;
	t239 = t247 * t256 + t249 * t289;
	t252 = sin(qJ(3));
	t248 = sin(pkin(6));
	t255 = cos(qJ(3));
	t290 = t248 * t255;
	t222 = -t239 * t252 - t249 * t290;
	t288 = t250 * t256;
	t264 = -t247 * t253 + t249 * t288;
	t234 = t264 * qJD(2);
	t204 = t222 * qJD(3) + t234 * t255;
	t235 = t239 * qJD(2);
	t251 = sin(qJ(4));
	t254 = cos(qJ(4));
	t292 = t248 * t252;
	t267 = -t239 * t255 + t249 * t292;
	t263 = t267 * t251;
	t187 = t235 * t251 + qJD(4) * t263 + (-qJD(4) * t264 + t204) * t254;
	t211 = -t251 * t264 - t267 * t254;
	t207 = t211 ^ 2;
	t243 = t250 * t252 + t253 * t290;
	t287 = t251 * t256;
	t265 = -t243 * t254 + t248 * t287;
	t226 = 0.1e1 / t265 ^ 2;
	t198 = t207 * t226 + 0.1e1;
	t196 = 0.1e1 / t198;
	t291 = t248 * t253;
	t242 = t250 * t255 - t252 * t291;
	t284 = qJD(2) * t256;
	t275 = t248 * t284;
	t229 = t242 * qJD(3) + t255 * t275;
	t285 = t254 * t256;
	t230 = -t243 * t251 - t248 * t285;
	t276 = qJD(2) * t291;
	t201 = t230 * qJD(4) + t229 * t254 + t251 * t276;
	t225 = 0.1e1 / t265;
	t298 = t211 * t226;
	t173 = (t187 * t225 + t201 * t298) * t196;
	t199 = atan2(-t211, -t265);
	t193 = sin(t199);
	t194 = cos(t199);
	t272 = t193 * t265 - t194 * t211;
	t169 = t272 * t173 - t193 * t187 + t194 * t201;
	t185 = -t193 * t211 - t194 * t265;
	t182 = 0.1e1 / t185;
	t183 = 0.1e1 / t185 ^ 2;
	t306 = t169 * t182 * t183;
	t277 = t247 * t289;
	t241 = t249 * t256 - t277;
	t224 = t241 * t255 + t247 * t292;
	t240 = t247 * t288 + t249 * t253;
	t214 = t224 * t254 + t240 * t251;
	t274 = 0.2e1 * t214 * t306;
	t268 = t222 * t225 + t242 * t298;
	t305 = t254 * t268;
	t300 = t201 * t225 * t226;
	t304 = -0.2e1 * (t187 * t298 + t207 * t300) / t198 ^ 2;
	t266 = -t241 * t252 + t247 * t290;
	t219 = 0.1e1 / t266;
	t220 = 0.1e1 / t266 ^ 2;
	t303 = t183 * t214;
	t302 = t193 * t214;
	t301 = t194 * t214;
	t236 = t240 * qJD(2);
	t205 = t224 * qJD(3) - t236 * t252;
	t221 = t219 * t220;
	t299 = t205 * t221;
	t213 = -t224 * t251 + t240 * t254;
	t297 = t213 * t220;
	t296 = t213 * t224;
	t295 = t266 * t254;
	t294 = t240 * t252;
	t293 = t240 * t255;
	t286 = t254 * t255;
	t283 = qJD(3) * t252;
	t282 = qJD(4) * t251;
	t281 = qJD(4) * t255;
	t209 = t214 ^ 2;
	t181 = t183 * t209 + 0.1e1;
	t206 = t266 * qJD(3) - t236 * t255;
	t237 = -qJD(2) * t277 + t249 * t284;
	t189 = t213 * qJD(4) + t206 * t254 + t237 * t251;
	t280 = 0.2e1 * (t189 * t303 - t209 * t306) / t181 ^ 2;
	t188 = -t214 * qJD(4) - t206 * t251 + t237 * t254;
	t208 = t213 ^ 2;
	t195 = t208 * t220 + 0.1e1;
	t279 = 0.2e1 * (t188 * t297 + t208 * t299) / t195 ^ 2;
	t273 = 0.2e1 * t211 * t300;
	t271 = -qJD(4) * t241 + t237 * t255;
	t210 = t254 * t264 - t263;
	t270 = -t210 * t225 + t230 * t298;
	t215 = t239 * t251 + t264 * t286;
	t232 = (t251 * t253 + t255 * t285) * t248;
	t269 = t215 * t225 + t232 * t298;
	t228 = -t243 * qJD(3) - t252 * t275;
	t217 = -t240 * t286 + t241 * t251;
	t216 = t241 * t254 + t251 * t293;
	t203 = t267 * qJD(3) - t234 * t252;
	t202 = ((qJD(2) - t281) * t287 + (-t256 * t283 + (-qJD(2) * t255 + qJD(4)) * t253) * t254) * t248;
	t200 = t265 * qJD(4) - t229 * t251 + t254 * t276;
	t191 = 0.1e1 / t195;
	t190 = (-t264 * t281 + t234) * t251 + (qJD(4) * t239 - t235 * t255 - t264 * t283) * t254;
	t186 = t211 * qJD(4) + t204 * t251 - t235 * t254;
	t179 = 0.1e1 / t181;
	t178 = t196 * t305;
	t177 = t269 * t196;
	t175 = t270 * t196;
	t172 = (-t193 * t222 + t194 * t242) * t254 + t272 * t178;
	t171 = t272 * t177 - t193 * t215 + t194 * t232;
	t170 = t272 * t175 + t193 * t210 + t194 * t230;
	t168 = t269 * t304 + (t232 * t273 + t190 * t225 + (t187 * t232 + t201 * t215 + t202 * t211) * t226) * t196;
	t166 = t270 * t304 + (t230 * t273 - t186 * t225 + (t187 * t230 + t200 * t211 - t201 * t210) * t226) * t196;
	t165 = t304 * t305 + (-t268 * t282 + (t242 * t273 + t203 * t225 + (t187 * t242 + t201 * t222 + t211 * t228) * t226) * t254) * t196;
	t1 = [0, t168, t165, t166, 0, 0; 0, (t171 * t303 - t182 * t217) * t280 + (t171 * t274 + (-t217 * t169 - t171 * t189 - (-t168 * t211 - t177 * t187 + t202 + (t177 * t265 - t215) * t173) * t301 - (t168 * t265 - t177 * t201 - t190 + (t177 * t211 - t232) * t173) * t302) * t183 + ((t240 * t281 - t236) * t251 + (t240 * t283 - t271) * t254) * t182) * t179, (t172 * t303 - t182 * t295) * t280 + ((-t205 * t254 - t266 * t282) * t182 + t172 * t274 + (-t172 * t189 - t295 * t169 - (-t242 * t282 - t165 * t211 - t178 * t187 + t228 * t254 + (t178 * t265 - t222 * t254) * t173) * t301 - (t222 * t282 + t165 * t265 - t178 * t201 - t203 * t254 + (t178 * t211 - t242 * t254) * t173) * t302) * t183) * t179, (t170 * t303 - t182 * t213) * t280 + (t170 * t274 + t188 * t182 + (-t213 * t169 - t170 * t189 - (-t166 * t211 - t175 * t187 + t200 + (t175 * t265 + t210) * t173) * t301 - (t166 * t265 - t175 * t201 + t186 + (t175 * t211 - t230) * t173) * t302) * t183) * t179, 0, 0; 0, (t216 * t219 - t294 * t297) * t279 + (-(-t236 * t254 + t271 * t251) * t219 + (-(-t251 * t283 + t254 * t281) * t219 + 0.2e1 * t252 * t213 * t299) * t240 + (t188 * t294 - t216 * t205 + (qJD(3) * t293 + t237 * t252) * t213) * t220) * t191, (-t219 * t251 * t266 + t220 * t296) * t279 + (qJD(4) * t219 * t295 + (-t188 * t224 - t206 * t213) * t220 + (-0.2e1 * t221 * t296 + (t220 * t266 - t219) * t251) * t205) * t191, -t214 * t219 * t279 + (t205 * t214 * t220 + t189 * t219) * t191, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end