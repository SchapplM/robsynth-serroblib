% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP2
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
%   Wie in S6PRRRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(11));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(11)) * cos(pkin(6));
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
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:57
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (756->55), mult. (2271->133), div. (423->14), fcn. (2956->11), ass. (0->65)
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t136 = sin(pkin(11));
	t166 = cos(pkin(6));
	t155 = t136 * t166;
	t165 = cos(pkin(11));
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
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:57
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (1157->57), mult. (2577->133), div. (441->14), fcn. (3313->11), ass. (0->69)
	t172 = sin(qJ(2));
	t173 = cos(qJ(2));
	t170 = sin(pkin(11));
	t203 = cos(pkin(6));
	t188 = t170 * t203;
	t202 = cos(pkin(11));
	t157 = -t172 * t188 + t202 * t173;
	t184 = t203 * t202;
	t153 = t170 * t172 - t173 * t184;
	t171 = sin(pkin(6));
	t192 = t171 * t173;
	t143 = atan2(-t153, -t192);
	t141 = sin(t143);
	t142 = cos(t143);
	t128 = -t141 * t153 - t142 * t192;
	t125 = 0.1e1 / t128;
	t169 = qJ(3) + qJ(4);
	t161 = sin(t169);
	t162 = cos(t169);
	t193 = t170 * t171;
	t140 = t157 * t162 + t161 * t193;
	t136 = 0.1e1 / t140;
	t166 = 0.1e1 / t173;
	t126 = 0.1e1 / t128 ^ 2;
	t137 = 0.1e1 / t140 ^ 2;
	t167 = 0.1e1 / t173 ^ 2;
	t155 = t170 * t173 + t172 * t184;
	t148 = t155 * qJD(2);
	t194 = t167 * t172;
	t189 = t153 * t194;
	t151 = t153 ^ 2;
	t164 = 0.1e1 / t171 ^ 2;
	t146 = t151 * t164 * t167 + 0.1e1;
	t144 = 0.1e1 / t146;
	t163 = 0.1e1 / t171;
	t196 = t144 * t163;
	t120 = (qJD(2) * t189 + t148 * t166) * t196;
	t182 = t141 * t192 - t142 * t153;
	t190 = t142 * t171 * t172;
	t117 = qJD(2) * t190 + t182 * t120 - t141 * t148;
	t201 = t117 * t125 * t126;
	t180 = -t202 * t172 - t173 * t188;
	t200 = t126 * t180;
	t139 = t157 * t161 - t162 * t193;
	t135 = t139 ^ 2;
	t132 = t135 * t137 + 0.1e1;
	t149 = t180 * qJD(2);
	t165 = qJD(3) + qJD(4);
	t185 = t165 * t193 + t149;
	t195 = t157 * t165;
	t133 = t185 * t161 + t162 * t195;
	t197 = t137 * t139;
	t134 = -t161 * t195 + t185 * t162;
	t198 = t134 * t136 * t137;
	t199 = 0.1e1 / t132 ^ 2 * (t133 * t197 - t135 * t198);
	t191 = -0.2e1 * t199;
	t183 = -t136 * t161 + t162 * t197;
	t181 = t155 * t166 + t189;
	t168 = t166 * t167;
	t152 = t180 ^ 2;
	t150 = t157 * qJD(2);
	t147 = t153 * qJD(2);
	t129 = 0.1e1 / t132;
	t124 = t126 * t152 + 0.1e1;
	t121 = t181 * t196;
	t118 = t182 * t121 - t141 * t155 + t190;
	t116 = (-0.2e1 * t181 / t146 ^ 2 * (qJD(2) * t151 * t168 * t172 + t148 * t153 * t167) * t164 + (t148 * t194 - t147 * t166 + (t155 * t194 + (0.2e1 * t168 * t172 ^ 2 + t166) * t153) * qJD(2)) * t144) * t163;
	t114 = t191 + 0.2e1 * (t129 * t133 * t137 + (-t129 * t198 - t137 * t199) * t139) * t139;
	t1 = [0, t116, 0, 0, 0, 0; 0, 0.2e1 * (-t118 * t200 - t125 * t157) / t124 ^ 2 * (-t150 * t200 - t152 * t201) + (t149 * t125 + (-t157 * t117 - t118 * t150) * t126 - (0.2e1 * t118 * t201 + (-(qJD(2) * t192 - t116 * t153 - t121 * t148 + (t121 * t192 - t155) * t120) * t142 - (t120 * t121 * t153 + t147 + (t116 * t173 + (-qJD(2) * t121 - t120) * t172) * t171) * t141) * t126) * t180) / t124, 0, 0, 0, 0; 0, -t183 * t180 * t191 + (t183 * t150 - ((-t136 * t165 - 0.2e1 * t139 * t198) * t162 + (t133 * t162 + (-t139 * t165 + t134) * t161) * t137) * t180) * t129, t114, t114, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:57
	% EndTime: 2019-10-09 23:03:58
	% DurationCPUTime: 1.72s
	% Computational Cost: add. (9025->114), mult. (13312->231), div. (822->12), fcn. (17117->13), ass. (0->113)
	t259 = sin(pkin(11));
	t261 = cos(pkin(11));
	t264 = sin(qJ(2));
	t262 = cos(pkin(6));
	t266 = cos(qJ(2));
	t293 = t262 * t266;
	t247 = -t259 * t264 + t261 * t293;
	t243 = t247 * qJD(2);
	t294 = t262 * t264;
	t248 = t259 * t266 + t261 * t294;
	t258 = qJ(3) + qJ(4);
	t255 = sin(t258);
	t257 = qJD(3) + qJD(4);
	t260 = sin(pkin(6));
	t297 = t260 * t261;
	t282 = t255 * t297;
	t256 = cos(t258);
	t299 = t256 * t257;
	t210 = t243 * t255 + t248 * t299 - t257 * t282;
	t232 = t248 * t255 + t256 * t297;
	t230 = t232 ^ 2;
	t296 = t260 * t264;
	t284 = t255 * t296;
	t240 = -t262 * t256 + t284;
	t238 = 0.1e1 / t240 ^ 2;
	t224 = t230 * t238 + 0.1e1;
	t222 = 0.1e1 / t224;
	t291 = qJD(2) * t266;
	t275 = t257 * t262 + t260 * t291;
	t283 = t256 * t296;
	t228 = t275 * t255 + t257 * t283;
	t237 = 0.1e1 / t240;
	t304 = t232 * t238;
	t194 = (-t210 * t237 + t228 * t304) * t222;
	t225 = atan2(-t232, t240);
	t220 = sin(t225);
	t221 = cos(t225);
	t278 = -t220 * t240 - t221 * t232;
	t190 = t278 * t194 - t220 * t210 + t221 * t228;
	t204 = -t220 * t232 + t221 * t240;
	t201 = 0.1e1 / t204;
	t202 = 0.1e1 / t204 ^ 2;
	t318 = t190 * t201 * t202;
	t285 = t259 * t294;
	t250 = t261 * t266 - t285;
	t298 = t259 * t260;
	t235 = t250 * t255 - t256 * t298;
	t317 = 0.2e1 * t235 * t318;
	t295 = t260 * t266;
	t274 = -t237 * t247 + t295 * t304;
	t316 = t255 * t274;
	t305 = t228 * t237 * t238;
	t315 = -0.2e1 * (t210 * t304 - t230 * t305) / t224 ^ 2;
	t236 = t250 * t256 + t255 * t298;
	t265 = cos(qJ(5));
	t249 = t259 * t293 + t261 * t264;
	t263 = sin(qJ(5));
	t302 = t249 * t263;
	t219 = t236 * t265 + t302;
	t215 = 0.1e1 / t219;
	t216 = 0.1e1 / t219 ^ 2;
	t245 = t249 * qJD(2);
	t280 = t257 * t298 - t245;
	t300 = t255 * t257;
	t213 = -t250 * t300 + t280 * t256;
	t246 = -qJD(2) * t285 + t261 * t291;
	t205 = t219 * qJD(5) + t213 * t263 - t246 * t265;
	t301 = t249 * t265;
	t218 = t236 * t263 - t301;
	t214 = t218 ^ 2;
	t209 = t214 * t216 + 0.1e1;
	t309 = t216 * t218;
	t290 = qJD(5) * t218;
	t206 = t213 * t265 + t246 * t263 - t290;
	t312 = t206 * t215 * t216;
	t314 = (t205 * t309 - t214 * t312) / t209 ^ 2;
	t313 = t202 * t235;
	t212 = t250 * t299 + t280 * t255;
	t311 = t212 * t202;
	t310 = t215 * t263;
	t308 = t218 * t265;
	t307 = t220 * t235;
	t306 = t221 * t235;
	t303 = t249 * t255;
	t292 = qJD(2) * t264;
	t231 = t235 ^ 2;
	t200 = t231 * t202 + 0.1e1;
	t289 = 0.2e1 * (-t231 * t318 + t235 * t311) / t200 ^ 2;
	t288 = -0.2e1 * t314;
	t286 = t218 * t312;
	t281 = -0.2e1 * t232 * t305;
	t279 = qJD(5) * t249 * t256 - t245;
	t277 = t216 * t308 - t310;
	t234 = t248 * t256 - t282;
	t241 = t262 * t255 + t283;
	t276 = -t234 * t237 + t241 * t304;
	t273 = qJD(5) * t250 - t246 * t256 + t249 * t300;
	t244 = t248 * qJD(2);
	t229 = t275 * t256 - t257 * t284;
	t227 = t250 * t263 - t256 * t301;
	t226 = -t250 * t265 - t256 * t302;
	t211 = -t248 * t300 + (-t257 * t297 + t243) * t256;
	t207 = 0.1e1 / t209;
	t198 = 0.1e1 / t200;
	t196 = t222 * t316;
	t195 = t276 * t222;
	t192 = (-t220 * t247 + t221 * t295) * t255 + t278 * t196;
	t191 = t278 * t195 - t220 * t234 + t221 * t241;
	t189 = t276 * t315 + (t241 * t281 - t211 * t237 + (t210 * t241 + t228 * t234 + t229 * t232) * t238) * t222;
	t187 = t315 * t316 + (t274 * t299 + (t281 * t295 + t237 * t244 + (t228 * t247 + (t210 * t266 - t232 * t292) * t260) * t238) * t255) * t222;
	t186 = t277 * t235 * t288 + (t277 * t212 + ((-qJD(5) * t215 - 0.2e1 * t286) * t265 + (t205 * t265 + (t206 - t290) * t263) * t216) * t235) * t207;
	t185 = (t191 * t313 - t201 * t236) * t289 + (t191 * t317 + t213 * t201 + (-t236 * t190 - t191 * t212 - (-t189 * t232 - t195 * t210 + t229 + (-t195 * t240 - t234) * t194) * t306 - (-t189 * t240 - t195 * t228 - t211 + (t195 * t232 - t241) * t194) * t307) * t202) * t198;
	t1 = [0, t187, t189, t189, 0, 0; 0, (t192 * t313 + t201 * t303) * t289 + ((-t246 * t255 - t249 * t299) * t201 + (-t311 + t317) * t192 + (t303 * t190 - (-t187 * t232 - t196 * t210 + (-t255 * t292 + t266 * t299) * t260 + (-t196 * t240 - t247 * t255) * t194) * t306 - (-t247 * t299 - t187 * t240 - t196 * t228 + t244 * t255 + (t196 * t232 - t255 * t295) * t194) * t307) * t202) * t198, t185, t185, 0, 0; 0, 0.2e1 * (-t215 * t226 + t227 * t309) * t314 + (0.2e1 * t227 * t286 - t279 * t215 * t265 + t273 * t310 + (-t279 * t218 * t263 - t227 * t205 - t226 * t206 - t273 * t308) * t216) * t207, t186, t186, t288 + 0.2e1 * (t205 * t216 * t207 + (-t207 * t312 - t216 * t314) * t218) * t218, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:57
	% EndTime: 2019-10-09 23:04:00
	% DurationCPUTime: 3.09s
	% Computational Cost: add. (16187->161), mult. (28326->308), div. (1062->12), fcn. (36524->13), ass. (0->132)
	t311 = qJ(3) + qJ(4);
	t309 = cos(t311);
	t314 = sin(qJ(2));
	t316 = cos(qJ(2));
	t376 = cos(pkin(11));
	t377 = cos(pkin(6));
	t338 = t377 * t376;
	t375 = sin(pkin(11));
	t325 = -t314 * t338 - t375 * t316;
	t336 = t375 * t314 - t316 * t338;
	t300 = t336 * qJD(2);
	t310 = qJD(3) + qJD(4);
	t312 = sin(pkin(6));
	t343 = t312 * t376;
	t333 = t310 * t343 + t300;
	t308 = sin(t311);
	t360 = t310 * t308;
	t265 = -t309 * t333 + t325 * t360;
	t285 = -t308 * t343 - t309 * t325;
	t315 = cos(qJ(5));
	t313 = sin(qJ(5));
	t328 = t336 * t313;
	t275 = t285 * t315 + t328;
	t323 = t325 * qJD(2);
	t247 = qJD(5) * t275 + t265 * t313 + t315 * t323;
	t327 = t336 * t315;
	t273 = t285 * t313 - t327;
	t268 = t273 ^ 2;
	t359 = t312 * t314;
	t345 = t309 * t359;
	t296 = t377 * t308 + t345;
	t357 = t315 * t316;
	t291 = t296 * t313 + t312 * t357;
	t289 = 0.1e1 / t291 ^ 2;
	t259 = t268 * t289 + 0.1e1;
	t257 = 0.1e1 / t259;
	t355 = qJD(2) * t312;
	t326 = t377 * t310 + t316 * t355;
	t346 = t308 * t359;
	t282 = t309 * t326 - t310 * t346;
	t358 = t313 * t316;
	t292 = t296 * t315 - t312 * t358;
	t344 = t314 * t355;
	t261 = qJD(5) * t292 + t282 * t313 - t315 * t344;
	t288 = 0.1e1 / t291;
	t366 = t273 * t289;
	t234 = (-t247 * t288 + t261 * t366) * t257;
	t260 = atan2(-t273, t291);
	t254 = sin(t260);
	t255 = cos(t260);
	t335 = -t254 * t291 - t255 * t273;
	t230 = t234 * t335 - t254 * t247 + t255 * t261;
	t246 = -t254 * t273 + t255 * t291;
	t243 = 0.1e1 / t246;
	t244 = 0.1e1 / t246 ^ 2;
	t381 = t230 * t243 * t244;
	t337 = t377 * t375;
	t305 = -t314 * t337 + t376 * t316;
	t342 = t312 * t375;
	t287 = t305 * t309 + t308 * t342;
	t324 = -t376 * t314 - t316 * t337;
	t276 = t287 * t313 + t315 * t324;
	t380 = -0.2e1 * t276;
	t339 = 0.2e1 * t276 * t381;
	t284 = t308 * t325 - t309 * t343;
	t295 = t377 * t309 - t346;
	t329 = -t284 * t288 + t295 * t366;
	t379 = t313 * t329;
	t368 = t261 * t288 * t289;
	t378 = -0.2e1 * (t247 * t366 - t268 * t368) / t259 ^ 2;
	t363 = t324 * t313;
	t277 = t287 * t315 - t363;
	t270 = 0.1e1 / t277;
	t271 = 0.1e1 / t277 ^ 2;
	t286 = -t305 * t308 + t309 * t342;
	t283 = t286 ^ 2;
	t365 = t283 * t271;
	t256 = 0.1e1 + t365;
	t301 = t324 * qJD(2);
	t332 = t310 * t342 + t301;
	t362 = t309 * t310;
	t266 = -t305 * t362 - t308 * t332;
	t267 = -t305 * t360 + t309 * t332;
	t302 = t305 * qJD(2);
	t250 = -qJD(5) * t276 + t267 * t315 + t302 * t313;
	t371 = t250 * t270 * t271;
	t348 = t283 * t371;
	t367 = t271 * t286;
	t374 = (t266 * t367 - t348) / t256 ^ 2;
	t373 = t244 * t276;
	t249 = qJD(5) * t277 + t267 * t313 - t302 * t315;
	t372 = t249 * t244;
	t370 = t254 * t276;
	t369 = t255 * t276;
	t364 = t286 * t313;
	t361 = t309 * t315;
	t356 = qJD(2) * t309;
	t354 = qJD(5) * t313;
	t353 = qJD(5) * t315;
	t352 = t309 * qJD(5);
	t269 = t276 ^ 2;
	t242 = t269 * t244 + 0.1e1;
	t351 = 0.2e1 * (-t269 * t381 + t276 * t372) / t242 ^ 2;
	t350 = 0.2e1 * t374;
	t347 = t286 * t371;
	t340 = -0.2e1 * t273 * t368;
	t331 = -t275 * t288 + t292 * t366;
	t278 = -t309 * t328 + t315 * t325;
	t293 = (t309 * t358 - t314 * t315) * t312;
	t330 = -t278 * t288 + t293 * t366;
	t281 = -t308 * t326 - t310 * t345;
	t280 = t305 * t313 + t324 * t361;
	t279 = -t305 * t315 + t309 * t363;
	t264 = t308 * t333 + t325 * t362;
	t263 = ((-qJD(2) + t352) * t357 + (-t316 * t360 + (qJD(5) - t356) * t314) * t313) * t312;
	t262 = -qJD(5) * t291 + t282 * t315 + t313 * t344;
	t252 = 0.1e1 / t256;
	t251 = t328 * t360 + (-t336 * t352 + t300) * t315 + (t313 * t356 - t354) * t325;
	t248 = qJD(5) * t327 + t265 * t315 - t285 * t354 - t313 * t323;
	t240 = 0.1e1 / t242;
	t239 = t257 * t379;
	t238 = t330 * t257;
	t237 = t331 * t257;
	t233 = (-t254 * t284 + t255 * t295) * t313 + t335 * t239;
	t232 = t238 * t335 - t254 * t278 + t255 * t293;
	t231 = t237 * t335 - t254 * t275 + t255 * t292;
	t229 = (t270 * t287 + t315 * t365) * t350 + (0.2e1 * t315 * t348 - t267 * t270 + (-0.2e1 * t266 * t286 * t315 + t250 * t287 + t283 * t354) * t271) * t252;
	t228 = t330 * t378 + (t293 * t340 - t251 * t288 + (t247 * t293 + t261 * t278 + t263 * t273) * t289) * t257;
	t226 = t331 * t378 + (t292 * t340 - t248 * t288 + (t247 * t292 + t261 * t275 + t262 * t273) * t289) * t257;
	t225 = t378 * t379 + (t329 * t353 + (t295 * t340 - t264 * t288 + (t247 * t295 + t261 * t284 + t273 * t281) * t289) * t313) * t257;
	t224 = (t233 * t373 - t243 * t364) * t351 + ((t266 * t313 + t286 * t353) * t243 + (-t372 + t339) * t233 + (-t364 * t230 - (t295 * t353 - t225 * t273 - t239 * t247 + t281 * t313 + (-t239 * t291 - t284 * t313) * t234) * t369 - (-t284 * t353 - t225 * t291 - t239 * t261 - t264 * t313 + (t239 * t273 - t295 * t313) * t234) * t370) * t244) * t240;
	t1 = [0, t228, t225, t225, t226, 0; 0, (t232 * t373 - t243 * t279) * t351 + (t232 * t339 + (-t279 * t230 - t232 * t249 - (-t228 * t273 - t238 * t247 + t263 + (-t238 * t291 - t278) * t234) * t369 - (-t228 * t291 - t238 * t261 - t251 + (t238 * t273 - t293) * t234) * t370) * t244 + ((t324 * t352 - t301) * t315 + (qJD(5) * t305 - t302 * t309 - t324 * t360) * t313) * t243) * t240, t224, t224, (t231 * t373 - t243 * t277) * t351 + (t231 * t339 + t250 * t243 + (-t277 * t230 - t231 * t249 - (-t226 * t273 - t237 * t247 + t262 + (-t237 * t291 - t275) * t234) * t369 - (-t226 * t291 - t237 * t261 - t248 + (t237 * t273 - t292) * t234) * t370) * t244) * t240, 0; 0, (t270 * t308 * t324 + t280 * t367) * t350 + (0.2e1 * t280 * t347 + (t302 * t308 - t324 * t362) * t270 + (-(t301 * t313 - t302 * t361 + t305 * t353) * t286 - t280 * t266 - (-t308 * t250 - (t313 * t352 + t315 * t360) * t286) * t324) * t271) * t252, t229, t229, t367 * t374 * t380 + (t347 * t380 + (t249 * t286 + t266 * t276) * t271) * t252, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end