% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR2
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
%   Wie in S6PRRRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
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
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.58s
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
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:36
	% DurationCPUTime: 1.59s
	% Computational Cost: add. (8455->106), mult. (12233->219), div. (792->12), fcn. (15777->13), ass. (0->109)
	t235 = sin(pkin(11));
	t238 = cos(pkin(11));
	t240 = sin(qJ(2));
	t239 = cos(pkin(6));
	t241 = cos(qJ(2));
	t265 = t239 * t241;
	t222 = -t235 * t240 + t238 * t265;
	t218 = t222 * qJD(2);
	t266 = t239 * t240;
	t223 = t235 * t241 + t238 * t266;
	t233 = qJ(3) + qJ(4);
	t230 = sin(t233);
	t232 = qJD(3) + qJD(4);
	t236 = sin(pkin(6));
	t269 = t236 * t238;
	t255 = t230 * t269;
	t231 = cos(t233);
	t271 = t231 * t232;
	t186 = t218 * t230 + t223 * t271 - t232 * t255;
	t208 = t223 * t230 + t231 * t269;
	t206 = t208 ^ 2;
	t268 = t236 * t240;
	t257 = t230 * t268;
	t216 = -t239 * t231 + t257;
	t214 = 0.1e1 / t216 ^ 2;
	t200 = t206 * t214 + 0.1e1;
	t198 = 0.1e1 / t200;
	t263 = qJD(2) * t241;
	t249 = t232 * t239 + t236 * t263;
	t256 = t231 * t268;
	t204 = t249 * t230 + t232 * t256;
	t213 = 0.1e1 / t216;
	t276 = t208 * t214;
	t170 = (-t186 * t213 + t204 * t276) * t198;
	t201 = atan2(-t208, t216);
	t196 = sin(t201);
	t197 = cos(t201);
	t252 = -t196 * t216 - t197 * t208;
	t166 = t252 * t170 - t196 * t186 + t197 * t204;
	t180 = -t196 * t208 + t197 * t216;
	t177 = 0.1e1 / t180;
	t178 = 0.1e1 / t180 ^ 2;
	t289 = t166 * t177 * t178;
	t258 = t235 * t266;
	t225 = t238 * t241 - t258;
	t270 = t235 * t236;
	t211 = t225 * t230 - t231 * t270;
	t288 = 0.2e1 * t211 * t289;
	t267 = t236 * t241;
	t248 = -t213 * t222 + t267 * t276;
	t287 = t230 * t248;
	t277 = t204 * t213 * t214;
	t286 = -0.2e1 * (t186 * t276 - t206 * t277) / t200 ^ 2;
	t212 = t225 * t231 + t230 * t270;
	t237 = cos(pkin(12));
	t224 = t235 * t265 + t238 * t240;
	t234 = sin(pkin(12));
	t274 = t224 * t234;
	t195 = t212 * t237 + t274;
	t191 = 0.1e1 / t195;
	t192 = 0.1e1 / t195 ^ 2;
	t285 = t178 * t211;
	t220 = t224 * qJD(2);
	t253 = t232 * t270 - t220;
	t272 = t230 * t232;
	t189 = -t225 * t272 + t253 * t231;
	t221 = -qJD(2) * t258 + t238 * t263;
	t185 = t189 * t237 + t221 * t234;
	t284 = t185 * t191 * t192;
	t188 = t225 * t271 + t253 * t230;
	t283 = t188 * t178;
	t282 = t191 * t234;
	t273 = t224 * t237;
	t194 = t212 * t234 - t273;
	t281 = t192 * t194;
	t280 = t194 * t237;
	t279 = t196 * t211;
	t278 = t197 * t211;
	t275 = t224 * t230;
	t264 = qJD(2) * t240;
	t207 = t211 ^ 2;
	t176 = t207 * t178 + 0.1e1;
	t262 = 0.2e1 * (-t207 * t289 + t211 * t283) / t176 ^ 2;
	t190 = t194 ^ 2;
	t183 = t190 * t192 + 0.1e1;
	t184 = t189 * t234 - t221 * t237;
	t261 = 0.2e1 * (t184 * t281 - t190 * t284) / t183 ^ 2;
	t259 = t194 * t284;
	t254 = -0.2e1 * t208 * t277;
	t210 = t223 * t231 - t255;
	t217 = t239 * t230 + t256;
	t251 = -t210 * t213 + t217 * t276;
	t250 = -t221 * t231 + t224 * t272;
	t219 = t223 * qJD(2);
	t205 = t249 * t231 - t232 * t257;
	t203 = t225 * t234 - t231 * t273;
	t202 = -t225 * t237 - t231 * t274;
	t187 = -t223 * t272 + (-t232 * t269 + t218) * t231;
	t181 = 0.1e1 / t183;
	t174 = 0.1e1 / t176;
	t172 = t198 * t287;
	t171 = t251 * t198;
	t168 = (-t196 * t222 + t197 * t267) * t230 + t252 * t172;
	t167 = t252 * t171 - t196 * t210 + t197 * t217;
	t165 = t251 * t286 + (t217 * t254 - t187 * t213 + (t186 * t217 + t204 * t210 + t205 * t208) * t214) * t198;
	t163 = t286 * t287 + (t248 * t271 + (t254 * t267 + t213 * t219 + (t204 * t222 + (t186 * t241 - t208 * t264) * t236) * t214) * t230) * t198;
	t162 = (-t192 * t280 + t282) * t211 * t261 + (-0.2e1 * t211 * t237 * t259 - t188 * t282 + (t188 * t280 + (t184 * t237 + t185 * t234) * t211) * t192) * t181;
	t161 = (t167 * t285 - t177 * t212) * t262 + (t167 * t288 + t189 * t177 + (-t212 * t166 - t167 * t188 - (-t165 * t208 - t171 * t186 + t205 + (-t171 * t216 - t210) * t170) * t278 - (-t165 * t216 - t171 * t204 - t187 + (t171 * t208 - t217) * t170) * t279) * t178) * t174;
	t1 = [0, t163, t165, t165, 0, 0; 0, (t168 * t285 + t177 * t275) * t262 + ((-t221 * t230 - t224 * t271) * t177 + (-t283 + t288) * t168 + (t275 * t166 - (-t163 * t208 - t172 * t186 + (-t230 * t264 + t241 * t271) * t236 + (-t172 * t216 - t222 * t230) * t170) * t278 - (-t222 * t271 - t163 * t216 - t172 * t204 + t219 * t230 + (t172 * t208 - t230 * t267) * t170) * t279) * t178) * t174, t161, t161, 0, 0; 0, (-t191 * t202 + t203 * t281) * t261 + ((t220 * t237 + t250 * t234) * t191 + 0.2e1 * t203 * t259 + (-t202 * t185 - (-t220 * t234 + t250 * t237) * t194 - t203 * t184) * t192) * t181, t162, t162, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:36
	% DurationCPUTime: 1.72s
	% Computational Cost: add. (9451->115), mult. (13312->231), div. (822->12), fcn. (17117->13), ass. (0->113)
	t270 = sin(pkin(11));
	t272 = cos(pkin(11));
	t274 = sin(qJ(2));
	t273 = cos(pkin(6));
	t275 = cos(qJ(2));
	t302 = t273 * t275;
	t255 = -t270 * t274 + t272 * t302;
	t251 = t255 * qJD(2);
	t303 = t273 * t274;
	t256 = t270 * t275 + t272 * t303;
	t269 = qJ(3) + qJ(4);
	t265 = sin(t269);
	t268 = qJD(3) + qJD(4);
	t271 = sin(pkin(6));
	t306 = t271 * t272;
	t291 = t265 * t306;
	t266 = cos(t269);
	t308 = t266 * t268;
	t218 = t251 * t265 + t256 * t308 - t268 * t291;
	t240 = t256 * t265 + t266 * t306;
	t238 = t240 ^ 2;
	t305 = t271 * t274;
	t293 = t265 * t305;
	t249 = -t266 * t273 + t293;
	t247 = 0.1e1 / t249 ^ 2;
	t232 = t238 * t247 + 0.1e1;
	t230 = 0.1e1 / t232;
	t300 = qJD(2) * t275;
	t284 = t268 * t273 + t271 * t300;
	t292 = t266 * t305;
	t236 = t265 * t284 + t268 * t292;
	t246 = 0.1e1 / t249;
	t312 = t240 * t247;
	t202 = (-t218 * t246 + t236 * t312) * t230;
	t233 = atan2(-t240, t249);
	t228 = sin(t233);
	t229 = cos(t233);
	t287 = -t228 * t249 - t229 * t240;
	t198 = t202 * t287 - t228 * t218 + t229 * t236;
	t212 = -t228 * t240 + t229 * t249;
	t209 = 0.1e1 / t212;
	t210 = 0.1e1 / t212 ^ 2;
	t326 = t198 * t209 * t210;
	t294 = t270 * t303;
	t258 = t272 * t275 - t294;
	t307 = t270 * t271;
	t243 = t258 * t265 - t266 * t307;
	t325 = 0.2e1 * t243 * t326;
	t304 = t271 * t275;
	t283 = -t246 * t255 + t304 * t312;
	t324 = t265 * t283;
	t313 = t236 * t246 * t247;
	t323 = -0.2e1 * (t218 * t312 - t238 * t313) / t232 ^ 2;
	t244 = t258 * t266 + t265 * t307;
	t257 = t270 * t302 + t272 * t274;
	t267 = pkin(12) + qJ(6);
	t263 = sin(t267);
	t264 = cos(t267);
	t227 = t244 * t264 + t257 * t263;
	t223 = 0.1e1 / t227;
	t224 = 0.1e1 / t227 ^ 2;
	t253 = t257 * qJD(2);
	t289 = t268 * t307 - t253;
	t309 = t265 * t268;
	t221 = -t258 * t309 + t266 * t289;
	t254 = -qJD(2) * t294 + t272 * t300;
	t213 = qJD(6) * t227 + t221 * t263 - t254 * t264;
	t226 = t244 * t263 - t257 * t264;
	t222 = t226 ^ 2;
	t217 = t222 * t224 + 0.1e1;
	t317 = t224 * t226;
	t299 = qJD(6) * t226;
	t214 = t221 * t264 + t254 * t263 - t299;
	t320 = t214 * t223 * t224;
	t322 = (t213 * t317 - t222 * t320) / t217 ^ 2;
	t321 = t210 * t243;
	t220 = t258 * t308 + t265 * t289;
	t319 = t220 * t210;
	t318 = t223 * t263;
	t316 = t226 * t264;
	t315 = t228 * t243;
	t314 = t229 * t243;
	t311 = t257 * t265;
	t310 = t257 * t266;
	t301 = qJD(2) * t274;
	t239 = t243 ^ 2;
	t208 = t210 * t239 + 0.1e1;
	t298 = 0.2e1 * (-t239 * t326 + t243 * t319) / t208 ^ 2;
	t297 = -0.2e1 * t322;
	t295 = t226 * t320;
	t290 = -0.2e1 * t240 * t313;
	t288 = qJD(6) * t310 - t253;
	t286 = t224 * t316 - t318;
	t242 = t256 * t266 - t291;
	t250 = t265 * t273 + t292;
	t285 = -t242 * t246 + t250 * t312;
	t282 = qJD(6) * t258 - t254 * t266 + t257 * t309;
	t252 = t256 * qJD(2);
	t237 = t266 * t284 - t268 * t293;
	t235 = t258 * t263 - t264 * t310;
	t234 = -t258 * t264 - t263 * t310;
	t219 = -t256 * t309 + (-t268 * t306 + t251) * t266;
	t215 = 0.1e1 / t217;
	t206 = 0.1e1 / t208;
	t204 = t230 * t324;
	t203 = t285 * t230;
	t200 = (-t228 * t255 + t229 * t304) * t265 + t287 * t204;
	t199 = t203 * t287 - t228 * t242 + t229 * t250;
	t197 = t285 * t323 + (t250 * t290 - t219 * t246 + (t218 * t250 + t236 * t242 + t237 * t240) * t247) * t230;
	t195 = t323 * t324 + (t283 * t308 + (t290 * t304 + t246 * t252 + (t236 * t255 + (t218 * t275 - t240 * t301) * t271) * t247) * t265) * t230;
	t194 = t286 * t243 * t297 + (t286 * t220 + ((-qJD(6) * t223 - 0.2e1 * t295) * t264 + (t213 * t264 + (t214 - t299) * t263) * t224) * t243) * t215;
	t193 = (t199 * t321 - t209 * t244) * t298 + (t199 * t325 + t221 * t209 + (-t244 * t198 - t199 * t220 - (-t197 * t240 - t203 * t218 + t237 + (-t203 * t249 - t242) * t202) * t314 - (-t197 * t249 - t203 * t236 - t219 + (t203 * t240 - t250) * t202) * t315) * t210) * t206;
	t1 = [0, t195, t197, t197, 0, 0; 0, (t200 * t321 + t209 * t311) * t298 + ((-t254 * t265 - t257 * t308) * t209 + (-t319 + t325) * t200 + (t311 * t198 - (-t195 * t240 - t204 * t218 + (-t265 * t301 + t275 * t308) * t271 + (-t204 * t249 - t255 * t265) * t202) * t314 - (-t255 * t308 - t195 * t249 - t204 * t236 + t252 * t265 + (t204 * t240 - t265 * t304) * t202) * t315) * t210) * t206, t193, t193, 0, 0; 0, 0.2e1 * (-t223 * t234 + t235 * t317) * t322 + (0.2e1 * t235 * t295 - t288 * t223 * t264 + t282 * t318 + (-t226 * t263 * t288 - t235 * t213 - t234 * t214 - t282 * t316) * t224) * t215, t194, t194, 0, t297 + 0.2e1 * (t213 * t224 * t215 + (-t215 * t320 - t224 * t322) * t226) * t226;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end