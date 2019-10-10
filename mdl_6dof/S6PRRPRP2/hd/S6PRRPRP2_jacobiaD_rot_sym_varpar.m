% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRP2
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
%   Wie in S6PRRPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
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
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
	% DurationCPUTime: 0.57s
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
	t148 = -t139 * t165 - t141 * t155;
	t119 = t148 * qJD(2);
	t103 = qJD(3) * t110 + t119 * t138;
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
	t88 = (t137 * t139 - t168) * t112 + (t159 * t91 - t125) * t111;
	t87 = (-t123 * t90 + t137 * t158) * t112 + (t159 * t90 - t118) * t111;
	t86 = (-0.2e1 * t150 * (t118 * t123 * t134 + t121 * t135 * t158) * t132 / t116 ^ 2 + (t118 * t161 - t117 * t133 + (t125 * t161 + (0.2e1 * t135 * t139 ^ 2 + t133) * t123) * qJD(2)) * t114) * t131;
	t1 = [0, t86, 0, 0, 0, 0; 0, 0.2e1 * (-t127 * t95 - t167 * t88) / t94 ^ 2 * (-t122 * t87 * t97 - t120 * t167) + (-t88 * t120 * t96 + t119 * t95 + (-0.2e1 * t148 * t88 * t97 - t127 * t96) * t87 - (-(-t118 * t91 - t123 * t86 - t125 * t90 + (t90 * t91 + qJD(2)) * t159) * t112 - (t90 * t168 + t117 + (t141 * t86 + (-qJD(2) * t91 - t90) * t139) * t137) * t111) * t167) / t94, 0, 0, 0, 0; 0, -t151 * t148 * t157 + (t151 * t120 - ((-qJD(3) * t106 + 0.2e1 * t149 * t164) * t140 + (t103 * t140 + (t104 + t170) * t138) * t107) * t148) * t100, t157 - 0.2e1 * (t100 * t103 * t107 - (-t100 * t164 - t107 * t169) * t149) * t149, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (941->56), mult. (2271->131), div. (423->14), fcn. (2956->11), ass. (0->66)
	t149 = sin(qJ(2));
	t150 = cos(qJ(2));
	t147 = sin(pkin(10));
	t176 = cos(pkin(6));
	t164 = t147 * t176;
	t175 = cos(pkin(10));
	t135 = -t149 * t164 + t175 * t150;
	t143 = qJ(3) + pkin(11);
	t139 = sin(t143);
	t140 = cos(t143);
	t148 = sin(pkin(6));
	t169 = t147 * t148;
	t158 = -t135 * t139 + t140 * t169;
	t180 = t158 * qJD(3);
	t161 = t176 * t175;
	t131 = t147 * t149 - t150 * t161;
	t168 = t148 * t150;
	t121 = atan2(-t131, -t168);
	t119 = sin(t121);
	t120 = cos(t121);
	t106 = -t119 * t131 - t120 * t168;
	t103 = 0.1e1 / t106;
	t118 = t135 * t140 + t139 * t169;
	t114 = 0.1e1 / t118;
	t144 = 0.1e1 / t150;
	t104 = 0.1e1 / t106 ^ 2;
	t115 = 0.1e1 / t118 ^ 2;
	t145 = 0.1e1 / t150 ^ 2;
	t133 = t147 * t150 + t149 * t161;
	t126 = t133 * qJD(2);
	t167 = qJD(2) * t149;
	t170 = t145 * t149;
	t165 = t131 * t170;
	t129 = t131 ^ 2;
	t142 = 0.1e1 / t148 ^ 2;
	t124 = t129 * t142 * t145 + 0.1e1;
	t122 = 0.1e1 / t124;
	t141 = 0.1e1 / t148;
	t171 = t122 * t141;
	t98 = (qJD(2) * t165 + t126 * t144) * t171;
	t95 = (-t131 * t98 + t148 * t167) * t120 + (t98 * t168 - t126) * t119;
	t179 = t103 * t104 * t95;
	t113 = t158 ^ 2;
	t109 = t113 * t115 + 0.1e1;
	t157 = -t175 * t149 - t150 * t164;
	t127 = t157 * qJD(2);
	t111 = t118 * qJD(3) + t127 * t139;
	t172 = t115 * t158;
	t112 = t127 * t140 + t180;
	t173 = t112 * t114 * t115;
	t178 = 0.1e1 / t109 ^ 2 * (-t111 * t172 - t113 * t173);
	t159 = t133 * t144 + t165;
	t99 = t159 * t171;
	t177 = t131 * t99;
	t174 = t104 * t157;
	t166 = -0.2e1 * t178;
	t160 = -t114 * t139 - t140 * t172;
	t146 = t144 * t145;
	t130 = t157 ^ 2;
	t128 = t135 * qJD(2);
	t125 = t131 * qJD(2);
	t107 = 0.1e1 / t109;
	t102 = t130 * t104 + 0.1e1;
	t96 = (t148 * t149 - t177) * t120 + (t99 * t168 - t133) * t119;
	t94 = (-0.2e1 * t159 / t124 ^ 2 * (t126 * t131 * t145 + t129 * t146 * t167) * t142 + (t126 * t170 - t125 * t144 + (t133 * t170 + (0.2e1 * t146 * t149 ^ 2 + t144) * t131) * qJD(2)) * t122) * t141;
	t1 = [0, t94, 0, 0, 0, 0; 0, 0.2e1 * (-t103 * t135 - t96 * t174) * (-t128 * t174 - t130 * t179) / t102 ^ 2 + (t127 * t103 + (-t96 * t128 - t135 * t95) * t104 - (0.2e1 * t96 * t179 + (-(-t126 * t99 - t131 * t94 - t133 * t98 + (t98 * t99 + qJD(2)) * t168) * t120 - (t98 * t177 + t125 + (t150 * t94 + (-qJD(2) * t99 - t98) * t149) * t148) * t119) * t104) * t157) / t102, 0, 0, 0, 0; 0, -t160 * t157 * t166 + (t160 * t128 - ((-qJD(3) * t114 + 0.2e1 * t158 * t173) * t140 + (t111 * t140 + (t112 + t180) * t139) * t115) * t157) * t107, t166 - 0.2e1 * (t107 * t111 * t115 - (-t107 * t173 - t115 * t178) * t158) * t158, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:16
	% DurationCPUTime: 1.40s
	% Computational Cost: add. (5804->110), mult. (9085->227), div. (559->12), fcn. (11668->13), ass. (0->106)
	t224 = sin(pkin(10));
	t226 = cos(pkin(10));
	t229 = sin(qJ(2));
	t227 = cos(pkin(6));
	t231 = cos(qJ(2));
	t257 = t227 * t231;
	t213 = -t224 * t229 + t226 * t257;
	t209 = t213 * qJD(2);
	t258 = t227 * t229;
	t214 = t224 * t231 + t226 * t258;
	t223 = qJ(3) + pkin(11);
	t221 = sin(t223);
	t225 = sin(pkin(6));
	t261 = t225 * t226;
	t247 = t221 * t261;
	t222 = cos(t223);
	t254 = qJD(3) * t222;
	t176 = -qJD(3) * t247 + t209 * t221 + t214 * t254;
	t198 = t214 * t221 + t222 * t261;
	t196 = t198 ^ 2;
	t260 = t225 * t229;
	t206 = t221 * t260 - t227 * t222;
	t204 = 0.1e1 / t206 ^ 2;
	t190 = t196 * t204 + 0.1e1;
	t188 = 0.1e1 / t190;
	t207 = t227 * t221 + t222 * t260;
	t255 = qJD(2) * t231;
	t246 = t225 * t255;
	t194 = t207 * qJD(3) + t221 * t246;
	t203 = 0.1e1 / t206;
	t266 = t198 * t204;
	t160 = (-t176 * t203 + t194 * t266) * t188;
	t191 = atan2(-t198, t206);
	t184 = sin(t191);
	t185 = cos(t191);
	t243 = -t184 * t206 - t185 * t198;
	t156 = t243 * t160 - t184 * t176 + t185 * t194;
	t170 = -t184 * t198 + t185 * t206;
	t167 = 0.1e1 / t170;
	t168 = 0.1e1 / t170 ^ 2;
	t280 = t156 * t167 * t168;
	t248 = t224 * t258;
	t216 = t226 * t231 - t248;
	t262 = t224 * t225;
	t240 = -t216 * t221 + t222 * t262;
	t279 = -0.2e1 * t240 * t280;
	t259 = t225 * t231;
	t239 = -t203 * t213 + t259 * t266;
	t278 = t221 * t239;
	t267 = t194 * t203 * t204;
	t277 = -0.2e1 * (t176 * t266 - t196 * t267) / t190 ^ 2;
	t202 = t216 * t222 + t221 * t262;
	t230 = cos(qJ(5));
	t215 = t224 * t257 + t226 * t229;
	t228 = sin(qJ(5));
	t264 = t215 * t228;
	t187 = t202 * t230 + t264;
	t181 = 0.1e1 / t187;
	t182 = 0.1e1 / t187 ^ 2;
	t211 = t215 * qJD(2);
	t179 = t240 * qJD(3) - t211 * t222;
	t212 = -qJD(2) * t248 + t226 * t255;
	t171 = t187 * qJD(5) + t179 * t228 - t212 * t230;
	t263 = t215 * t230;
	t186 = t202 * t228 - t263;
	t180 = t186 ^ 2;
	t175 = t180 * t182 + 0.1e1;
	t271 = t182 * t186;
	t253 = qJD(5) * t186;
	t172 = t179 * t230 + t212 * t228 - t253;
	t274 = t172 * t181 * t182;
	t276 = (t171 * t271 - t180 * t274) / t175 ^ 2;
	t275 = t168 * t240;
	t178 = t202 * qJD(3) - t211 * t221;
	t273 = t178 * t168;
	t272 = t181 * t228;
	t270 = t184 * t240;
	t269 = t185 * t240;
	t268 = t186 * t230;
	t265 = t215 * t221;
	t256 = qJD(2) * t229;
	t197 = t240 ^ 2;
	t166 = t197 * t168 + 0.1e1;
	t252 = 0.2e1 * (-t197 * t280 - t240 * t273) / t166 ^ 2;
	t251 = -0.2e1 * t276;
	t249 = t186 * t274;
	t245 = -0.2e1 * t198 * t267;
	t244 = qJD(5) * t215 * t222 - t211;
	t242 = t182 * t268 - t272;
	t200 = t214 * t222 - t247;
	t241 = -t200 * t203 + t207 * t266;
	t238 = qJD(3) * t265 + qJD(5) * t216 - t212 * t222;
	t210 = t214 * qJD(2);
	t195 = -t206 * qJD(3) + t222 * t246;
	t193 = t216 * t228 - t222 * t263;
	t192 = -t216 * t230 - t222 * t264;
	t177 = -t198 * qJD(3) + t209 * t222;
	t173 = 0.1e1 / t175;
	t163 = 0.1e1 / t166;
	t162 = t188 * t278;
	t161 = t241 * t188;
	t158 = (-t184 * t213 + t185 * t259) * t221 + t243 * t162;
	t157 = t243 * t161 - t184 * t200 + t185 * t207;
	t155 = t241 * t277 + (t207 * t245 - t177 * t203 + (t176 * t207 + t194 * t200 + t195 * t198) * t204) * t188;
	t153 = t277 * t278 + (t239 * t254 + (t245 * t259 + t203 * t210 + (t194 * t213 + (t176 * t231 - t198 * t256) * t225) * t204) * t221) * t188;
	t1 = [0, t153, t155, 0, 0, 0; 0, (-t158 * t275 + t167 * t265) * t252 + ((-t212 * t221 - t215 * t254) * t167 + (-t273 + t279) * t158 + (t265 * t156 + (-t153 * t198 - t162 * t176 + (-t221 * t256 + t231 * t254) * t225 + (-t162 * t206 - t213 * t221) * t160) * t269 + (-t213 * t254 - t153 * t206 - t162 * t194 + t210 * t221 + (t162 * t198 - t221 * t259) * t160) * t270) * t168) * t163, (-t157 * t275 - t167 * t202) * t252 + (t157 * t279 + t179 * t167 + (-t202 * t156 - t157 * t178 + (-t155 * t198 - t161 * t176 + t195 + (-t161 * t206 - t200) * t160) * t269 + (-t155 * t206 - t161 * t194 - t177 + (t161 * t198 - t207) * t160) * t270) * t168) * t163, 0, 0, 0; 0, 0.2e1 * (-t181 * t192 + t193 * t271) * t276 + (0.2e1 * t193 * t249 - t244 * t181 * t230 + t238 * t272 + (-t244 * t186 * t228 - t193 * t171 - t192 * t172 - t238 * t268) * t182) * t173, -t242 * t240 * t251 + (t242 * t178 - ((-qJD(5) * t181 - 0.2e1 * t249) * t230 + (t171 * t230 + (t172 - t253) * t228) * t182) * t240) * t173, 0, t251 + 0.2e1 * (t171 * t182 * t173 + (-t173 * t274 - t182 * t276) * t186) * t186, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:17
	% DurationCPUTime: 2.56s
	% Computational Cost: add. (11602->157), mult. (21125->306), div. (784->12), fcn. (27225->13), ass. (0->123)
	t264 = qJ(3) + pkin(11);
	t262 = sin(t264);
	t263 = cos(t264);
	t267 = sin(qJ(2));
	t269 = cos(qJ(2));
	t325 = cos(pkin(10));
	t326 = cos(pkin(6));
	t290 = t326 * t325;
	t324 = sin(pkin(10));
	t278 = -t267 * t290 - t324 * t269;
	t265 = sin(pkin(6));
	t295 = t265 * t325;
	t238 = t262 * t278 - t263 * t295;
	t288 = t324 * t267 - t269 * t290;
	t254 = t288 * qJD(2);
	t219 = t238 * qJD(3) - t254 * t263;
	t268 = cos(qJ(5));
	t280 = t262 * t295 + t263 * t278;
	t266 = sin(qJ(5));
	t282 = t288 * t266;
	t229 = -t268 * t280 + t282;
	t276 = qJD(2) * t278;
	t201 = t229 * qJD(5) + t219 * t266 + t268 * t276;
	t281 = t288 * t268;
	t227 = -t266 * t280 - t281;
	t222 = t227 ^ 2;
	t311 = t265 * t267;
	t250 = t326 * t262 + t263 * t311;
	t309 = t268 * t269;
	t245 = t250 * t266 + t265 * t309;
	t243 = 0.1e1 / t245 ^ 2;
	t213 = t222 * t243 + 0.1e1;
	t211 = 0.1e1 / t213;
	t249 = -t262 * t311 + t326 * t263;
	t307 = qJD(2) * t265;
	t296 = t269 * t307;
	t236 = t249 * qJD(3) + t263 * t296;
	t310 = t266 * t269;
	t246 = t250 * t268 - t265 * t310;
	t297 = t267 * t307;
	t215 = t246 * qJD(5) + t236 * t266 - t268 * t297;
	t242 = 0.1e1 / t245;
	t315 = t227 * t243;
	t188 = (-t201 * t242 + t215 * t315) * t211;
	t214 = atan2(-t227, t245);
	t208 = sin(t214);
	t209 = cos(t214);
	t287 = -t208 * t245 - t209 * t227;
	t184 = t287 * t188 - t201 * t208 + t209 * t215;
	t200 = -t208 * t227 + t209 * t245;
	t197 = 0.1e1 / t200;
	t198 = 0.1e1 / t200 ^ 2;
	t330 = t184 * t197 * t198;
	t289 = t326 * t324;
	t277 = -t325 * t267 - t269 * t289;
	t259 = -t267 * t289 + t325 * t269;
	t294 = t265 * t324;
	t279 = -t259 * t263 - t262 * t294;
	t230 = -t266 * t279 + t268 * t277;
	t329 = -0.2e1 * t230;
	t292 = 0.2e1 * t230 * t330;
	t283 = -t238 * t242 + t249 * t315;
	t328 = t266 * t283;
	t318 = t215 * t242 * t243;
	t327 = -0.2e1 * (t201 * t315 - t222 * t318) / t213 ^ 2;
	t313 = t277 * t266;
	t231 = -t268 * t279 - t313;
	t224 = 0.1e1 / t231;
	t225 = 0.1e1 / t231 ^ 2;
	t240 = -t259 * t262 + t263 * t294;
	t237 = t240 ^ 2;
	t317 = t225 * t237;
	t210 = 0.1e1 + t317;
	t255 = t277 * qJD(2);
	t220 = t279 * qJD(3) - t255 * t262;
	t221 = t240 * qJD(3) + t255 * t263;
	t256 = t259 * qJD(2);
	t204 = -t230 * qJD(5) + t221 * t268 + t256 * t266;
	t321 = t204 * t224 * t225;
	t299 = t237 * t321;
	t316 = t225 * t240;
	t323 = (t220 * t316 - t299) / t210 ^ 2;
	t322 = t198 * t230;
	t320 = t208 * t230;
	t319 = t209 * t230;
	t314 = t240 * t266;
	t312 = t263 * t268;
	t308 = qJD(2) * t263;
	t306 = qJD(3) * t262;
	t305 = qJD(5) * t263;
	t304 = qJD(5) * t266;
	t303 = qJD(5) * t268;
	t223 = t230 ^ 2;
	t196 = t198 * t223 + 0.1e1;
	t203 = t231 * qJD(5) + t221 * t266 - t256 * t268;
	t302 = 0.2e1 * (t203 * t322 - t223 * t330) / t196 ^ 2;
	t301 = 0.2e1 * t323;
	t298 = t240 * t321;
	t291 = -0.2e1 * t227 * t318;
	t285 = -t229 * t242 + t246 * t315;
	t232 = -t263 * t282 + t268 * t278;
	t247 = (t263 * t310 - t267 * t268) * t265;
	t284 = -t232 * t242 + t247 * t315;
	t235 = -t250 * qJD(3) - t262 * t296;
	t234 = t259 * t266 + t277 * t312;
	t233 = -t259 * t268 + t263 * t313;
	t218 = t280 * qJD(3) + t254 * t262;
	t217 = ((-qJD(2) + t305) * t309 + (-t269 * t306 + (qJD(5) - t308) * t267) * t266) * t265;
	t216 = -t245 * qJD(5) + t236 * t268 + t266 * t297;
	t206 = 0.1e1 / t210;
	t205 = -t288 * t263 * t303 + t254 * t268 + t282 * t306 + (t266 * t308 - t304) * t278;
	t202 = qJD(5) * t281 + t219 * t268 - t266 * t276 + t280 * t304;
	t194 = 0.1e1 / t196;
	t193 = t211 * t328;
	t192 = t284 * t211;
	t190 = t285 * t211;
	t187 = (-t208 * t238 + t209 * t249) * t266 + t287 * t193;
	t186 = t287 * t192 - t208 * t232 + t209 * t247;
	t185 = t287 * t190 - t208 * t229 + t209 * t246;
	t183 = t284 * t327 + (t247 * t291 - t205 * t242 + (t201 * t247 + t215 * t232 + t217 * t227) * t243) * t211;
	t181 = t285 * t327 + (t246 * t291 - t202 * t242 + (t201 * t246 + t215 * t229 + t216 * t227) * t243) * t211;
	t180 = t327 * t328 + (t283 * t303 + (t249 * t291 - t218 * t242 + (t201 * t249 + t215 * t238 + t227 * t235) * t243) * t266) * t211;
	t1 = [0, t183, t180, 0, t181, 0; 0, (t186 * t322 - t197 * t233) * t302 + (t186 * t292 + (-t233 * t184 - t186 * t203 - (-t183 * t227 - t192 * t201 + t217 + (-t192 * t245 - t232) * t188) * t319 - (-t183 * t245 - t192 * t215 - t205 + (t192 * t227 - t247) * t188) * t320) * t198 + ((t277 * t305 - t255) * t268 + (qJD(5) * t259 - t256 * t263 - t277 * t306) * t266) * t197) * t194, (t187 * t322 - t197 * t314) * t302 + ((t220 * t266 + t240 * t303) * t197 + t187 * t292 + (-t187 * t203 - t314 * t184 - (t249 * t303 - t180 * t227 - t193 * t201 + t235 * t266 + (-t193 * t245 - t238 * t266) * t188) * t319 - (-t238 * t303 - t180 * t245 - t193 * t215 - t218 * t266 + (t193 * t227 - t249 * t266) * t188) * t320) * t198) * t194, 0, (t185 * t322 - t197 * t231) * t302 + (t185 * t292 + t204 * t197 + (-t231 * t184 - t185 * t203 - (-t181 * t227 - t190 * t201 + t216 + (-t190 * t245 - t229) * t188) * t319 - (-t181 * t245 - t190 * t215 - t202 + (t190 * t227 - t246) * t188) * t320) * t198) * t194, 0; 0, (t224 * t262 * t277 + t234 * t316) * t301 + (0.2e1 * t234 * t298 + (-qJD(3) * t263 * t277 + t256 * t262) * t224 + (-(t255 * t266 - t256 * t312 + t259 * t303) * t240 - t234 * t220 - (-t262 * t204 - (t263 * t304 + t268 * t306) * t240) * t277) * t225) * t206, (-t224 * t279 + t268 * t317) * t301 + (0.2e1 * t268 * t299 - t221 * t224 + (-0.2e1 * t220 * t240 * t268 - t204 * t279 + t237 * t304) * t225) * t206, 0, t316 * t323 * t329 + (t298 * t329 + (t203 * t240 + t220 * t230) * t225) * t206, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end