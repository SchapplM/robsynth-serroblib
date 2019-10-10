% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR1
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
%   Wie in S6PRRRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
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
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.56s
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
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
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
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (1419->57), mult. (2577->133), div. (441->14), fcn. (3313->11), ass. (0->69)
	t177 = sin(qJ(2));
	t178 = cos(qJ(2));
	t175 = sin(pkin(11));
	t208 = cos(pkin(6));
	t193 = t175 * t208;
	t207 = cos(pkin(11));
	t162 = -t177 * t193 + t207 * t178;
	t189 = t208 * t207;
	t158 = t175 * t177 - t178 * t189;
	t176 = sin(pkin(6));
	t197 = t176 * t178;
	t148 = atan2(-t158, -t197);
	t146 = sin(t148);
	t147 = cos(t148);
	t133 = -t146 * t158 - t147 * t197;
	t130 = 0.1e1 / t133;
	t168 = qJ(3) + qJ(4) + pkin(12);
	t166 = sin(t168);
	t167 = cos(t168);
	t198 = t175 * t176;
	t145 = t162 * t167 + t166 * t198;
	t141 = 0.1e1 / t145;
	t172 = 0.1e1 / t178;
	t131 = 0.1e1 / t133 ^ 2;
	t142 = 0.1e1 / t145 ^ 2;
	t173 = 0.1e1 / t178 ^ 2;
	t160 = t175 * t178 + t177 * t189;
	t153 = t160 * qJD(2);
	t199 = t173 * t177;
	t194 = t158 * t199;
	t156 = t158 ^ 2;
	t170 = 0.1e1 / t176 ^ 2;
	t151 = t156 * t170 * t173 + 0.1e1;
	t149 = 0.1e1 / t151;
	t169 = 0.1e1 / t176;
	t201 = t149 * t169;
	t125 = (qJD(2) * t194 + t153 * t172) * t201;
	t187 = t146 * t197 - t147 * t158;
	t195 = t147 * t176 * t177;
	t122 = qJD(2) * t195 + t187 * t125 - t146 * t153;
	t206 = t122 * t130 * t131;
	t185 = -t207 * t177 - t178 * t193;
	t205 = t131 * t185;
	t144 = t162 * t166 - t167 * t198;
	t140 = t144 ^ 2;
	t136 = t140 * t142 + 0.1e1;
	t154 = t185 * qJD(2);
	t171 = qJD(3) + qJD(4);
	t190 = t171 * t198 + t154;
	t200 = t162 * t171;
	t137 = t190 * t166 + t167 * t200;
	t202 = t142 * t144;
	t138 = -t166 * t200 + t190 * t167;
	t203 = t138 * t141 * t142;
	t204 = 0.1e1 / t136 ^ 2 * (t137 * t202 - t140 * t203);
	t196 = -0.2e1 * t204;
	t188 = -t141 * t166 + t167 * t202;
	t186 = t160 * t172 + t194;
	t174 = t172 * t173;
	t157 = t185 ^ 2;
	t155 = t162 * qJD(2);
	t152 = t158 * qJD(2);
	t134 = 0.1e1 / t136;
	t129 = t131 * t157 + 0.1e1;
	t126 = t186 * t201;
	t123 = t187 * t126 - t146 * t160 + t195;
	t121 = (-0.2e1 * t186 / t151 ^ 2 * (qJD(2) * t156 * t174 * t177 + t153 * t158 * t173) * t170 + (t153 * t199 - t152 * t172 + (t160 * t199 + (0.2e1 * t174 * t177 ^ 2 + t172) * t158) * qJD(2)) * t149) * t169;
	t119 = t196 + 0.2e1 * (t134 * t137 * t142 + (-t134 * t203 - t142 * t204) * t144) * t144;
	t1 = [0, t121, 0, 0, 0, 0; 0, 0.2e1 * (-t123 * t205 - t130 * t162) / t129 ^ 2 * (-t155 * t205 - t157 * t206) + (t154 * t130 + (-t162 * t122 - t123 * t155) * t131 - (0.2e1 * t123 * t206 + (-(qJD(2) * t197 - t121 * t158 - t126 * t153 + (t126 * t197 - t160) * t125) * t147 - (t125 * t126 * t158 + t152 + (t121 * t178 + (-qJD(2) * t126 - t125) * t177) * t176) * t146) * t131) * t185) / t129, 0, 0, 0, 0; 0, -t188 * t185 * t196 + (t188 * t155 - ((-t141 * t171 - 0.2e1 * t144 * t203) * t167 + (t137 * t167 + (-t144 * t171 + t138) * t166) * t142) * t185) * t134, t119, t119, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:45
	% DurationCPUTime: 1.76s
	% Computational Cost: add. (13183->114), mult. (13312->231), div. (822->12), fcn. (17117->13), ass. (0->113)
	t263 = sin(pkin(11));
	t265 = cos(pkin(11));
	t268 = sin(qJ(2));
	t266 = cos(pkin(6));
	t270 = cos(qJ(2));
	t297 = t266 * t270;
	t251 = -t263 * t268 + t265 * t297;
	t247 = t251 * qJD(2);
	t298 = t266 * t268;
	t252 = t263 * t270 + t265 * t298;
	t261 = qJ(3) + qJ(4) + pkin(12);
	t259 = sin(t261);
	t262 = qJD(3) + qJD(4);
	t264 = sin(pkin(6));
	t301 = t264 * t265;
	t286 = t259 * t301;
	t260 = cos(t261);
	t303 = t260 * t262;
	t214 = t247 * t259 + t252 * t303 - t262 * t286;
	t236 = t252 * t259 + t260 * t301;
	t234 = t236 ^ 2;
	t300 = t264 * t268;
	t288 = t259 * t300;
	t244 = -t266 * t260 + t288;
	t242 = 0.1e1 / t244 ^ 2;
	t222 = t234 * t242 + 0.1e1;
	t220 = 0.1e1 / t222;
	t295 = qJD(2) * t270;
	t279 = t262 * t266 + t264 * t295;
	t287 = t260 * t300;
	t232 = t279 * t259 + t262 * t287;
	t241 = 0.1e1 / t244;
	t308 = t236 * t242;
	t198 = (-t214 * t241 + t232 * t308) * t220;
	t223 = atan2(-t236, t244);
	t218 = sin(t223);
	t219 = cos(t223);
	t282 = -t218 * t244 - t219 * t236;
	t194 = t282 * t198 - t218 * t214 + t219 * t232;
	t208 = -t218 * t236 + t219 * t244;
	t205 = 0.1e1 / t208;
	t206 = 0.1e1 / t208 ^ 2;
	t322 = t194 * t205 * t206;
	t289 = t263 * t298;
	t254 = t265 * t270 - t289;
	t302 = t263 * t264;
	t239 = t254 * t259 - t260 * t302;
	t321 = 0.2e1 * t239 * t322;
	t299 = t264 * t270;
	t278 = -t241 * t251 + t299 * t308;
	t320 = t259 * t278;
	t309 = t232 * t241 * t242;
	t319 = -0.2e1 * (t214 * t308 - t234 * t309) / t222 ^ 2;
	t240 = t254 * t260 + t259 * t302;
	t269 = cos(qJ(6));
	t253 = t263 * t297 + t265 * t268;
	t267 = sin(qJ(6));
	t306 = t253 * t267;
	t229 = t240 * t269 + t306;
	t225 = 0.1e1 / t229;
	t226 = 0.1e1 / t229 ^ 2;
	t249 = t253 * qJD(2);
	t284 = t262 * t302 - t249;
	t304 = t259 * t262;
	t217 = -t254 * t304 + t284 * t260;
	t250 = -qJD(2) * t289 + t265 * t295;
	t209 = t229 * qJD(6) + t217 * t267 - t250 * t269;
	t305 = t253 * t269;
	t228 = t240 * t267 - t305;
	t224 = t228 ^ 2;
	t213 = t224 * t226 + 0.1e1;
	t311 = t226 * t228;
	t294 = qJD(6) * t228;
	t210 = t217 * t269 + t250 * t267 - t294;
	t316 = t210 * t225 * t226;
	t318 = (t209 * t311 - t224 * t316) / t213 ^ 2;
	t317 = t206 * t239;
	t216 = t254 * t303 + t284 * t259;
	t315 = t216 * t206;
	t314 = t218 * t239;
	t313 = t219 * t239;
	t312 = t225 * t267;
	t310 = t228 * t269;
	t307 = t253 * t259;
	t296 = qJD(2) * t268;
	t235 = t239 ^ 2;
	t204 = t235 * t206 + 0.1e1;
	t293 = 0.2e1 * (-t235 * t322 + t239 * t315) / t204 ^ 2;
	t292 = -0.2e1 * t318;
	t290 = t228 * t316;
	t285 = -0.2e1 * t236 * t309;
	t283 = qJD(6) * t253 * t260 - t249;
	t281 = t226 * t310 - t312;
	t238 = t252 * t260 - t286;
	t245 = t266 * t259 + t287;
	t280 = -t238 * t241 + t245 * t308;
	t277 = qJD(6) * t254 - t250 * t260 + t253 * t304;
	t248 = t252 * qJD(2);
	t233 = t279 * t260 - t262 * t288;
	t231 = t254 * t267 - t260 * t305;
	t230 = -t254 * t269 - t260 * t306;
	t215 = -t252 * t304 + (-t262 * t301 + t247) * t260;
	t211 = 0.1e1 / t213;
	t202 = 0.1e1 / t204;
	t200 = t220 * t320;
	t199 = t280 * t220;
	t196 = (-t218 * t251 + t219 * t299) * t259 + t282 * t200;
	t195 = t282 * t199 - t218 * t238 + t219 * t245;
	t192 = t280 * t319 + (t245 * t285 - t215 * t241 + (t214 * t245 + t232 * t238 + t233 * t236) * t242) * t220;
	t191 = t319 * t320 + (t278 * t303 + (t285 * t299 + t241 * t248 + (t232 * t251 + (t214 * t270 - t236 * t296) * t264) * t242) * t259) * t220;
	t190 = t281 * t239 * t292 + (t281 * t216 + ((-qJD(6) * t225 - 0.2e1 * t290) * t269 + (t209 * t269 + (t210 - t294) * t267) * t226) * t239) * t211;
	t189 = (t195 * t317 - t205 * t240) * t293 + (t195 * t321 + t217 * t205 + (-t240 * t194 - t195 * t216 - (-t192 * t236 - t199 * t214 + t233 + (-t199 * t244 - t238) * t198) * t313 - (-t192 * t244 - t199 * t232 - t215 + (t199 * t236 - t245) * t198) * t314) * t206) * t202;
	t1 = [0, t191, t192, t192, 0, 0; 0, (t196 * t317 + t205 * t307) * t293 + ((-t250 * t259 - t253 * t303) * t205 + (-t315 + t321) * t196 + (t307 * t194 - (-t191 * t236 - t200 * t214 + (-t259 * t296 + t270 * t303) * t264 + (-t200 * t244 - t251 * t259) * t198) * t313 - (-t251 * t303 - t191 * t244 - t200 * t232 + t248 * t259 + (t200 * t236 - t259 * t299) * t198) * t314) * t206) * t202, t189, t189, 0, 0; 0, 0.2e1 * (-t225 * t230 + t231 * t311) * t318 + (0.2e1 * t231 * t290 - t283 * t225 * t269 + t277 * t312 + (-t283 * t228 * t267 - t231 * t209 - t230 * t210 - t277 * t310) * t226) * t211, t190, t190, 0, t292 + 0.2e1 * (t209 * t226 * t211 + (-t211 * t316 - t226 * t318) * t228) * t228;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end