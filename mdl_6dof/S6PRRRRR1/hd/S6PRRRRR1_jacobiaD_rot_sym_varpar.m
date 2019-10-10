% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR1
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
%   Wie in S6PRRRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(12));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(12)) * cos(pkin(6));
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
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (756->55), mult. (2271->133), div. (423->14), fcn. (2956->11), ass. (0->65)
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t136 = sin(pkin(12));
	t166 = cos(pkin(6));
	t155 = t136 * t166;
	t165 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (1157->57), mult. (2577->133), div. (441->14), fcn. (3313->11), ass. (0->69)
	t172 = sin(qJ(2));
	t173 = cos(qJ(2));
	t170 = sin(pkin(12));
	t203 = cos(pkin(6));
	t188 = t170 * t203;
	t202 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (1740->57), mult. (2883->133), div. (459->14), fcn. (3670->11), ass. (0->69)
	t176 = sin(qJ(2));
	t177 = cos(qJ(2));
	t174 = sin(pkin(12));
	t207 = cos(pkin(6));
	t192 = t174 * t207;
	t206 = cos(pkin(12));
	t161 = -t176 * t192 + t206 * t177;
	t188 = t207 * t206;
	t157 = t174 * t176 - t177 * t188;
	t175 = sin(pkin(6));
	t196 = t175 * t177;
	t147 = atan2(-t157, -t196);
	t145 = sin(t147);
	t146 = cos(t147);
	t132 = -t145 * t157 - t146 * t196;
	t129 = 0.1e1 / t132;
	t168 = qJ(3) + qJ(4) + qJ(5);
	t165 = sin(t168);
	t166 = cos(t168);
	t197 = t174 * t175;
	t144 = t161 * t166 + t165 * t197;
	t140 = 0.1e1 / t144;
	t171 = 0.1e1 / t177;
	t130 = 0.1e1 / t132 ^ 2;
	t141 = 0.1e1 / t144 ^ 2;
	t172 = 0.1e1 / t177 ^ 2;
	t159 = t174 * t177 + t176 * t188;
	t152 = t159 * qJD(2);
	t198 = t172 * t176;
	t193 = t157 * t198;
	t155 = t157 ^ 2;
	t170 = 0.1e1 / t175 ^ 2;
	t150 = t155 * t170 * t172 + 0.1e1;
	t148 = 0.1e1 / t150;
	t169 = 0.1e1 / t175;
	t200 = t148 * t169;
	t124 = (qJD(2) * t193 + t152 * t171) * t200;
	t186 = t145 * t196 - t146 * t157;
	t194 = t146 * t175 * t176;
	t121 = qJD(2) * t194 + t186 * t124 - t145 * t152;
	t205 = t121 * t129 * t130;
	t143 = t161 * t165 - t166 * t197;
	t139 = t143 ^ 2;
	t135 = t139 * t141 + 0.1e1;
	t184 = -t206 * t176 - t177 * t192;
	t153 = t184 * qJD(2);
	t167 = qJD(3) + qJD(4) + qJD(5);
	t189 = t167 * t197 + t153;
	t199 = t161 * t167;
	t136 = t189 * t165 + t166 * t199;
	t201 = t141 * t143;
	t137 = -t165 * t199 + t189 * t166;
	t202 = t137 * t140 * t141;
	t204 = (t136 * t201 - t139 * t202) / t135 ^ 2;
	t203 = t130 * t184;
	t195 = -0.2e1 * t204;
	t187 = -t140 * t165 + t166 * t201;
	t185 = t159 * t171 + t193;
	t173 = t171 * t172;
	t156 = t184 ^ 2;
	t154 = t161 * qJD(2);
	t151 = qJD(2) * t157;
	t133 = 0.1e1 / t135;
	t128 = t156 * t130 + 0.1e1;
	t125 = t185 * t200;
	t122 = t186 * t125 - t145 * t159 + t194;
	t120 = (-0.2e1 * t185 / t150 ^ 2 * (qJD(2) * t155 * t173 * t176 + t152 * t157 * t172) * t170 + (t152 * t198 - t151 * t171 + (t159 * t198 + (0.2e1 * t173 * t176 ^ 2 + t171) * t157) * qJD(2)) * t148) * t169;
	t118 = t195 + 0.2e1 * (t133 * t136 * t141 + (-t133 * t202 - t141 * t204) * t143) * t143;
	t1 = [0, t120, 0, 0, 0, 0; 0, 0.2e1 * (-t122 * t203 - t129 * t161) / t128 ^ 2 * (-t154 * t203 - t156 * t205) + (t153 * t129 + (-t161 * t121 - t122 * t154) * t130 - (0.2e1 * t122 * t205 + (-(qJD(2) * t196 - t120 * t157 - t125 * t152 + (t125 * t196 - t159) * t124) * t146 - (t124 * t125 * t157 + t151 + (t120 * t177 + (-qJD(2) * t125 - t124) * t176) * t175) * t145) * t130) * t184) / t128, 0, 0, 0, 0; 0, -t187 * t184 * t195 + (t187 * t154 - ((-t140 * t167 - 0.2e1 * t143 * t202) * t166 + (t136 * t166 + (-t143 * t167 + t137) * t165) * t141) * t184) * t133, t118, t118, t118, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:47
	% DurationCPUTime: 2.04s
	% Computational Cost: add. (18064->114), mult. (17539->231), div. (1085->12), fcn. (22566->13), ass. (0->113)
	t266 = sin(pkin(12));
	t268 = cos(pkin(12));
	t271 = sin(qJ(2));
	t269 = cos(pkin(6));
	t273 = cos(qJ(2));
	t300 = t269 * t273;
	t255 = -t266 * t271 + t268 * t300;
	t250 = t255 * qJD(2);
	t301 = t269 * t271;
	t256 = t266 * t273 + t268 * t301;
	t265 = qJ(3) + qJ(4) + qJ(5);
	t262 = sin(t265);
	t264 = qJD(3) + qJD(4) + qJD(5);
	t267 = sin(pkin(6));
	t304 = t267 * t268;
	t289 = t262 * t304;
	t263 = cos(t265);
	t306 = t263 * t264;
	t217 = t250 * t262 + t256 * t306 - t264 * t289;
	t239 = t256 * t262 + t263 * t304;
	t237 = t239 ^ 2;
	t303 = t267 * t271;
	t291 = t262 * t303;
	t247 = -t269 * t263 + t291;
	t245 = 0.1e1 / t247 ^ 2;
	t225 = t237 * t245 + 0.1e1;
	t223 = 0.1e1 / t225;
	t298 = qJD(2) * t273;
	t282 = t264 * t269 + t267 * t298;
	t290 = t263 * t303;
	t235 = t282 * t262 + t264 * t290;
	t244 = 0.1e1 / t247;
	t311 = t239 * t245;
	t201 = (-t217 * t244 + t235 * t311) * t223;
	t230 = atan2(-t239, t247);
	t221 = sin(t230);
	t222 = cos(t230);
	t285 = -t221 * t247 - t222 * t239;
	t197 = t285 * t201 - t221 * t217 + t222 * t235;
	t211 = -t221 * t239 + t222 * t247;
	t208 = 0.1e1 / t211;
	t209 = 0.1e1 / t211 ^ 2;
	t325 = t197 * t208 * t209;
	t292 = t266 * t301;
	t258 = t268 * t273 - t292;
	t305 = t266 * t267;
	t242 = t258 * t262 - t263 * t305;
	t324 = 0.2e1 * t242 * t325;
	t302 = t267 * t273;
	t281 = -t244 * t255 + t302 * t311;
	t323 = t262 * t281;
	t312 = t235 * t244 * t245;
	t322 = -0.2e1 * (t217 * t311 - t237 * t312) / t225 ^ 2;
	t243 = t258 * t263 + t262 * t305;
	t272 = cos(qJ(6));
	t257 = t266 * t300 + t268 * t271;
	t270 = sin(qJ(6));
	t309 = t257 * t270;
	t232 = t243 * t272 + t309;
	t227 = 0.1e1 / t232;
	t228 = 0.1e1 / t232 ^ 2;
	t252 = t257 * qJD(2);
	t287 = t264 * t305 - t252;
	t307 = t262 * t264;
	t220 = -t258 * t307 + t287 * t263;
	t253 = -qJD(2) * t292 + t268 * t298;
	t212 = t232 * qJD(6) + t220 * t270 - t253 * t272;
	t308 = t257 * t272;
	t231 = t243 * t270 - t308;
	t226 = t231 ^ 2;
	t216 = t226 * t228 + 0.1e1;
	t314 = t228 * t231;
	t297 = qJD(6) * t231;
	t213 = t220 * t272 + t253 * t270 - t297;
	t319 = t213 * t227 * t228;
	t321 = (t212 * t314 - t226 * t319) / t216 ^ 2;
	t320 = t209 * t242;
	t219 = t258 * t306 + t287 * t262;
	t318 = t219 * t209;
	t317 = t221 * t242;
	t316 = t222 * t242;
	t315 = t227 * t270;
	t313 = t231 * t272;
	t310 = t257 * t262;
	t299 = qJD(2) * t271;
	t238 = t242 ^ 2;
	t207 = t238 * t209 + 0.1e1;
	t296 = 0.2e1 * (-t238 * t325 + t242 * t318) / t207 ^ 2;
	t295 = -0.2e1 * t321;
	t293 = t231 * t319;
	t288 = -0.2e1 * t239 * t312;
	t286 = qJD(6) * t257 * t263 - t252;
	t284 = t228 * t313 - t315;
	t241 = t256 * t263 - t289;
	t248 = t269 * t262 + t290;
	t283 = -t241 * t244 + t248 * t311;
	t280 = qJD(6) * t258 - t253 * t263 + t257 * t307;
	t251 = t256 * qJD(2);
	t236 = t282 * t263 - t264 * t291;
	t234 = t258 * t270 - t263 * t308;
	t233 = -t258 * t272 - t263 * t309;
	t218 = -t256 * t307 + (-t264 * t304 + t250) * t263;
	t214 = 0.1e1 / t216;
	t205 = 0.1e1 / t207;
	t203 = t223 * t323;
	t202 = t283 * t223;
	t199 = (-t221 * t255 + t222 * t302) * t262 + t285 * t203;
	t198 = t285 * t202 - t221 * t241 + t222 * t248;
	t195 = t283 * t322 + (t248 * t288 - t218 * t244 + (t217 * t248 + t235 * t241 + t236 * t239) * t245) * t223;
	t194 = t322 * t323 + (t281 * t306 + (t288 * t302 + t244 * t251 + (t235 * t255 + (t217 * t273 - t239 * t299) * t267) * t245) * t262) * t223;
	t193 = t284 * t242 * t295 + (t284 * t219 + ((-qJD(6) * t227 - 0.2e1 * t293) * t272 + (t212 * t272 + (t213 - t297) * t270) * t228) * t242) * t214;
	t192 = (t198 * t320 - t208 * t243) * t296 + (t198 * t324 + t220 * t208 + (-t243 * t197 - t198 * t219 - (-t195 * t239 - t202 * t217 + t236 + (-t202 * t247 - t241) * t201) * t316 - (-t195 * t247 - t202 * t235 - t218 + (t202 * t239 - t248) * t201) * t317) * t209) * t205;
	t1 = [0, t194, t195, t195, t195, 0; 0, (t199 * t320 + t208 * t310) * t296 + ((-t253 * t262 - t257 * t306) * t208 + (-t318 + t324) * t199 + (t310 * t197 - (-t194 * t239 - t203 * t217 + (-t262 * t299 + t273 * t306) * t267 + (-t203 * t247 - t255 * t262) * t201) * t316 - (-t255 * t306 - t194 * t247 - t203 * t235 + t251 * t262 + (t203 * t239 - t262 * t302) * t201) * t317) * t209) * t205, t192, t192, t192, 0; 0, 0.2e1 * (-t227 * t233 + t234 * t314) * t321 + (0.2e1 * t234 * t293 - t286 * t227 * t272 + t280 * t315 + (-t231 * t270 * t286 - t234 * t212 - t233 * t213 - t280 * t313) * t228) * t214, t193, t193, t193, t295 + 0.2e1 * (t212 * t228 * t214 + (-t214 * t319 - t228 * t321) * t231) * t231;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end