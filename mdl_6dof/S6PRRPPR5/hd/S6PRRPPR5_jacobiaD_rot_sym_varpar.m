% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR5
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
%   Wie in S6PRRPPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
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
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:36
	% DurationCPUTime: 0.56s
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
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:36
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (2404->89), mult. (7372->192), div. (524->12), fcn. (9540->11), ass. (0->86)
	t161 = sin(pkin(10));
	t163 = cos(pkin(10));
	t166 = sin(qJ(2));
	t164 = cos(pkin(6));
	t168 = cos(qJ(2));
	t189 = t164 * t168;
	t151 = -t161 * t166 + t163 * t189;
	t141 = t151 * qJD(2);
	t190 = t164 * t166;
	t152 = t161 * t168 + t163 * t190;
	t165 = sin(qJ(3));
	t162 = sin(pkin(6));
	t193 = t162 * t165;
	t181 = t163 * t193;
	t167 = cos(qJ(3));
	t186 = qJD(3) * t167;
	t118 = -qJD(3) * t181 + t141 * t165 + t152 * t186;
	t192 = t162 * t167;
	t134 = t152 * t165 + t163 * t192;
	t131 = t134 ^ 2;
	t155 = -t164 * t167 + t166 * t193;
	t149 = 0.1e1 / t155 ^ 2;
	t129 = t131 * t149 + 0.1e1;
	t126 = 0.1e1 / t129;
	t156 = t164 * t165 + t166 * t192;
	t187 = qJD(2) * t168;
	t180 = t162 * t187;
	t139 = t156 * qJD(3) + t165 * t180;
	t148 = 0.1e1 / t155;
	t197 = t134 * t149;
	t106 = (-t118 * t148 + t139 * t197) * t126;
	t130 = atan2(-t134, t155);
	t122 = sin(t130);
	t123 = cos(t130);
	t178 = -t122 * t155 - t123 * t134;
	t103 = t178 * t106 - t122 * t118 + t123 * t139;
	t117 = -t122 * t134 + t123 * t155;
	t114 = 0.1e1 / t117;
	t115 = 0.1e1 / t117 ^ 2;
	t205 = t103 * t114 * t115;
	t182 = t161 * t190;
	t154 = t163 * t168 - t182;
	t176 = -t154 * t165 + t161 * t192;
	t204 = -0.2e1 * t176 * t205;
	t191 = t162 * t168;
	t175 = -t148 * t151 + t191 * t197;
	t203 = t165 * t175;
	t195 = t139 * t148 * t149;
	t202 = -0.2e1 * (t118 * t197 - t131 * t195) / t129 ^ 2;
	t153 = t161 * t189 + t163 * t166;
	t145 = 0.1e1 / t153;
	t146 = 0.1e1 / t153 ^ 2;
	t201 = t115 * t176;
	t138 = t154 * t167 + t161 * t193;
	t143 = t153 * qJD(2);
	t120 = t138 * qJD(3) - t143 * t165;
	t200 = t120 * t115;
	t199 = t122 * t176;
	t198 = t123 * t176;
	t196 = t138 * t154;
	t194 = t153 * t165;
	t188 = qJD(2) * t166;
	t132 = t176 ^ 2;
	t112 = t132 * t115 + 0.1e1;
	t185 = 0.2e1 * (-t132 * t205 - t176 * t200) / t112 ^ 2;
	t121 = t176 * qJD(3) - t143 * t167;
	t133 = t138 ^ 2;
	t128 = t133 * t146 + 0.1e1;
	t144 = -qJD(2) * t182 + t163 * t187;
	t147 = t145 * t146;
	t184 = 0.2e1 * (t138 * t146 * t121 - t133 * t147 * t144) / t128 ^ 2;
	t179 = -0.2e1 * t134 * t195;
	t136 = t152 * t167 - t181;
	t177 = -t136 * t148 + t156 * t197;
	t142 = t152 * qJD(2);
	t140 = -t155 * qJD(3) + t167 * t180;
	t124 = 0.1e1 / t128;
	t119 = -t134 * qJD(3) + t141 * t167;
	t109 = 0.1e1 / t112;
	t108 = t126 * t203;
	t107 = t177 * t126;
	t105 = (-t122 * t151 + t123 * t191) * t165 + t178 * t108;
	t104 = t178 * t107 - t122 * t136 + t123 * t156;
	t102 = t177 * t202 + (t156 * t179 - t119 * t148 + (t118 * t156 + t134 * t140 + t136 * t139) * t149) * t126;
	t100 = t202 * t203 + (t175 * t186 + (t179 * t191 + t142 * t148 + (t139 * t151 + (t118 * t168 - t134 * t188) * t162) * t149) * t165) * t126;
	t1 = [0, t100, t102, 0, 0, 0; 0, (-t105 * t201 + t114 * t194) * t185 + ((-t144 * t165 - t153 * t186) * t114 + (-t200 + t204) * t105 + (t194 * t103 + (-t100 * t134 - t108 * t118 + (-t165 * t188 + t168 * t186) * t162 + (-t108 * t155 - t151 * t165) * t106) * t198 + (-t151 * t186 - t100 * t155 - t108 * t139 + t142 * t165 + (t108 * t134 - t165 * t191) * t106) * t199) * t115) * t109, (-t104 * t201 - t114 * t138) * t185 + (t104 * t204 + t121 * t114 + (-t138 * t103 - t104 * t120 + (-t102 * t134 - t107 * t118 + t140 + (-t107 * t155 - t136) * t106) * t198 + (-t102 * t155 - t107 * t139 - t119 + (t107 * t134 - t156) * t106) * t199) * t115) * t109, 0, 0, 0; 0, (t145 * t153 * t167 + t146 * t196) * t184 + (qJD(3) * t145 * t194 + (-t121 * t154 + t138 * t143) * t146 + (0.2e1 * t147 * t196 + (t146 * t153 - t145) * t167) * t144) * t124, -t176 * t145 * t184 + (-t144 * t146 * t176 - t120 * t145) * t124, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:36
	% DurationCPUTime: 1.32s
	% Computational Cost: add. (2688->99), mult. (8196->209), div. (535->12), fcn. (10572->13), ass. (0->98)
	t190 = cos(pkin(6));
	t194 = cos(qJ(2));
	t238 = cos(pkin(10));
	t212 = t238 * t194;
	t187 = sin(pkin(10));
	t192 = sin(qJ(2));
	t226 = t187 * t192;
	t177 = t190 * t212 - t226;
	t170 = t177 * qJD(2);
	t193 = cos(qJ(3));
	t191 = sin(qJ(3));
	t213 = t238 * t192;
	t225 = t187 * t194;
	t202 = -t190 * t213 - t225;
	t188 = sin(pkin(6));
	t214 = t188 * t238;
	t240 = t202 * t191 - t193 * t214;
	t143 = qJD(3) * t240 + t170 * t193;
	t163 = -t191 * t214 - t202 * t193;
	t160 = t163 ^ 2;
	t223 = t188 * t193;
	t181 = t190 * t191 + t192 * t223;
	t175 = 0.1e1 / t181 ^ 2;
	t156 = t160 * t175 + 0.1e1;
	t154 = 0.1e1 / t156;
	t224 = t188 * t191;
	t180 = t190 * t193 - t192 * t224;
	t222 = t188 * t194;
	t215 = qJD(2) * t222;
	t168 = t180 * qJD(3) + t193 * t215;
	t174 = 0.1e1 / t181;
	t230 = t163 * t175;
	t126 = (-t143 * t174 + t168 * t230) * t154;
	t157 = atan2(-t163, t181);
	t152 = sin(t157);
	t153 = cos(t157);
	t208 = -t152 * t181 - t153 * t163;
	t122 = t208 * t126 - t152 * t143 + t153 * t168;
	t136 = -t152 * t163 + t153 * t181;
	t133 = 0.1e1 / t136;
	t134 = 0.1e1 / t136 ^ 2;
	t243 = t122 * t133 * t134;
	t179 = -t190 * t226 + t212;
	t166 = t179 * t193 + t187 * t224;
	t242 = 0.2e1 * t166 * t243;
	t204 = -t174 * t177 + t222 * t230;
	t241 = t193 * t204;
	t229 = t168 * t174 * t175;
	t239 = -0.2e1 * (t143 * t230 - t160 * t229) / t156 ^ 2;
	t186 = sin(pkin(11));
	t189 = cos(pkin(11));
	t203 = -t190 * t225 - t213;
	t206 = -t179 * t191 + t187 * t223;
	t151 = -t186 * t206 - t189 * t203;
	t147 = 0.1e1 / t151;
	t148 = 0.1e1 / t151 ^ 2;
	t237 = t134 * t166;
	t172 = t203 * qJD(2);
	t144 = t166 * qJD(3) + t172 * t191;
	t173 = t179 * qJD(2);
	t141 = t144 * t186 + t173 * t189;
	t236 = t141 * t147 * t148;
	t235 = t147 * t189;
	t150 = -t186 * t203 + t189 * t206;
	t234 = t148 * t150;
	t233 = t150 * t186;
	t232 = t152 * t166;
	t231 = t153 * t166;
	t228 = t203 * t191;
	t227 = t203 * t193;
	t221 = qJD(2) * t192;
	t220 = qJD(3) * t191;
	t161 = t166 ^ 2;
	t132 = t134 * t161 + 0.1e1;
	t145 = t206 * qJD(3) + t172 * t193;
	t219 = 0.2e1 * (t145 * t237 - t161 * t243) / t132 ^ 2;
	t146 = t150 ^ 2;
	t139 = t146 * t148 + 0.1e1;
	t140 = -t144 * t189 + t173 * t186;
	t218 = 0.2e1 * (t140 * t234 - t146 * t236) / t139 ^ 2;
	t211 = 0.2e1 * t150 * t236;
	t210 = -0.2e1 * t163 * t229;
	t207 = -t174 * t240 + t180 * t230;
	t205 = -qJD(3) * t227 + t173 * t191;
	t171 = t202 * qJD(2);
	t167 = -t181 * qJD(3) - t191 * t215;
	t159 = t179 * t189 + t186 * t228;
	t158 = t179 * t186 - t189 * t228;
	t142 = t163 * qJD(3) + t170 * t191;
	t137 = 0.1e1 / t139;
	t129 = 0.1e1 / t132;
	t128 = t154 * t241;
	t127 = t207 * t154;
	t124 = (-t152 * t177 + t153 * t222) * t193 + t208 * t128;
	t123 = t208 * t127 - t152 * t240 + t153 * t180;
	t121 = t207 * t239 + (t180 * t210 + t142 * t174 + (t143 * t180 + t163 * t167 + t168 * t240) * t175) * t154;
	t119 = t239 * t241 + (-t204 * t220 + (t210 * t222 - t171 * t174 + (t168 * t177 + (t143 * t194 - t163 * t221) * t188) * t175) * t193) * t154;
	t1 = [0, t119, t121, 0, 0, 0; 0, (t124 * t237 - t133 * t227) * t219 + ((-t173 * t193 - t203 * t220) * t133 + t124 * t242 + (-t124 * t145 - t227 * t122 - (-t119 * t163 - t128 * t143 + (-t193 * t221 - t194 * t220) * t188 + (-t128 * t181 - t177 * t193) * t126) * t231 - (t177 * t220 - t119 * t181 - t128 * t168 - t171 * t193 + (t128 * t163 - t193 * t222) * t126) * t232) * t134) * t129, (t123 * t237 - t133 * t206) * t219 + (t123 * t242 - t144 * t133 + (-t206 * t122 - t123 * t145 - (-t121 * t163 - t127 * t143 + t167 + (-t127 * t181 - t240) * t126) * t231 - (-t121 * t181 - t127 * t168 + t142 + (t127 * t163 - t180) * t126) * t232) * t134) * t129, 0, 0, 0; 0, (-t147 * t158 + t159 * t234) * t218 + ((t172 * t186 + t205 * t189) * t147 + t159 * t211 + (-t158 * t141 - (t172 * t189 - t205 * t186) * t150 - t159 * t140) * t148) * t137, (t148 * t233 + t235) * t166 * t218 + (t166 * t186 * t211 - t145 * t235 + (-t145 * t233 + (-t140 * t186 + t141 * t189) * t166) * t148) * t137, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:37
	% DurationCPUTime: 1.43s
	% Computational Cost: add. (3313->108), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->104)
	t220 = cos(pkin(6));
	t224 = cos(qJ(2));
	t273 = cos(pkin(10));
	t245 = t273 * t224;
	t218 = sin(pkin(10));
	t222 = sin(qJ(2));
	t259 = t218 * t222;
	t206 = t220 * t245 - t259;
	t199 = t206 * qJD(2);
	t223 = cos(qJ(3));
	t221 = sin(qJ(3));
	t246 = t273 * t222;
	t258 = t218 * t224;
	t233 = -t220 * t246 - t258;
	t219 = sin(pkin(6));
	t247 = t219 * t273;
	t275 = t233 * t221 - t223 * t247;
	t177 = t275 * qJD(3) + t199 * t223;
	t192 = -t221 * t247 - t233 * t223;
	t189 = t192 ^ 2;
	t256 = t219 * t223;
	t210 = t220 * t221 + t222 * t256;
	t204 = 0.1e1 / t210 ^ 2;
	t184 = t189 * t204 + 0.1e1;
	t182 = 0.1e1 / t184;
	t257 = t219 * t221;
	t209 = t220 * t223 - t222 * t257;
	t255 = t219 * t224;
	t248 = qJD(2) * t255;
	t197 = t209 * qJD(3) + t223 * t248;
	t203 = 0.1e1 / t210;
	t263 = t192 * t204;
	t154 = (-t177 * t203 + t197 * t263) * t182;
	t185 = atan2(-t192, t210);
	t180 = sin(t185);
	t181 = cos(t185);
	t240 = -t180 * t210 - t181 * t192;
	t150 = t240 * t154 - t180 * t177 + t181 * t197;
	t166 = -t180 * t192 + t181 * t210;
	t163 = 0.1e1 / t166;
	t164 = 0.1e1 / t166 ^ 2;
	t279 = t150 * t163 * t164;
	t217 = pkin(11) + qJ(6);
	t215 = sin(t217);
	t216 = cos(t217);
	t234 = -t220 * t258 - t246;
	t208 = -t220 * t259 + t245;
	t236 = -t208 * t221 + t218 * t256;
	t239 = t215 * t234 - t216 * t236;
	t278 = t239 * qJD(6);
	t195 = t208 * t223 + t218 * t257;
	t277 = 0.2e1 * t195 * t279;
	t235 = -t203 * t206 + t255 * t263;
	t276 = t223 * t235;
	t262 = t197 * t203 * t204;
	t274 = -0.2e1 * (t177 * t263 - t189 * t262) / t184 ^ 2;
	t175 = -t215 * t236 - t216 * t234;
	t171 = 0.1e1 / t175;
	t172 = 0.1e1 / t175 ^ 2;
	t201 = t234 * qJD(2);
	t178 = t195 * qJD(3) + t201 * t221;
	t202 = t208 * qJD(2);
	t161 = t175 * qJD(6) - t178 * t216 + t202 * t215;
	t170 = t239 ^ 2;
	t169 = t170 * t172 + 0.1e1;
	t268 = t172 * t239;
	t162 = t178 * t215 + t202 * t216 + t278;
	t271 = t162 * t171 * t172;
	t272 = (-t161 * t268 - t170 * t271) / t169 ^ 2;
	t270 = t164 * t195;
	t269 = t171 * t216;
	t267 = t239 * t215;
	t179 = t236 * qJD(3) + t201 * t223;
	t266 = t179 * t164;
	t265 = t180 * t195;
	t264 = t181 * t195;
	t261 = t234 * t221;
	t260 = t234 * t223;
	t254 = qJD(2) * t222;
	t253 = qJD(3) * t221;
	t190 = t195 ^ 2;
	t160 = t190 * t164 + 0.1e1;
	t252 = 0.2e1 * (-t190 * t279 + t195 * t266) / t160 ^ 2;
	t251 = 0.2e1 * t272;
	t244 = -0.2e1 * t239 * t271;
	t243 = -0.2e1 * t192 * t262;
	t241 = qJD(6) * t261 + t201;
	t238 = -t172 * t267 + t269;
	t237 = -t203 * t275 + t209 * t263;
	t232 = -qJD(3) * t260 + qJD(6) * t208 + t202 * t221;
	t200 = t233 * qJD(2);
	t196 = -t210 * qJD(3) - t221 * t248;
	t187 = t208 * t216 + t215 * t261;
	t186 = t208 * t215 - t216 * t261;
	t176 = t192 * qJD(3) + t199 * t221;
	t167 = 0.1e1 / t169;
	t157 = 0.1e1 / t160;
	t156 = t182 * t276;
	t155 = t237 * t182;
	t152 = (-t180 * t206 + t181 * t255) * t223 + t240 * t156;
	t151 = t240 * t155 - t180 * t275 + t181 * t209;
	t149 = t237 * t274 + (t209 * t243 + t176 * t203 + (t177 * t209 + t192 * t196 + t197 * t275) * t204) * t182;
	t147 = t274 * t276 + (-t235 * t253 + (t243 * t255 - t200 * t203 + (t197 * t206 + (t177 * t224 - t192 * t254) * t219) * t204) * t223) * t182;
	t1 = [0, t147, t149, 0, 0, 0; 0, (t152 * t270 - t163 * t260) * t252 + ((-t202 * t223 - t234 * t253) * t163 + (-t266 + t277) * t152 + (-t260 * t150 - (-t147 * t192 - t156 * t177 + (-t223 * t254 - t224 * t253) * t219 + (-t156 * t210 - t206 * t223) * t154) * t264 - (t206 * t253 - t147 * t210 - t156 * t197 - t200 * t223 + (t156 * t192 - t223 * t255) * t154) * t265) * t164) * t157, (t151 * t270 - t163 * t236) * t252 + (t151 * t277 - t178 * t163 + (-t236 * t150 - t151 * t179 - (-t149 * t192 - t155 * t177 + t196 + (-t155 * t210 - t275) * t154) * t264 - (-t149 * t210 - t155 * t197 + t176 + (t155 * t192 - t209) * t154) * t265) * t164) * t157, 0, 0, 0; 0, (-t171 * t186 - t187 * t268) * t251 + (t187 * t244 + t241 * t171 * t215 + t232 * t269 + (t216 * t239 * t241 - t187 * t161 - t186 * t162 - t232 * t267) * t172) * t167, t238 * t195 * t251 + (-t238 * t179 + ((qJD(6) * t171 + t244) * t215 + (-t161 * t215 + (t162 + t278) * t216) * t172) * t195) * t167, 0, 0, -0.2e1 * t272 - 0.2e1 * (t161 * t172 * t167 - (-t167 * t271 - t172 * t272) * t239) * t239;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end