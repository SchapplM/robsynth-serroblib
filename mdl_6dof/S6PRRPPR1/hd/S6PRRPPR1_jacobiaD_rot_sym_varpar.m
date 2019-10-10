% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR1
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
%   Wie in S6PRRPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
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
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:16
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
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:16
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
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:17
	% DurationCPUTime: 1.45s
	% Computational Cost: add. (5351->102), mult. (8196->213), div. (535->12), fcn. (10572->13), ass. (0->102)
	t202 = sin(pkin(10));
	t205 = cos(pkin(10));
	t207 = sin(qJ(2));
	t206 = cos(pkin(6));
	t208 = cos(qJ(2));
	t231 = t206 * t208;
	t190 = -t202 * t207 + t205 * t231;
	t186 = t190 * qJD(2);
	t232 = t206 * t207;
	t191 = t202 * t208 + t205 * t232;
	t200 = qJ(3) + pkin(11);
	t198 = sin(t200);
	t203 = sin(pkin(6));
	t235 = t203 * t205;
	t222 = t198 * t235;
	t199 = cos(t200);
	t228 = qJD(3) * t199;
	t154 = -qJD(3) * t222 + t186 * t198 + t191 * t228;
	t176 = t191 * t198 + t199 * t235;
	t174 = t176 ^ 2;
	t234 = t203 * t207;
	t184 = t198 * t234 - t206 * t199;
	t182 = 0.1e1 / t184 ^ 2;
	t168 = t174 * t182 + 0.1e1;
	t166 = 0.1e1 / t168;
	t185 = t206 * t198 + t199 * t234;
	t229 = qJD(2) * t208;
	t221 = t203 * t229;
	t172 = t185 * qJD(3) + t198 * t221;
	t181 = 0.1e1 / t184;
	t240 = t176 * t182;
	t138 = (-t154 * t181 + t172 * t240) * t166;
	t169 = atan2(-t176, t184);
	t164 = sin(t169);
	t165 = cos(t169);
	t219 = -t164 * t184 - t165 * t176;
	t134 = t219 * t138 - t164 * t154 + t165 * t172;
	t148 = -t164 * t176 + t165 * t184;
	t145 = 0.1e1 / t148;
	t146 = 0.1e1 / t148 ^ 2;
	t254 = t134 * t145 * t146;
	t223 = t202 * t232;
	t193 = t205 * t208 - t223;
	t236 = t202 * t203;
	t180 = t193 * t199 + t198 * t236;
	t201 = sin(pkin(12));
	t192 = t202 * t231 + t205 * t207;
	t204 = cos(pkin(12));
	t237 = t192 * t204;
	t162 = t180 * t201 - t237;
	t188 = t192 * qJD(2);
	t217 = -t193 * t198 + t199 * t236;
	t157 = t217 * qJD(3) - t188 * t199;
	t189 = -qJD(2) * t223 + t205 * t229;
	t153 = t157 * t204 + t189 * t201;
	t238 = t192 * t201;
	t163 = t180 * t204 + t238;
	t159 = 0.1e1 / t163;
	t160 = 0.1e1 / t163 ^ 2;
	t248 = t153 * t159 * t160;
	t253 = 0.2e1 * t162 * t248;
	t252 = -0.2e1 * t217 * t254;
	t233 = t203 * t208;
	t215 = -t181 * t190 + t233 * t240;
	t251 = t198 * t215;
	t241 = t172 * t181 * t182;
	t250 = -0.2e1 * (t154 * t240 - t174 * t241) / t168 ^ 2;
	t249 = t146 * t217;
	t156 = t180 * qJD(3) - t188 * t198;
	t247 = t156 * t146;
	t246 = t159 * t201;
	t245 = t160 * t162;
	t244 = t162 * t204;
	t243 = t164 * t217;
	t242 = t165 * t217;
	t239 = t192 * t198;
	t230 = qJD(2) * t207;
	t175 = t217 ^ 2;
	t144 = t175 * t146 + 0.1e1;
	t227 = 0.2e1 * (-t175 * t254 - t217 * t247) / t144 ^ 2;
	t158 = t162 ^ 2;
	t151 = t158 * t160 + 0.1e1;
	t152 = t157 * t201 - t189 * t204;
	t226 = 0.2e1 * (t152 * t245 - t158 * t248) / t151 ^ 2;
	t220 = -0.2e1 * t176 * t241;
	t178 = t191 * t199 - t222;
	t218 = -t178 * t181 + t185 * t240;
	t216 = qJD(3) * t239 - t189 * t199;
	t187 = t191 * qJD(2);
	t173 = -t184 * qJD(3) + t199 * t221;
	t171 = t193 * t201 - t199 * t237;
	t170 = -t193 * t204 - t199 * t238;
	t155 = -t176 * qJD(3) + t186 * t199;
	t149 = 0.1e1 / t151;
	t141 = 0.1e1 / t144;
	t140 = t166 * t251;
	t139 = t218 * t166;
	t136 = (-t164 * t190 + t165 * t233) * t198 + t219 * t140;
	t135 = t219 * t139 - t164 * t178 + t165 * t185;
	t133 = t218 * t250 + (t185 * t220 - t155 * t181 + (t154 * t185 + t172 * t178 + t173 * t176) * t182) * t166;
	t131 = t250 * t251 + (t215 * t228 + (t220 * t233 + t181 * t187 + (t172 * t190 + (t154 * t208 - t176 * t230) * t203) * t182) * t198) * t166;
	t1 = [0, t131, t133, 0, 0, 0; 0, (-t136 * t249 + t145 * t239) * t227 + ((-t189 * t198 - t192 * t228) * t145 + (-t247 + t252) * t136 + (t239 * t134 + (-t131 * t176 - t140 * t154 + (-t198 * t230 + t208 * t228) * t203 + (-t140 * t184 - t190 * t198) * t138) * t242 + (-t190 * t228 - t131 * t184 - t140 * t172 + t187 * t198 + (t140 * t176 - t198 * t233) * t138) * t243) * t146) * t141, (-t135 * t249 - t145 * t180) * t227 + (t135 * t252 + t157 * t145 + (-t180 * t134 - t135 * t156 + (-t133 * t176 - t139 * t154 + t173 + (-t139 * t184 - t178) * t138) * t242 + (-t133 * t184 - t139 * t172 - t155 + (t139 * t176 - t185) * t138) * t243) * t146) * t141, 0, 0, 0; 0, (-t159 * t170 + t171 * t245) * t226 + ((t188 * t204 + t216 * t201) * t159 + t171 * t253 + (-t170 * t153 - (-t188 * t201 + t216 * t204) * t162 - t171 * t152) * t160) * t149, -(-t160 * t244 + t246) * t217 * t226 + (t217 * t204 * t253 - t156 * t246 + (t156 * t244 - (t152 * t204 + t153 * t201) * t217) * t160) * t149, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:16
	% EndTime: 2019-10-09 22:07:17
	% DurationCPUTime: 1.39s
	% Computational Cost: add. (6115->111), mult. (9085->227), div. (559->12), fcn. (11668->13), ass. (0->106)
	t233 = sin(pkin(10));
	t235 = cos(pkin(10));
	t237 = sin(qJ(2));
	t236 = cos(pkin(6));
	t238 = cos(qJ(2));
	t264 = t236 * t238;
	t219 = -t233 * t237 + t235 * t264;
	t215 = t219 * qJD(2);
	t265 = t236 * t237;
	t220 = t233 * t238 + t235 * t265;
	t232 = qJ(3) + pkin(11);
	t228 = sin(t232);
	t234 = sin(pkin(6));
	t268 = t234 * t235;
	t254 = t228 * t268;
	t230 = cos(t232);
	t261 = qJD(3) * t230;
	t188 = -qJD(3) * t254 + t215 * t228 + t220 * t261;
	t204 = t220 * t228 + t230 * t268;
	t202 = t204 ^ 2;
	t267 = t234 * t237;
	t213 = t228 * t267 - t230 * t236;
	t211 = 0.1e1 / t213 ^ 2;
	t196 = t202 * t211 + 0.1e1;
	t194 = 0.1e1 / t196;
	t214 = t228 * t236 + t230 * t267;
	t262 = qJD(2) * t238;
	t253 = t234 * t262;
	t200 = t214 * qJD(3) + t228 * t253;
	t210 = 0.1e1 / t213;
	t272 = t204 * t211;
	t166 = (-t188 * t210 + t200 * t272) * t194;
	t197 = atan2(-t204, t213);
	t192 = sin(t197);
	t193 = cos(t197);
	t250 = -t192 * t213 - t193 * t204;
	t162 = t250 * t166 - t188 * t192 + t193 * t200;
	t176 = -t192 * t204 + t193 * t213;
	t173 = 0.1e1 / t176;
	t174 = 0.1e1 / t176 ^ 2;
	t286 = t162 * t173 * t174;
	t255 = t233 * t265;
	t222 = t235 * t238 - t255;
	t269 = t233 * t234;
	t247 = -t222 * t228 + t230 * t269;
	t285 = -0.2e1 * t247 * t286;
	t266 = t234 * t238;
	t246 = -t210 * t219 + t266 * t272;
	t284 = t228 * t246;
	t273 = t200 * t210 * t211;
	t283 = -0.2e1 * (t188 * t272 - t202 * t273) / t196 ^ 2;
	t208 = t222 * t230 + t228 * t269;
	t221 = t233 * t264 + t235 * t237;
	t231 = pkin(12) + qJ(6);
	t227 = sin(t231);
	t229 = cos(t231);
	t187 = t208 * t229 + t221 * t227;
	t183 = 0.1e1 / t187;
	t184 = 0.1e1 / t187 ^ 2;
	t217 = t221 * qJD(2);
	t191 = t247 * qJD(3) - t217 * t230;
	t218 = -qJD(2) * t255 + t235 * t262;
	t177 = t187 * qJD(6) + t191 * t227 - t218 * t229;
	t186 = t208 * t227 - t221 * t229;
	t182 = t186 ^ 2;
	t181 = t182 * t184 + 0.1e1;
	t278 = t184 * t186;
	t260 = qJD(6) * t186;
	t178 = t191 * t229 + t218 * t227 - t260;
	t280 = t178 * t183 * t184;
	t282 = (t177 * t278 - t182 * t280) / t181 ^ 2;
	t281 = t174 * t247;
	t279 = t183 * t227;
	t277 = t186 * t229;
	t190 = t208 * qJD(3) - t217 * t228;
	t276 = t190 * t174;
	t275 = t192 * t247;
	t274 = t193 * t247;
	t271 = t221 * t228;
	t270 = t221 * t230;
	t263 = qJD(2) * t237;
	t203 = t247 ^ 2;
	t172 = t174 * t203 + 0.1e1;
	t259 = 0.2e1 * (-t203 * t286 - t247 * t276) / t172 ^ 2;
	t258 = -0.2e1 * t282;
	t256 = t186 * t280;
	t252 = -0.2e1 * t204 * t273;
	t251 = qJD(6) * t270 - t217;
	t249 = t184 * t277 - t279;
	t206 = t220 * t230 - t254;
	t248 = -t206 * t210 + t214 * t272;
	t245 = qJD(3) * t271 + qJD(6) * t222 - t218 * t230;
	t216 = t220 * qJD(2);
	t201 = -t213 * qJD(3) + t230 * t253;
	t199 = t222 * t227 - t229 * t270;
	t198 = -t222 * t229 - t227 * t270;
	t189 = -t204 * qJD(3) + t215 * t230;
	t179 = 0.1e1 / t181;
	t169 = 0.1e1 / t172;
	t168 = t194 * t284;
	t167 = t248 * t194;
	t164 = (-t192 * t219 + t193 * t266) * t228 + t250 * t168;
	t163 = t250 * t167 - t192 * t206 + t193 * t214;
	t161 = t248 * t283 + (t214 * t252 - t189 * t210 + (t188 * t214 + t200 * t206 + t201 * t204) * t211) * t194;
	t159 = t283 * t284 + (t246 * t261 + (t252 * t266 + t210 * t216 + (t200 * t219 + (t188 * t238 - t204 * t263) * t234) * t211) * t228) * t194;
	t1 = [0, t159, t161, 0, 0, 0; 0, (-t164 * t281 + t173 * t271) * t259 + ((-t218 * t228 - t221 * t261) * t173 + (-t276 + t285) * t164 + (t271 * t162 + (-t159 * t204 - t168 * t188 + (-t228 * t263 + t238 * t261) * t234 + (-t168 * t213 - t219 * t228) * t166) * t274 + (-t219 * t261 - t159 * t213 - t168 * t200 + t216 * t228 + (t168 * t204 - t228 * t266) * t166) * t275) * t174) * t169, (-t163 * t281 - t173 * t208) * t259 + (t163 * t285 + t191 * t173 + (-t208 * t162 - t163 * t190 + (-t161 * t204 - t167 * t188 + t201 + (-t167 * t213 - t206) * t166) * t274 + (-t161 * t213 - t167 * t200 - t189 + (t167 * t204 - t214) * t166) * t275) * t174) * t169, 0, 0, 0; 0, 0.2e1 * (-t183 * t198 + t199 * t278) * t282 + (0.2e1 * t199 * t256 - t251 * t183 * t229 + t245 * t279 + (-t251 * t186 * t227 - t199 * t177 - t198 * t178 - t245 * t277) * t184) * t179, -t249 * t247 * t258 + (t249 * t190 - ((-qJD(6) * t183 - 0.2e1 * t256) * t229 + (t177 * t229 + (t178 - t260) * t227) * t184) * t247) * t179, 0, 0, t258 + 0.2e1 * (t177 * t179 * t184 + (-t179 * t280 - t184 * t282) * t186) * t186;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end