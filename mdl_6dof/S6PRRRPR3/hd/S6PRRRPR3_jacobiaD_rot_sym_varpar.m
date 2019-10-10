% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR3
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
%   Wie in S6PRRRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:50
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
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
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:25
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
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:25
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
	% StartTime: 2019-10-09 22:50:25
	% EndTime: 2019-10-09 22:50:26
	% DurationCPUTime: 1.40s
	% Computational Cost: add. (7756->94), mult. (10920->196), div. (767->12), fcn. (14140->11), ass. (0->96)
	t204 = sin(pkin(11));
	t206 = cos(pkin(11));
	t208 = sin(qJ(2));
	t207 = cos(pkin(6));
	t209 = cos(qJ(2));
	t231 = t207 * t209;
	t192 = -t204 * t208 + t206 * t231;
	t185 = t192 * qJD(2);
	t232 = t207 * t208;
	t193 = t204 * t209 + t206 * t232;
	t203 = qJ(3) + qJ(4);
	t200 = sin(t203);
	t202 = qJD(3) + qJD(4);
	t205 = sin(pkin(6));
	t235 = t205 * t206;
	t222 = t200 * t235;
	t201 = cos(t203);
	t237 = t201 * t202;
	t157 = t185 * t200 + t193 * t237 - t202 * t222;
	t175 = t193 * t200 + t201 * t235;
	t172 = t175 ^ 2;
	t234 = t205 * t208;
	t224 = t200 * t234;
	t183 = -t207 * t201 + t224;
	t181 = 0.1e1 / t183 ^ 2;
	t165 = t172 * t181 + 0.1e1;
	t163 = 0.1e1 / t165;
	t229 = qJD(2) * t209;
	t217 = t202 * t207 + t205 * t229;
	t223 = t201 * t234;
	t170 = t200 * t217 + t202 * t223;
	t180 = 0.1e1 / t183;
	t242 = t175 * t181;
	t145 = (-t157 * t180 + t170 * t242) * t163;
	t166 = atan2(-t175, t183);
	t161 = sin(t166);
	t162 = cos(t166);
	t219 = -t161 * t183 - t162 * t175;
	t141 = t145 * t219 - t161 * t157 + t162 * t170;
	t156 = -t161 * t175 + t162 * t183;
	t153 = 0.1e1 / t156;
	t154 = 0.1e1 / t156 ^ 2;
	t251 = t141 * t153 * t154;
	t225 = t204 * t232;
	t195 = t206 * t209 - t225;
	t236 = t204 * t205;
	t178 = t195 * t200 - t201 * t236;
	t250 = 0.2e1 * t178 * t251;
	t233 = t205 * t209;
	t216 = -t180 * t192 + t233 * t242;
	t249 = t200 * t216;
	t243 = t170 * t180 * t181;
	t248 = -0.2e1 * (t157 * t242 - t172 * t243) / t165 ^ 2;
	t194 = t204 * t231 + t206 * t208;
	t189 = 0.1e1 / t194;
	t190 = 0.1e1 / t194 ^ 2;
	t247 = t154 * t178;
	t187 = t194 * qJD(2);
	t220 = t202 * t236 - t187;
	t159 = t195 * t237 + t200 * t220;
	t246 = t159 * t154;
	t245 = t161 * t178;
	t244 = t162 * t178;
	t179 = t195 * t201 + t200 * t236;
	t241 = t179 * t195;
	t240 = t189 * t194;
	t239 = t194 * t200;
	t238 = t200 * t202;
	t230 = qJD(2) * t208;
	t173 = t178 ^ 2;
	t151 = t173 * t154 + 0.1e1;
	t228 = 0.2e1 * (-t173 * t251 + t178 * t246) / t151 ^ 2;
	t160 = -t195 * t238 + t201 * t220;
	t174 = t179 ^ 2;
	t169 = t174 * t190 + 0.1e1;
	t188 = -qJD(2) * t225 + t206 * t229;
	t191 = t189 * t190;
	t227 = 0.2e1 * (t179 * t190 * t160 - t174 * t191 * t188) / t169 ^ 2;
	t221 = -0.2e1 * t175 * t243;
	t177 = t193 * t201 - t222;
	t184 = t207 * t200 + t223;
	t218 = -t177 * t180 + t184 * t242;
	t186 = t193 * qJD(2);
	t171 = t201 * t217 - t202 * t224;
	t167 = 0.1e1 / t169;
	t158 = -t193 * t238 + (-t202 * t235 + t185) * t201;
	t149 = 0.1e1 / t151;
	t147 = t163 * t249;
	t146 = t218 * t163;
	t144 = t178 * t189 * t227 + (t178 * t188 * t190 - t159 * t189) * t167;
	t143 = (-t161 * t192 + t162 * t233) * t200 + t219 * t147;
	t142 = t146 * t219 - t161 * t177 + t162 * t184;
	t140 = t218 * t248 + (t184 * t221 - t158 * t180 + (t157 * t184 + t170 * t177 + t171 * t175) * t181) * t163;
	t138 = t248 * t249 + (t216 * t237 + (t221 * t233 + t180 * t186 + (t170 * t192 + (t157 * t209 - t175 * t230) * t205) * t181) * t200) * t163;
	t137 = (t142 * t247 - t153 * t179) * t228 + (t142 * t250 + t160 * t153 + (-t179 * t141 - t142 * t159 - (-t140 * t175 - t146 * t157 + t171 + (-t146 * t183 - t177) * t145) * t244 - (-t140 * t183 - t146 * t170 - t158 + (t146 * t175 - t184) * t145) * t245) * t154) * t149;
	t1 = [0, t138, t140, t140, 0, 0; 0, (t143 * t247 + t153 * t239) * t228 + ((-t188 * t200 - t194 * t237) * t153 + (-t246 + t250) * t143 + (t239 * t141 - (-t138 * t175 - t147 * t157 + (-t200 * t230 + t209 * t237) * t205 + (-t147 * t183 - t192 * t200) * t145) * t244 - (-t192 * t237 - t138 * t183 - t147 * t170 + t186 * t200 + (t147 * t175 - t200 * t233) * t145) * t245) * t154) * t149, t137, t137, 0, 0; 0, (t190 * t241 + t201 * t240) * t227 + (t238 * t240 + (-t160 * t195 + t179 * t187) * t190 + (0.2e1 * t191 * t241 + (t190 * t194 - t189) * t201) * t188) * t167, t144, t144, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:25
	% EndTime: 2019-10-09 22:50:27
	% DurationCPUTime: 1.75s
	% Computational Cost: add. (9025->112), mult. (13312->224), div. (822->12), fcn. (17117->13), ass. (0->114)
	t264 = qJ(3) + qJ(4);
	t262 = cos(t264);
	t263 = qJD(3) + qJD(4);
	t261 = sin(t264);
	t267 = cos(pkin(6));
	t269 = sin(qJ(2));
	t325 = cos(pkin(11));
	t295 = t325 * t269;
	t265 = sin(pkin(11));
	t271 = cos(qJ(2));
	t306 = t265 * t271;
	t281 = -t267 * t295 - t306;
	t279 = t281 * t261;
	t294 = t325 * t271;
	t307 = t265 * t269;
	t254 = t267 * t294 - t307;
	t266 = sin(pkin(6));
	t296 = t266 * t325;
	t327 = t254 * qJD(2) - t263 * t296;
	t217 = t327 * t262 + t263 * t279;
	t278 = t281 * t262;
	t240 = -t261 * t296 - t278;
	t237 = t240 ^ 2;
	t305 = t266 * t269;
	t297 = t262 * t305;
	t249 = t267 * t261 + t297;
	t246 = 0.1e1 / t249 ^ 2;
	t230 = t237 * t246 + 0.1e1;
	t228 = 0.1e1 / t230;
	t304 = t266 * t271;
	t284 = qJD(2) * t304 + t263 * t267;
	t298 = t261 * t305;
	t235 = t262 * t284 - t263 * t298;
	t245 = 0.1e1 / t249;
	t314 = t240 * t246;
	t200 = (-t217 * t245 + t235 * t314) * t228;
	t231 = atan2(-t240, t249);
	t226 = sin(t231);
	t227 = cos(t231);
	t288 = -t226 * t249 - t227 * t240;
	t196 = t200 * t288 - t226 * t217 + t227 * t235;
	t210 = -t226 * t240 + t227 * t249;
	t207 = 0.1e1 / t210;
	t208 = 0.1e1 / t210 ^ 2;
	t331 = t196 * t207 * t208;
	t256 = -t267 * t307 + t294;
	t308 = t265 * t266;
	t242 = t256 * t261 - t262 * t308;
	t270 = cos(qJ(6));
	t268 = sin(qJ(6));
	t282 = -t267 * t306 - t295;
	t312 = t282 * t268;
	t287 = t242 * t270 + t312;
	t330 = qJD(6) * t287;
	t283 = -t245 * t254 + t304 * t314;
	t329 = t262 * t283;
	t243 = t256 * t262 + t261 * t308;
	t328 = 0.2e1 * t243 * t331;
	t315 = t235 * t245 * t246;
	t326 = -0.2e1 * (t217 * t314 - t237 * t315) / t230 ^ 2;
	t311 = t282 * t270;
	t225 = t242 * t268 - t311;
	t221 = 0.1e1 / t225;
	t222 = 0.1e1 / t225 ^ 2;
	t252 = t282 * qJD(2);
	t290 = t263 * t308 + t252;
	t309 = t262 * t263;
	t218 = t256 * t309 + t261 * t290;
	t253 = t256 * qJD(2);
	t211 = qJD(6) * t225 - t218 * t270 + t253 * t268;
	t220 = t287 ^ 2;
	t215 = t220 * t222 + 0.1e1;
	t319 = t222 * t287;
	t212 = t218 * t268 + t253 * t270 + t330;
	t322 = t212 * t221 * t222;
	t324 = (-t211 * t319 - t220 * t322) / t215 ^ 2;
	t323 = t208 * t243;
	t310 = t261 * t263;
	t219 = -t256 * t310 + t262 * t290;
	t321 = t219 * t208;
	t320 = t221 * t270;
	t318 = t287 * t268;
	t317 = t226 * t243;
	t316 = t227 * t243;
	t313 = t282 * t262;
	t303 = qJD(2) * t269;
	t238 = t243 ^ 2;
	t206 = t238 * t208 + 0.1e1;
	t302 = 0.2e1 * (-t238 * t331 + t243 * t321) / t206 ^ 2;
	t301 = 0.2e1 * t324;
	t293 = -0.2e1 * t287 * t322;
	t292 = -0.2e1 * t240 * t315;
	t289 = qJD(6) * t261 * t282 + t252;
	t286 = -t222 * t318 + t320;
	t239 = t262 * t296 - t279;
	t248 = t267 * t262 - t298;
	t285 = t239 * t245 + t248 * t314;
	t280 = qJD(6) * t256 + t253 * t261 - t282 * t309;
	t251 = t281 * qJD(2);
	t234 = -t261 * t284 - t263 * t297;
	t233 = t256 * t270 + t261 * t312;
	t232 = t256 * t268 - t261 * t311;
	t216 = t327 * t261 - t263 * t278;
	t213 = 0.1e1 / t215;
	t204 = 0.1e1 / t206;
	t202 = t228 * t329;
	t201 = t285 * t228;
	t198 = (-t226 * t254 + t227 * t304) * t262 + t288 * t202;
	t197 = t201 * t288 + t226 * t239 + t227 * t248;
	t195 = t285 * t326 + (t248 * t292 + t216 * t245 + (t217 * t248 + t234 * t240 - t235 * t239) * t246) * t228;
	t193 = t326 * t329 + (-t283 * t310 + (t292 * t304 - t245 * t251 + (t235 * t254 + (t217 * t271 - t240 * t303) * t266) * t246) * t262) * t228;
	t192 = t286 * t243 * t301 + (-t286 * t219 + ((qJD(6) * t221 + t293) * t268 + (-t211 * t268 + (t212 + t330) * t270) * t222) * t243) * t213;
	t191 = (t197 * t323 + t207 * t242) * t302 + (t197 * t328 - t218 * t207 + (t242 * t196 - t197 * t219 - (-t195 * t240 - t201 * t217 + t234 + (-t201 * t249 + t239) * t200) * t316 - (-t195 * t249 - t201 * t235 + t216 + (t201 * t240 - t248) * t200) * t317) * t208) * t204;
	t1 = [0, t193, t195, t195, 0, 0; 0, (t198 * t323 - t207 * t313) * t302 + ((-t253 * t262 - t282 * t310) * t207 + (-t321 + t328) * t198 + (-t313 * t196 - (-t193 * t240 - t202 * t217 + (-t262 * t303 - t271 * t310) * t266 + (-t202 * t249 - t254 * t262) * t200) * t316 - (t254 * t310 - t193 * t249 - t202 * t235 - t251 * t262 + (t202 * t240 - t262 * t304) * t200) * t317) * t208) * t204, t191, t191, 0, 0; 0, (-t221 * t232 - t233 * t319) * t301 + (t233 * t293 + t289 * t221 * t268 + t280 * t320 + (t270 * t287 * t289 - t233 * t211 - t232 * t212 - t280 * t318) * t222) * t213, t192, t192, 0, -0.2e1 * t324 - 0.2e1 * (t211 * t222 * t213 - (-t213 * t322 - t222 * t324) * t287) * t287;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end