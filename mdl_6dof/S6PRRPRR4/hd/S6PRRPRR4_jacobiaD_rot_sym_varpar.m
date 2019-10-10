% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR4
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
%   Wie in S6PRRPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
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
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:30
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
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:30
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (2420->88), mult. (7430->188), div. (524->12), fcn. (9601->11), ass. (0->89)
	t164 = sin(pkin(11));
	t166 = cos(pkin(11));
	t169 = sin(qJ(2));
	t167 = cos(pkin(6));
	t171 = cos(qJ(2));
	t193 = t167 * t171;
	t154 = -t164 * t169 + t166 * t193;
	t146 = t154 * qJD(2);
	t194 = t167 * t169;
	t155 = t164 * t171 + t166 * t194;
	t168 = sin(qJ(3));
	t165 = sin(pkin(6));
	t197 = t165 * t168;
	t186 = t166 * t197;
	t170 = cos(qJ(3));
	t191 = qJD(3) * t170;
	t121 = -qJD(3) * t186 + t146 * t168 + t155 * t191;
	t196 = t165 * t170;
	t139 = t155 * t168 + t166 * t196;
	t134 = t139 ^ 2;
	t158 = -t167 * t170 + t169 * t197;
	t152 = 0.1e1 / t158 ^ 2;
	t132 = t134 * t152 + 0.1e1;
	t129 = 0.1e1 / t132;
	t159 = t167 * t168 + t169 * t196;
	t195 = t165 * t171;
	t185 = qJD(2) * t195;
	t144 = t159 * qJD(3) + t168 * t185;
	t151 = 0.1e1 / t158;
	t202 = t139 * t152;
	t109 = (-t121 * t151 + t144 * t202) * t129;
	t133 = atan2(-t139, t158);
	t125 = sin(t133);
	t126 = cos(t133);
	t183 = -t125 * t158 - t126 * t139;
	t106 = t183 * t109 - t125 * t121 + t126 * t144;
	t120 = -t125 * t139 + t126 * t158;
	t117 = 0.1e1 / t120;
	t118 = 0.1e1 / t120 ^ 2;
	t214 = t106 * t117 * t118;
	t179 = t164 * t194 - t166 * t171;
	t181 = t164 * t196 + t168 * t179;
	t213 = -0.2e1 * t181;
	t212 = 0.2e1 * t170;
	t211 = t213 * t214;
	t178 = -t151 * t154 + t195 * t202;
	t210 = t168 * t178;
	t201 = t144 * t151 * t152;
	t209 = -0.2e1 * (t121 * t202 - t134 * t201) / t132 ^ 2;
	t143 = t164 * t197 - t170 * t179;
	t136 = 0.1e1 / t143;
	t137 = 0.1e1 / t143 ^ 2;
	t180 = t164 * t193 + t166 * t169;
	t150 = t180 ^ 2;
	t199 = t150 * t137;
	t131 = 0.1e1 + t199;
	t148 = t180 * qJD(2);
	t124 = t181 * qJD(3) - t148 * t170;
	t205 = t124 * t136 * t137;
	t187 = t150 * t205;
	t149 = t179 * qJD(2);
	t200 = t149 * t180;
	t208 = (-t137 * t200 - t187) / t131 ^ 2;
	t207 = t118 * t181;
	t123 = t143 * qJD(3) - t148 * t168;
	t206 = t123 * t118;
	t204 = t125 * t181;
	t203 = t126 * t181;
	t198 = t180 * t168;
	t192 = qJD(2) * t169;
	t135 = t181 ^ 2;
	t115 = t135 * t118 + 0.1e1;
	t190 = 0.2e1 * (-t135 * t214 - t181 * t206) / t115 ^ 2;
	t188 = t180 * t213;
	t184 = -0.2e1 * t139 * t201;
	t141 = t155 * t170 - t186;
	t182 = -t141 * t151 + t159 * t202;
	t147 = t155 * qJD(2);
	t145 = -t158 * qJD(3) + t170 * t185;
	t127 = 0.1e1 / t131;
	t122 = -t139 * qJD(3) + t146 * t170;
	t112 = 0.1e1 / t115;
	t111 = t129 * t210;
	t110 = t182 * t129;
	t108 = (-t125 * t154 + t126 * t195) * t168 + t183 * t111;
	t107 = t183 * t110 - t125 * t141 + t126 * t159;
	t105 = t182 * t209 + (t159 * t184 - t122 * t151 + (t121 * t159 + t139 * t145 + t141 * t144) * t152) * t129;
	t103 = t209 * t210 + (t178 * t191 + (t184 * t195 + t147 * t151 + (t144 * t154 + (t121 * t171 - t139 * t192) * t165) * t152) * t168) * t129;
	t1 = [0, t103, t105, 0, 0, 0; 0, (-t108 * t207 + t117 * t198) * t190 + ((t149 * t168 - t180 * t191) * t117 + (-t206 + t211) * t108 + (t198 * t106 + (-t103 * t139 - t111 * t121 + (-t168 * t192 + t171 * t191) * t165 + (-t111 * t158 - t154 * t168) * t109) * t203 + (-t154 * t191 - t103 * t158 - t111 * t144 + t147 * t168 + (t111 * t139 - t168 * t195) * t109) * t204) * t118) * t112, (-t107 * t207 - t117 * t143) * t190 + (t107 * t211 + t124 * t117 + (-t143 * t106 - t107 * t123 + (-t105 * t139 - t110 * t121 + t145 + (-t110 * t158 - t141) * t109) * t203 + (-t105 * t158 - t110 * t144 - t122 + (t110 * t139 - t159) * t109) * t204) * t118) * t112, 0, 0, 0; 0, 0.2e1 * (-t136 * t179 + t170 * t199) * t208 + (t187 * t212 + t136 * t148 + (qJD(3) * t150 * t168 - t124 * t179 + t200 * t212) * t137) * t127, t137 * t188 * t208 + (t188 * t205 + (-t123 * t180 - t149 * t181) * t137) * t127, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:30
	% DurationCPUTime: 0.93s
	% Computational Cost: add. (1361->69), mult. (4214->164), div. (442->14), fcn. (5379->13), ass. (0->80)
	t235 = -qJD(3) + qJD(5);
	t192 = sin(qJ(3));
	t195 = cos(qJ(3));
	t187 = sin(pkin(11));
	t189 = cos(pkin(11));
	t196 = cos(qJ(2));
	t190 = cos(pkin(6));
	t193 = sin(qJ(2));
	t221 = t190 * t193;
	t203 = t187 * t221 - t189 * t196;
	t188 = sin(pkin(6));
	t223 = t187 * t188;
	t162 = t192 * t223 - t195 * t203;
	t191 = sin(qJ(5));
	t194 = cos(qJ(5));
	t205 = t192 * t203 + t195 * t223;
	t155 = t162 * t194 - t191 * t205;
	t220 = t190 * t196;
	t204 = t187 * t220 + t189 * t193;
	t171 = t204 * qJD(2);
	t156 = qJD(3) * t162 - t171 * t192;
	t157 = qJD(3) * t205 - t171 * t195;
	t133 = qJD(5) * t155 - t156 * t194 + t157 * t191;
	t150 = 0.1e1 / t155 ^ 2;
	t210 = -t162 * t191 - t194 * t205;
	t227 = t150 * t210;
	t234 = t133 * t227;
	t233 = t235 * t195;
	t232 = t235 * t192;
	t134 = qJD(5) * t210 + t156 * t191 + t157 * t194;
	t175 = t187 * t193 - t189 * t220;
	t222 = t188 * t196;
	t167 = atan2(t175, t222);
	t163 = sin(t167);
	t164 = cos(t167);
	t146 = t163 * t175 + t222 * t164;
	t143 = 0.1e1 / t146;
	t149 = 0.1e1 / t155;
	t184 = 0.1e1 / t196;
	t144 = 0.1e1 / t146 ^ 2;
	t185 = 0.1e1 / t196 ^ 2;
	t148 = t210 ^ 2;
	t137 = t148 * t150 + 0.1e1;
	t151 = t149 * t150;
	t229 = t134 * t151;
	t231 = (-t148 * t229 - t234) / t137 ^ 2;
	t176 = t187 * t196 + t189 * t221;
	t170 = t176 * qJD(2);
	t224 = t185 * t193;
	t212 = t175 * t224;
	t173 = t175 ^ 2;
	t183 = 0.1e1 / t188 ^ 2;
	t168 = t173 * t183 * t185 + 0.1e1;
	t165 = 0.1e1 / t168;
	t182 = 0.1e1 / t188;
	t225 = t165 * t182;
	t138 = (qJD(2) * t212 + t170 * t184) * t225;
	t207 = -t163 * t222 + t164 * t175;
	t213 = t164 * t188 * t193;
	t131 = -qJD(2) * t213 + t138 * t207 + t163 * t170;
	t230 = t131 * t143 * t144;
	t228 = t144 * t204;
	t208 = t191 * t192 + t194 * t195;
	t159 = t208 * t204;
	t226 = t150 * t159;
	t215 = 0.2e1 * t231;
	t214 = -0.2e1 * t151 * t210;
	t209 = t191 * t195 - t192 * t194;
	t206 = t176 * t184 + t212;
	t186 = t184 * t185;
	t174 = t204 ^ 2;
	t172 = t203 * qJD(2);
	t169 = qJD(2) * t175;
	t158 = t209 * t204;
	t142 = t174 * t144 + 0.1e1;
	t139 = t206 * t225;
	t135 = 0.1e1 / t137;
	t132 = t139 * t207 + t163 * t176 - t213;
	t130 = (-0.2e1 * t206 / t168 ^ 2 * (qJD(2) * t173 * t186 * t193 + t170 * t175 * t185) * t183 + (t170 * t224 - t169 * t184 + (t176 * t224 + (0.2e1 * t186 * t193 ^ 2 + t184) * t175) * qJD(2)) * t165) * t182;
	t1 = [0, t130, 0, 0, 0, 0; 0, 0.2e1 * (-t132 * t228 - t143 * t203) / t142 ^ 2 * (-t172 * t228 - t174 * t230) + (t171 * t143 + (-t131 * t203 - t132 * t172) * t144 - (0.2e1 * t132 * t230 + (-(-qJD(2) * t222 + t130 * t175 + t139 * t170 + (-t139 * t222 + t176) * t138) * t164 - (-t138 * t139 * t175 - t169 + (-t130 * t196 + (qJD(2) * t139 + t138) * t193) * t188) * t163) * t144) * t204) / t142, 0, 0, 0, 0; 0, (t149 * t158 + t210 * t226) * t215 + (t133 * t226 + (t158 * t150 - t159 * t214) * t134 + (t149 * t209 + t208 * t227) * t172 - ((t233 * t149 + t232 * t227) * t194 + (t232 * t149 - t233 * t227) * t191) * t204) * t135, (t149 * t155 + t210 * t227) * t215 + (0.2e1 * t234 + (t150 * t155 - t210 * t214 - t149) * t134) * t135, 0, -0.2e1 * t231 - 0.2e1 * (t133 * t150 * t135 - (-t135 * t229 - t150 * t231) * t210) * t210, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:30
	% EndTime: 2019-10-09 22:31:34
	% DurationCPUTime: 4.14s
	% Computational Cost: add. (10809->132), mult. (30690->255), div. (801->12), fcn. (39581->15), ass. (0->122)
	t351 = sin(qJ(3));
	t354 = cos(qJ(3));
	t355 = cos(qJ(2));
	t417 = cos(pkin(11));
	t418 = cos(pkin(6));
	t386 = t418 * t417;
	t416 = sin(pkin(11));
	t419 = sin(qJ(2));
	t370 = -t416 * t355 - t419 * t386;
	t348 = sin(pkin(6));
	t393 = t348 * t417;
	t329 = -t351 * t393 - t370 * t354;
	t350 = sin(qJ(5));
	t353 = cos(qJ(5));
	t364 = -t370 * t351 + t354 * t393;
	t308 = t329 * t353 + t364 * t350;
	t366 = qJD(3) * t370;
	t340 = t355 * t386 - t416 * t419;
	t424 = qJD(2) * t340 - qJD(3) * t393;
	t315 = t351 * t366 + t424 * t354;
	t362 = t424 * t351 - t354 * t366;
	t275 = t308 * qJD(5) + t315 * t350 - t362 * t353;
	t442 = -t329 * t350 + t364 * t353;
	t276 = t442 * qJD(5) + t315 * t353 + t362 * t350;
	t303 = t442 ^ 2;
	t396 = t351 * t419;
	t343 = t348 * t396 - t418 * t354;
	t395 = t354 * t419;
	t344 = t348 * t395 + t418 * t351;
	t380 = t343 * t353 - t344 * t350;
	t322 = 0.1e1 / t380 ^ 2;
	t289 = t303 * t322 + 0.1e1;
	t287 = 0.1e1 / t289;
	t326 = t343 * t350 + t344 * t353;
	t403 = t348 * t355;
	t394 = qJD(2) * t403;
	t332 = t344 * qJD(3) + t351 * t394;
	t333 = -t343 * qJD(3) + t354 * t394;
	t292 = t326 * qJD(5) - t332 * t353 + t333 * t350;
	t293 = t380 * qJD(5) + t332 * t350 + t333 * t353;
	t321 = 0.1e1 / t380;
	t406 = t442 * t322;
	t375 = t308 * t321 - t326 * t406;
	t408 = t292 * t321 * t322;
	t389 = -0.2e1 * t442 * t408;
	t422 = -0.2e1 * (-t275 * t406 + t303 * t408) / t289 ^ 2;
	t252 = (t326 * t389 + t276 * t321 + (t275 * t326 + t292 * t308 - t293 * t442) * t322) * t287 + t375 * t422;
	t260 = (t275 * t321 - t292 * t406) * t287;
	t290 = atan2(t442, -t380);
	t285 = sin(t290);
	t286 = cos(t290);
	t384 = t285 * t380 + t286 * t442;
	t255 = t384 * t260 - t285 * t275 + t286 * t292;
	t444 = t375 * t287;
	t257 = -t285 * t308 + t286 * t326 + t384 * t444;
	t271 = t285 * t442 - t286 * t380;
	t269 = 0.1e1 / t271 ^ 2;
	t385 = t418 * t416;
	t368 = -t417 * t355 + t419 * t385;
	t392 = t348 * t416;
	t331 = t351 * t392 - t354 * t368;
	t371 = t351 * t368 + t354 * t392;
	t382 = -t331 * t350 - t353 * t371;
	t304 = t382 ^ 2;
	t267 = t269 * t304 + 0.1e1;
	t265 = 0.1e1 / t267;
	t268 = 0.1e1 / t271;
	t312 = t331 * t353 - t350 * t371;
	t369 = t355 * t385 + t417 * t419;
	t338 = t369 * qJD(2);
	t316 = t331 * qJD(3) - t338 * t351;
	t317 = t371 * qJD(3) - t338 * t354;
	t278 = t312 * qJD(5) - t316 * t353 + t317 * t350;
	t279 = t382 * qJD(5) + t316 * t350 + t317 * t353;
	t415 = t255 * t268 * t269;
	t391 = -0.2e1 * t382 * t415;
	t414 = t269 * t382;
	t399 = 0.2e1 * (-t278 * t414 - t304 * t415) / t267 ^ 2;
	t409 = t286 * t382;
	t410 = t285 * t382;
	t456 = ((t312 * t255 + t257 * t278 - (t252 * t442 - t444 * t275 + t293 + (t380 * t444 - t308) * t260) * t409 - (t252 * t380 - (t442 * t444 + t326) * t260 - t444 * t292 - t276) * t410) * t269 - t257 * t391 - t279 * t268) * t265 + (t257 * t414 + t268 * t312) * t399;
	t349 = sin(qJ(6));
	t352 = cos(qJ(6));
	t299 = t312 * t352 - t349 * t369;
	t339 = t368 * qJD(2);
	t272 = t299 * qJD(6) + t279 * t349 - t339 * t352;
	t298 = t312 * t349 + t352 * t369;
	t400 = qJD(6) * t298;
	t273 = t279 * t352 + t339 * t349 - t400;
	t294 = t298 ^ 2;
	t296 = 0.1e1 / t299 ^ 2;
	t282 = t294 * t296 + 0.1e1;
	t280 = 0.1e1 / t282;
	t295 = 0.1e1 / t299;
	t407 = t296 * t298;
	t377 = -t295 * t349 + t352 * t407;
	t412 = t273 * t295 * t296;
	t437 = 0.2e1 * t298;
	t390 = t412 * t437;
	t443 = (t377 * t278 + (((-t273 + t400) * t349 - t272 * t352) * t296 + (qJD(6) * t295 + t390) * t352) * t382) * t280;
	t430 = t377 * t382;
	t379 = t350 * t354 - t351 * t353;
	t425 = t379 * t369;
	t423 = qJD(5) - qJD(3);
	t411 = 0.1e1 / t282 ^ 2 * (t272 * t407 - t294 * t412);
	t398 = -0.2e1 * t411;
	t397 = 0.2e1 * t411;
	t378 = t350 * t351 + t353 * t354;
	t320 = t378 * t369;
	t302 = -t320 * t352 + t349 * t368;
	t301 = -t320 * t349 - t352 * t368;
	t318 = t379 * t340;
	t334 = t379 * t403;
	t374 = t318 * t321 - t334 * t406;
	t365 = t423 * t378;
	t300 = ((-t350 * t395 + t353 * t396) * qJD(2) + t365 * t355) * t348;
	t284 = t378 * t339 + t423 * t425;
	t283 = t379 * t370 * qJD(2) + t365 * t340;
	t264 = t374 * t287;
	t258 = t384 * t264 - t285 * t318 + t286 * t334;
	t254 = t374 * t422 + (t334 * t389 + t283 * t321 + (t275 * t334 + t292 * t318 - t300 * t442) * t322) * t287;
	t1 = [0, t254, -t252, 0, t252, 0; 0, (-t258 * t414 + t268 * t425) * t399 + (t258 * t391 + (t425 * t255 - t258 * t278 + (t254 * t442 - t264 * t275 + t300 + (t264 * t380 - t318) * t260) * t409 + (t254 * t380 - t264 * t292 - t283 + (-t264 * t442 - t334) * t260) * t410) * t269 + (t379 * t339 - t365 * t369) * t268) * t265, t456, 0, -t456, 0; 0, (-t295 * t301 + t302 * t407) * t397 + ((t302 * qJD(6) + t284 * t349 - t338 * t352) * t295 + t302 * t390 + (-t301 * t273 - (-t301 * qJD(6) + t284 * t352 + t338 * t349) * t298 - t302 * t272) * t296) * t280, -t397 * t430 - t443, 0, -t398 * t430 + t443, t398 + (t272 * t296 * t280 + (-t280 * t412 - t296 * t411) * t298) * t437;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end