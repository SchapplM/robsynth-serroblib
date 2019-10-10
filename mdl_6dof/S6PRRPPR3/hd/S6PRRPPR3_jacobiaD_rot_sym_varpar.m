% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR3
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
%   Wie in S6PRRPPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
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
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:55
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
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:56
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (2420->88), mult. (7430->188), div. (524->12), fcn. (9601->11), ass. (0->89)
	t164 = sin(pkin(10));
	t166 = cos(pkin(10));
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
	t144 = qJD(3) * t159 + t168 * t185;
	t151 = 0.1e1 / t158;
	t202 = t139 * t152;
	t109 = (-t121 * t151 + t144 * t202) * t129;
	t133 = atan2(-t139, t158);
	t125 = sin(t133);
	t126 = cos(t133);
	t183 = -t125 * t158 - t126 * t139;
	t106 = t109 * t183 - t125 * t121 + t126 * t144;
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
	t124 = qJD(3) * t181 - t148 * t170;
	t205 = t124 * t136 * t137;
	t187 = t150 * t205;
	t149 = t179 * qJD(2);
	t200 = t149 * t180;
	t208 = (-t137 * t200 - t187) / t131 ^ 2;
	t207 = t118 * t181;
	t123 = qJD(3) * t143 - t148 * t168;
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
	t145 = -qJD(3) * t158 + t170 * t185;
	t127 = 0.1e1 / t131;
	t122 = -qJD(3) * t139 + t146 * t170;
	t112 = 0.1e1 / t115;
	t111 = t129 * t210;
	t110 = t182 * t129;
	t108 = (-t125 * t154 + t126 * t195) * t168 + t183 * t111;
	t107 = t110 * t183 - t125 * t141 + t126 * t159;
	t105 = t182 * t209 + (t159 * t184 - t122 * t151 + (t121 * t159 + t139 * t145 + t141 * t144) * t152) * t129;
	t103 = t209 * t210 + (t178 * t191 + (t184 * t195 + t147 * t151 + (t144 * t154 + (t121 * t171 - t139 * t192) * t165) * t152) * t168) * t129;
	t1 = [0, t103, t105, 0, 0, 0; 0, (-t108 * t207 + t117 * t198) * t190 + ((t149 * t168 - t180 * t191) * t117 + (-t206 + t211) * t108 + (t198 * t106 + (-t103 * t139 - t111 * t121 + (-t168 * t192 + t171 * t191) * t165 + (-t111 * t158 - t154 * t168) * t109) * t203 + (-t154 * t191 - t103 * t158 - t111 * t144 + t147 * t168 + (t111 * t139 - t168 * t195) * t109) * t204) * t118) * t112, (-t107 * t207 - t117 * t143) * t190 + (t107 * t211 + t124 * t117 + (-t143 * t106 - t107 * t123 + (-t105 * t139 - t110 * t121 + t145 + (-t110 * t158 - t141) * t109) * t203 + (-t105 * t158 - t110 * t144 - t122 + (t110 * t139 - t159) * t109) * t204) * t118) * t112, 0, 0, 0; 0, 0.2e1 * (-t136 * t179 + t170 * t199) * t208 + (t187 * t212 + t136 * t148 + (qJD(3) * t150 * t168 - t124 * t179 + t200 * t212) * t137) * t127, t137 * t188 * t208 + (t188 * t205 + (-t123 * t180 - t149 * t181) * t137) * t127, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:55
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (694->54), mult. (2307->132), div. (427->14), fcn. (2998->11), ass. (0->67)
	t145 = sin(pkin(10));
	t147 = cos(pkin(10));
	t150 = sin(qJ(2));
	t148 = cos(pkin(6));
	t152 = cos(qJ(2));
	t170 = t148 * t152;
	t133 = t145 * t150 - t147 * t170;
	t146 = sin(pkin(6));
	t172 = t146 * t152;
	t125 = atan2(t133, t172);
	t121 = sin(t125);
	t122 = cos(t125);
	t108 = t121 * t133 + t122 * t172;
	t105 = 0.1e1 / t108;
	t149 = sin(qJ(3));
	t151 = cos(qJ(3));
	t171 = t148 * t150;
	t159 = t145 * t171 - t147 * t152;
	t173 = t145 * t146;
	t161 = t149 * t159 + t151 * t173;
	t115 = 0.1e1 / t161;
	t142 = 0.1e1 / t152;
	t106 = 0.1e1 / t108 ^ 2;
	t116 = 0.1e1 / t161 ^ 2;
	t143 = 0.1e1 / t152 ^ 2;
	t134 = t145 * t152 + t147 * t171;
	t128 = t134 * qJD(2);
	t174 = t143 * t150;
	t165 = t133 * t174;
	t131 = t133 ^ 2;
	t141 = 0.1e1 / t146 ^ 2;
	t126 = t131 * t141 * t143 + 0.1e1;
	t123 = 0.1e1 / t126;
	t140 = 0.1e1 / t146;
	t175 = t123 * t140;
	t100 = (qJD(2) * t165 + t128 * t142) * t175;
	t163 = -t121 * t172 + t122 * t133;
	t166 = t122 * t146 * t150;
	t97 = -qJD(2) * t166 + t163 * t100 + t121 * t128;
	t180 = t105 * t106 * t97;
	t160 = t145 * t170 + t147 * t150;
	t179 = t106 * t160;
	t120 = t149 * t173 - t151 * t159;
	t118 = t120 ^ 2;
	t178 = t116 * t118;
	t177 = t116 * t120;
	t117 = t115 * t116;
	t176 = t117 * t118;
	t169 = qJD(3) * t120;
	t112 = 0.1e1 + t178;
	t129 = t160 * qJD(2);
	t113 = -t129 * t149 + t169;
	t114 = t161 * qJD(3) - t129 * t151;
	t167 = t114 * t177;
	t168 = 0.2e1 / t112 ^ 2 * (t113 * t176 + t167);
	t164 = t115 * t151 + t149 * t177;
	t162 = t134 * t142 + t165;
	t144 = t142 * t143;
	t132 = t160 ^ 2;
	t130 = t159 * qJD(2);
	t127 = t133 * qJD(2);
	t110 = 0.1e1 / t112;
	t104 = t106 * t132 + 0.1e1;
	t101 = t162 * t175;
	t98 = t163 * t101 + t121 * t134 - t166;
	t96 = (-0.2e1 * t162 / t126 ^ 2 * (qJD(2) * t131 * t144 * t150 + t128 * t133 * t143) * t141 + (t128 * t174 - t127 * t142 + (t134 * t174 + (0.2e1 * t144 * t150 ^ 2 + t142) * t133) * qJD(2)) * t123) * t140;
	t1 = [0, t96, 0, 0, 0, 0; 0, 0.2e1 * (-t105 * t159 - t98 * t179) * (-t130 * t179 - t132 * t180) / t104 ^ 2 + (t129 * t105 + (-t98 * t130 - t159 * t97) * t106 - (0.2e1 * t98 * t180 + (-(-qJD(2) * t172 + t101 * t128 + t133 * t96 + (-t101 * t172 + t134) * t100) * t122 - (-t100 * t101 * t133 - t127 + (-t152 * t96 + (qJD(2) * t101 + t100) * t150) * t146) * t121) * t106) * t160) / t104, 0, 0, 0, 0; 0, -t164 * t160 * t168 + (-t164 * t130 - ((-0.2e1 * t113 * t117 * t120 + qJD(3) * t115) * t149 + (-t114 * t149 + (-t113 - t169) * t151) * t116) * t160) * t110, (t115 * t161 + t178) * t168 + (-0.2e1 * t167 + (-t116 * t161 + t115 - 0.2e1 * t176) * t113) * t110, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:55
	% EndTime: 2019-10-09 22:10:56
	% DurationCPUTime: 1.41s
	% Computational Cost: add. (3002->107), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->104)
	t211 = cos(pkin(6));
	t217 = cos(qJ(2));
	t267 = cos(pkin(10));
	t238 = t267 * t217;
	t209 = sin(pkin(10));
	t214 = sin(qJ(2));
	t252 = t209 * t214;
	t200 = t211 * t238 - t252;
	t193 = t200 * qJD(2);
	t216 = cos(qJ(3));
	t213 = sin(qJ(3));
	t239 = t267 * t214;
	t251 = t209 * t217;
	t226 = -t211 * t239 - t251;
	t210 = sin(pkin(6));
	t240 = t210 * t267;
	t269 = t226 * t213 - t216 * t240;
	t165 = t269 * qJD(3) + t193 * t216;
	t185 = -t213 * t240 - t226 * t216;
	t182 = t185 ^ 2;
	t249 = t210 * t216;
	t204 = t211 * t213 + t214 * t249;
	t198 = 0.1e1 / t204 ^ 2;
	t178 = t182 * t198 + 0.1e1;
	t176 = 0.1e1 / t178;
	t250 = t210 * t213;
	t203 = t211 * t216 - t214 * t250;
	t248 = t210 * t217;
	t241 = qJD(2) * t248;
	t190 = t203 * qJD(3) + t216 * t241;
	t197 = 0.1e1 / t204;
	t257 = t185 * t198;
	t148 = (-t165 * t197 + t190 * t257) * t176;
	t179 = atan2(-t185, t204);
	t174 = sin(t179);
	t175 = cos(t179);
	t233 = -t174 * t204 - t175 * t185;
	t144 = t233 * t148 - t174 * t165 + t175 * t190;
	t158 = -t174 * t185 + t175 * t204;
	t155 = 0.1e1 / t158;
	t156 = 0.1e1 / t158 ^ 2;
	t272 = t144 * t155 * t156;
	t228 = t211 * t252 - t238;
	t188 = t209 * t250 - t216 * t228;
	t271 = 0.2e1 * t188 * t272;
	t229 = -t197 * t200 + t248 * t257;
	t270 = t216 * t229;
	t256 = t190 * t197 * t198;
	t268 = -0.2e1 * (t165 * t257 - t182 * t256) / t178 ^ 2;
	t212 = sin(qJ(6));
	t215 = cos(qJ(6));
	t227 = t211 * t251 + t239;
	t230 = t209 * t249 + t213 * t228;
	t173 = -t212 * t227 - t215 * t230;
	t169 = 0.1e1 / t173;
	t170 = 0.1e1 / t173 ^ 2;
	t195 = t227 * qJD(2);
	t166 = t188 * qJD(3) - t195 * t213;
	t196 = t228 * qJD(2);
	t159 = t173 * qJD(6) + t166 * t212 - t196 * t215;
	t254 = t227 * t215;
	t172 = -t212 * t230 + t254;
	t168 = t172 ^ 2;
	t163 = t168 * t170 + 0.1e1;
	t261 = t170 * t172;
	t245 = qJD(6) * t172;
	t160 = t166 * t215 + t196 * t212 - t245;
	t264 = t160 * t169 * t170;
	t266 = (t159 * t261 - t168 * t264) / t163 ^ 2;
	t265 = t156 * t188;
	t167 = t230 * qJD(3) - t195 * t216;
	t263 = t167 * t156;
	t262 = t169 * t212;
	t260 = t172 * t215;
	t259 = t174 * t188;
	t258 = t175 * t188;
	t255 = t227 * t213;
	t253 = t227 * t216;
	t247 = qJD(2) * t214;
	t246 = qJD(3) * t213;
	t183 = t188 ^ 2;
	t154 = t183 * t156 + 0.1e1;
	t244 = 0.2e1 * (-t183 * t272 + t188 * t263) / t154 ^ 2;
	t243 = 0.2e1 * t266;
	t237 = 0.2e1 * t172 * t264;
	t236 = -0.2e1 * t185 * t256;
	t234 = -qJD(6) * t255 - t195;
	t232 = t170 * t260 - t262;
	t231 = -t197 * t269 + t203 * t257;
	t225 = -qJD(3) * t253 + qJD(6) * t228 + t196 * t213;
	t194 = t226 * qJD(2);
	t189 = -t204 * qJD(3) - t213 * t241;
	t181 = t212 * t228 - t213 * t254;
	t180 = -t212 * t255 - t215 * t228;
	t164 = t185 * qJD(3) + t193 * t213;
	t161 = 0.1e1 / t163;
	t151 = 0.1e1 / t154;
	t150 = t176 * t270;
	t149 = t231 * t176;
	t146 = (-t174 * t200 + t175 * t248) * t216 + t233 * t150;
	t145 = t233 * t149 - t174 * t269 + t175 * t203;
	t143 = t231 * t268 + (t203 * t236 + t164 * t197 + (t165 * t203 + t185 * t189 + t190 * t269) * t198) * t176;
	t141 = t268 * t270 + (-t229 * t246 + (t236 * t248 - t194 * t197 + (t190 * t200 + (t165 * t217 - t185 * t247) * t210) * t198) * t216) * t176;
	t1 = [0, t141, t143, 0, 0, 0; 0, (t146 * t265 + t155 * t253) * t244 + ((t196 * t216 + t227 * t246) * t155 + (-t263 + t271) * t146 + (t253 * t144 - (-t141 * t185 - t150 * t165 + (-t216 * t247 - t217 * t246) * t210 + (-t150 * t204 - t200 * t216) * t148) * t258 - (t200 * t246 - t141 * t204 - t150 * t190 - t194 * t216 + (t150 * t185 - t216 * t248) * t148) * t259) * t156) * t151, (t145 * t265 - t155 * t230) * t244 + (t145 * t271 - t166 * t155 + (-t230 * t144 - t145 * t167 - (-t143 * t185 - t149 * t165 + t189 + (-t149 * t204 - t269) * t148) * t258 - (-t143 * t204 - t149 * t190 + t164 + (t149 * t185 - t203) * t148) * t259) * t156) * t151, 0, 0, 0; 0, (-t169 * t180 + t181 * t261) * t243 + (t181 * t237 + t234 * t169 * t215 + t225 * t262 + (t234 * t172 * t212 - t181 * t159 - t180 * t160 - t225 * t260) * t170) * t161, t232 * t188 * t243 + (-t232 * t167 + ((qJD(6) * t169 + t237) * t215 + (-t159 * t215 + (-t160 + t245) * t212) * t170) * t188) * t161, 0, 0, -0.2e1 * t266 + 0.2e1 * (t159 * t170 * t161 + (-t161 * t264 - t170 * t266) * t172) * t172;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end