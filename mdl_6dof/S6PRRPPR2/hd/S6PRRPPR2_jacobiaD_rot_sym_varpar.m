% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR2
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
%   Wie in S6PRRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:07
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
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:07
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
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:07
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (4935->90), mult. (7372->193), div. (524->12), fcn. (9540->11), ass. (0->87)
	t176 = sin(pkin(10));
	t178 = cos(pkin(10));
	t180 = sin(qJ(2));
	t179 = cos(pkin(6));
	t181 = cos(qJ(2));
	t202 = t179 * t181;
	t165 = -t176 * t180 + t178 * t202;
	t158 = t165 * qJD(2);
	t203 = t179 * t180;
	t166 = t176 * t181 + t178 * t203;
	t175 = qJ(3) + pkin(11);
	t173 = sin(t175);
	t177 = sin(pkin(6));
	t206 = t177 * t178;
	t194 = t173 * t206;
	t174 = cos(t175);
	t199 = qJD(3) * t174;
	t130 = -qJD(3) * t194 + t158 * t173 + t166 * t199;
	t148 = t166 * t173 + t174 * t206;
	t145 = t148 ^ 2;
	t205 = t177 * t180;
	t156 = t173 * t205 - t174 * t179;
	t154 = 0.1e1 / t156 ^ 2;
	t138 = t145 * t154 + 0.1e1;
	t136 = 0.1e1 / t138;
	t157 = t173 * t179 + t174 * t205;
	t200 = qJD(2) * t181;
	t193 = t177 * t200;
	t143 = t157 * qJD(3) + t173 * t193;
	t153 = 0.1e1 / t156;
	t210 = t148 * t154;
	t118 = (-t130 * t153 + t143 * t210) * t136;
	t139 = atan2(-t148, t156);
	t134 = sin(t139);
	t135 = cos(t139);
	t191 = -t134 * t156 - t135 * t148;
	t115 = t191 * t118 - t130 * t134 + t135 * t143;
	t129 = -t134 * t148 + t135 * t156;
	t126 = 0.1e1 / t129;
	t127 = 0.1e1 / t129 ^ 2;
	t218 = t115 * t126 * t127;
	t195 = t176 * t203;
	t168 = t178 * t181 - t195;
	t207 = t176 * t177;
	t189 = -t168 * t173 + t174 * t207;
	t217 = -0.2e1 * t189 * t218;
	t204 = t177 * t181;
	t188 = -t153 * t165 + t204 * t210;
	t216 = t173 * t188;
	t211 = t143 * t153 * t154;
	t215 = -0.2e1 * (t130 * t210 - t145 * t211) / t138 ^ 2;
	t167 = t176 * t202 + t178 * t180;
	t162 = 0.1e1 / t167;
	t163 = 0.1e1 / t167 ^ 2;
	t214 = t127 * t189;
	t213 = t134 * t189;
	t212 = t135 * t189;
	t152 = t168 * t174 + t173 * t207;
	t209 = t152 * t168;
	t208 = t167 * t173;
	t201 = qJD(2) * t180;
	t146 = t189 ^ 2;
	t124 = t127 * t146 + 0.1e1;
	t160 = t167 * qJD(2);
	t132 = t152 * qJD(3) - t160 * t173;
	t198 = 0.2e1 * (-t132 * t214 - t146 * t218) / t124 ^ 2;
	t133 = t189 * qJD(3) - t160 * t174;
	t147 = t152 ^ 2;
	t142 = t147 * t163 + 0.1e1;
	t161 = -qJD(2) * t195 + t178 * t200;
	t164 = t162 * t163;
	t197 = 0.2e1 * (t133 * t152 * t163 - t147 * t161 * t164) / t142 ^ 2;
	t192 = -0.2e1 * t148 * t211;
	t150 = t166 * t174 - t194;
	t190 = -t150 * t153 + t157 * t210;
	t159 = t166 * qJD(2);
	t144 = -t156 * qJD(3) + t174 * t193;
	t140 = 0.1e1 / t142;
	t131 = -t148 * qJD(3) + t158 * t174;
	t121 = 0.1e1 / t124;
	t120 = t136 * t216;
	t119 = t190 * t136;
	t117 = (-t134 * t165 + t135 * t204) * t173 + t191 * t120;
	t116 = t191 * t119 - t134 * t150 + t135 * t157;
	t114 = t190 * t215 + (t157 * t192 - t131 * t153 + (t130 * t157 + t143 * t150 + t144 * t148) * t154) * t136;
	t112 = t215 * t216 + (t188 * t199 + (t192 * t204 + t153 * t159 + (t143 * t165 + (t130 * t181 - t148 * t201) * t177) * t154) * t173) * t136;
	t1 = [0, t112, t114, 0, 0, 0; 0, (-t117 * t214 + t126 * t208) * t198 + ((-t161 * t173 - t167 * t199) * t126 + t117 * t217 + (-t117 * t132 + t208 * t115 + (-t112 * t148 - t120 * t130 + (-t173 * t201 + t181 * t199) * t177 + (-t120 * t156 - t165 * t173) * t118) * t212 + (-t165 * t199 - t112 * t156 - t120 * t143 + t159 * t173 + (t120 * t148 - t173 * t204) * t118) * t213) * t127) * t121, (-t116 * t214 - t126 * t152) * t198 + (t116 * t217 + t133 * t126 + (-t152 * t115 - t116 * t132 + (-t114 * t148 - t119 * t130 + t144 + (-t119 * t156 - t150) * t118) * t212 + (-t114 * t156 - t119 * t143 - t131 + (t119 * t148 - t157) * t118) * t213) * t127) * t121, 0, 0, 0; 0, (t162 * t167 * t174 + t163 * t209) * t197 + (qJD(3) * t162 * t208 + (-t133 * t168 + t152 * t160) * t163 + (0.2e1 * t164 * t209 + (t163 * t167 - t162) * t174) * t161) * t140, -t189 * t162 * t197 + (-t161 * t163 * t189 - t132 * t162) * t140, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:07
	% EndTime: 2019-10-09 22:09:08
	% DurationCPUTime: 1.49s
	% Computational Cost: add. (5804->108), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->105)
	t228 = cos(pkin(6));
	t232 = cos(qJ(2));
	t282 = cos(pkin(10));
	t253 = t282 * t232;
	t226 = sin(pkin(10));
	t230 = sin(qJ(2));
	t266 = t226 * t230;
	t216 = t228 * t253 - t266;
	t212 = t216 * qJD(2);
	t225 = qJ(3) + pkin(11);
	t224 = cos(t225);
	t223 = sin(t225);
	t254 = t282 * t230;
	t265 = t226 * t232;
	t241 = -t228 * t254 - t265;
	t227 = sin(pkin(6));
	t255 = t227 * t282;
	t284 = t241 * t223 - t224 * t255;
	t179 = t284 * qJD(3) + t212 * t224;
	t202 = -t223 * t255 - t241 * t224;
	t199 = t202 ^ 2;
	t264 = t227 * t230;
	t211 = t228 * t223 + t224 * t264;
	t208 = 0.1e1 / t211 ^ 2;
	t192 = t199 * t208 + 0.1e1;
	t190 = 0.1e1 / t192;
	t210 = -t223 * t264 + t228 * t224;
	t263 = t227 * t232;
	t256 = qJD(2) * t263;
	t198 = t210 * qJD(3) + t224 * t256;
	t207 = 0.1e1 / t211;
	t271 = t202 * t208;
	t162 = (-t179 * t207 + t198 * t271) * t190;
	t193 = atan2(-t202, t211);
	t186 = sin(t193);
	t187 = cos(t193);
	t248 = -t186 * t211 - t187 * t202;
	t158 = t248 * t162 - t186 * t179 + t187 * t198;
	t172 = -t186 * t202 + t187 * t211;
	t169 = 0.1e1 / t172;
	t170 = 0.1e1 / t172 ^ 2;
	t288 = t158 * t169 * t170;
	t231 = cos(qJ(6));
	t218 = -t228 * t266 + t253;
	t267 = t226 * t227;
	t244 = -t218 * t223 + t224 * t267;
	t229 = sin(qJ(6));
	t242 = -t228 * t265 - t254;
	t269 = t242 * t229;
	t247 = -t231 * t244 + t269;
	t287 = t247 * qJD(6);
	t205 = t218 * t224 + t223 * t267;
	t286 = 0.2e1 * t205 * t288;
	t243 = -t207 * t216 + t263 * t271;
	t285 = t224 * t243;
	t272 = t198 * t207 * t208;
	t283 = -0.2e1 * (t179 * t271 - t199 * t272) / t192 ^ 2;
	t268 = t242 * t231;
	t189 = -t229 * t244 - t268;
	t183 = 0.1e1 / t189;
	t184 = 0.1e1 / t189 ^ 2;
	t214 = t242 * qJD(2);
	t180 = t205 * qJD(3) + t214 * t223;
	t215 = t218 * qJD(2);
	t173 = t189 * qJD(6) - t180 * t231 + t215 * t229;
	t182 = t247 ^ 2;
	t177 = t182 * t184 + 0.1e1;
	t276 = t184 * t247;
	t174 = t180 * t229 + t215 * t231 + t287;
	t279 = t174 * t183 * t184;
	t281 = (-t173 * t276 - t182 * t279) / t177 ^ 2;
	t280 = t170 * t205;
	t181 = t244 * qJD(3) + t214 * t224;
	t278 = t181 * t170;
	t277 = t183 * t231;
	t275 = t186 * t205;
	t274 = t187 * t205;
	t273 = t247 * t229;
	t270 = t242 * t224;
	t262 = qJD(2) * t230;
	t261 = qJD(3) * t223;
	t200 = t205 ^ 2;
	t168 = t200 * t170 + 0.1e1;
	t260 = 0.2e1 * (-t200 * t288 + t205 * t278) / t168 ^ 2;
	t259 = 0.2e1 * t281;
	t252 = -0.2e1 * t247 * t279;
	t251 = -0.2e1 * t202 * t272;
	t249 = qJD(6) * t223 * t242 + t214;
	t246 = -t184 * t273 + t277;
	t245 = -t207 * t284 + t210 * t271;
	t240 = -qJD(3) * t270 + qJD(6) * t218 + t215 * t223;
	t213 = t241 * qJD(2);
	t197 = -t211 * qJD(3) - t223 * t256;
	t195 = t218 * t231 + t223 * t269;
	t194 = t218 * t229 - t223 * t268;
	t178 = t202 * qJD(3) + t212 * t223;
	t175 = 0.1e1 / t177;
	t165 = 0.1e1 / t168;
	t164 = t190 * t285;
	t163 = t245 * t190;
	t160 = (-t186 * t216 + t187 * t263) * t224 + t248 * t164;
	t159 = t248 * t163 - t186 * t284 + t187 * t210;
	t157 = t245 * t283 + (t210 * t251 + t178 * t207 + (t179 * t210 + t197 * t202 + t198 * t284) * t208) * t190;
	t155 = t283 * t285 + (-t243 * t261 + (t251 * t263 - t207 * t213 + (t198 * t216 + (t179 * t232 - t202 * t262) * t227) * t208) * t224) * t190;
	t1 = [0, t155, t157, 0, 0, 0; 0, (t160 * t280 - t169 * t270) * t260 + ((-t215 * t224 - t242 * t261) * t169 + (-t278 + t286) * t160 + (-t270 * t158 - (-t155 * t202 - t164 * t179 + (-t224 * t262 - t232 * t261) * t227 + (-t164 * t211 - t216 * t224) * t162) * t274 - (t216 * t261 - t155 * t211 - t164 * t198 - t213 * t224 + (t164 * t202 - t224 * t263) * t162) * t275) * t170) * t165, (t159 * t280 - t169 * t244) * t260 + (t159 * t286 - t180 * t169 + (-t244 * t158 - t159 * t181 - (-t157 * t202 - t163 * t179 + t197 + (-t163 * t211 - t284) * t162) * t274 - (-t157 * t211 - t163 * t198 + t178 + (t163 * t202 - t210) * t162) * t275) * t170) * t165, 0, 0, 0; 0, (-t183 * t194 - t195 * t276) * t259 + (t195 * t252 + t249 * t183 * t229 + t240 * t277 + (t231 * t247 * t249 - t195 * t173 - t194 * t174 - t240 * t273) * t184) * t175, t246 * t205 * t259 + (-t246 * t181 + ((qJD(6) * t183 + t252) * t229 + (-t173 * t229 + (t174 + t287) * t231) * t184) * t205) * t175, 0, 0, -0.2e1 * t281 - 0.2e1 * (t173 * t184 * t175 - (-t175 * t279 - t184 * t281) * t247) * t247;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end