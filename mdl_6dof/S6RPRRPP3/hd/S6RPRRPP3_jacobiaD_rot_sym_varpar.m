% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP3
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
%   Wie in S6RPRRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:11
	% EndTime: 2019-10-10 01:13:12
	% DurationCPUTime: 1.05s
	% Computational Cost: add. (1976->91), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->94)
	t123 = qJ(1) + pkin(9);
	t121 = sin(t123);
	t119 = t121 ^ 2;
	t130 = sin(qJ(3));
	t125 = t130 ^ 2;
	t132 = cos(qJ(3));
	t127 = 0.1e1 / t132 ^ 2;
	t175 = t125 * t127;
	t115 = t119 * t175 + 0.1e1;
	t124 = t130 * t125;
	t126 = 0.1e1 / t132;
	t140 = qJD(3) * (t124 * t127 + t130) * t126;
	t122 = cos(t123);
	t171 = qJD(1) * t122;
	t180 = t121 * t125;
	t145 = t171 * t180;
	t187 = 0.1e1 / t115 ^ 2 * (t119 * t140 + t127 * t145);
	t198 = -0.2e1 * t187;
	t129 = sin(qJ(4));
	t131 = cos(qJ(4));
	t173 = t131 * t132;
	t109 = t121 * t129 + t122 * t173;
	t103 = 0.1e1 / t109;
	t104 = 0.1e1 / t109 ^ 2;
	t174 = t129 * t132;
	t108 = -t121 * t131 + t122 * t174;
	t183 = t104 * t108;
	t142 = -t103 * t129 + t131 * t183;
	t102 = t108 ^ 2;
	t101 = t102 * t104 + 0.1e1;
	t98 = 0.1e1 / t101;
	t197 = t142 * t98;
	t169 = qJD(3) * t121;
	t113 = 0.1e1 / t115;
	t156 = t127 * t169;
	t170 = qJD(1) * t130;
	t157 = t122 * t170;
	t167 = qJD(3) * t132;
	t87 = (-(-t121 * t167 - t157) * t126 + t125 * t156) * t113;
	t149 = t87 - t169;
	t152 = 0.1e1 + t175;
	t196 = t121 * t152;
	t150 = -t121 * t87 + qJD(3);
	t148 = qJD(4) * t132 - qJD(1);
	t168 = qJD(3) * t130;
	t195 = t148 * t129 + t131 * t168;
	t178 = t121 * t130;
	t112 = atan2(-t178, -t132);
	t111 = cos(t112);
	t110 = sin(t112);
	t161 = t110 * t178;
	t96 = -t111 * t132 - t161;
	t93 = 0.1e1 / t96;
	t94 = 0.1e1 / t96 ^ 2;
	t194 = t113 - 0.1e1;
	t120 = t122 ^ 2;
	t162 = t94 * t167;
	t181 = t111 * t130;
	t82 = t150 * t181 + (t149 * t132 - t157) * t110;
	t192 = t82 * t93 * t94;
	t92 = t120 * t125 * t94 + 0.1e1;
	t193 = (-t94 * t145 + (-t125 * t192 + t130 * t162) * t120) / t92 ^ 2;
	t147 = -qJD(1) * t132 + qJD(4);
	t143 = t147 * t131;
	t89 = t121 * t143 - t195 * t122;
	t188 = t103 * t104 * t89;
	t141 = t121 * t174 + t122 * t131;
	t155 = t129 * t168;
	t88 = t141 * qJD(1) - t109 * qJD(4) + t122 * t155;
	t191 = (-t102 * t188 - t88 * t183) / t101 ^ 2;
	t90 = 0.1e1 / t92;
	t190 = t90 * t94;
	t189 = t93 * t90;
	t185 = t122 * t94;
	t182 = t110 * t132;
	t177 = t122 * t129;
	t176 = t122 * t130;
	t172 = qJD(1) * t121;
	t166 = 0.2e1 * t192;
	t165 = -0.2e1 * t191;
	t164 = t93 * t193;
	t163 = t108 * t188;
	t160 = t113 * t125 * t126;
	t158 = t121 * t170;
	t153 = 0.2e1 * t94 * t193;
	t151 = t126 * t198;
	t146 = t121 * t160;
	t144 = t152 * t122;
	t107 = -t121 * t173 + t177;
	t100 = t113 * t196;
	t86 = (t194 * t130 * t110 - t111 * t146) * t122;
	t85 = -t121 * t182 + t181 + (-t111 * t178 + t182) * t100;
	t83 = t196 * t198 + (qJD(1) * t144 + 0.2e1 * t121 * t140) * t113;
	t1 = [t151 * t176 + (qJD(3) * t144 - t126 * t158) * t113, 0, t83, 0, 0, 0; (-t167 * t189 + (0.2e1 * t164 + (qJD(1) * t86 + t82) * t190) * t130) * t121 + (t86 * t153 * t130 + (-t86 * t162 + (t86 * t166 + ((0.2e1 * t130 * t187 - t87 * t146 - t194 * t167) * t110 + (t151 * t180 + t130 * t87 + (t124 * t156 - (t87 - 0.2e1 * t169) * t130) * t113) * t111) * t185) * t130 + (-t93 + (-(t119 - t120) * t111 * t160 + t194 * t161) * t94) * t170) * t90) * t122, 0, (-t172 * t189 + (-0.2e1 * t164 + (-qJD(3) * t85 - t82) * t190) * t122) * t132 + (t85 * t122 * t153 + (-t122 * qJD(3) * t93 - ((-t100 * t171 - t121 * t83) * t111 + (-t150 * t100 - t149) * t110) * t94 * t176 + (t122 * t166 + t94 * t172) * t85 - ((t83 - t171) * t110 + (t149 * t100 + t150) * t111) * t132 * t185) * t90) * t130, 0, 0, 0; 0.2e1 * (t103 * t141 + t107 * t183) * t191 + (0.2e1 * t107 * t163 + (t141 * t89 + t107 * t88 + (-t195 * t121 - t122 * t143) * t108) * t104 + (t147 * t177 + (-t148 * t131 + t155) * t121) * t103) * t98, 0, -t158 * t197 + (t167 * t197 + (t142 * t165 + ((-qJD(4) * t103 - 0.2e1 * t163) * t131 + (-t131 * t88 + (-qJD(4) * t108 + t89) * t129) * t104) * t98) * t130) * t122, t165 + 0.2e1 * (-t104 * t88 * t98 + (-t104 * t191 - t98 * t188) * t108) * t108, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:11
	% EndTime: 2019-10-10 01:13:12
	% DurationCPUTime: 1.51s
	% Computational Cost: add. (3939->118), mult. (5988->261), div. (1156->14), fcn. (7606->9), ass. (0->107)
	t146 = qJ(1) + pkin(9);
	t144 = sin(t146);
	t145 = cos(t146);
	t156 = cos(qJ(4));
	t154 = sin(qJ(4));
	t157 = cos(qJ(3));
	t208 = t154 * t157;
	t126 = t144 * t208 + t145 * t156;
	t155 = sin(qJ(3));
	t207 = t155 * t154;
	t121 = atan2(-t126, t207);
	t117 = sin(t121);
	t118 = cos(t121);
	t123 = t126 ^ 2;
	t148 = 0.1e1 / t154 ^ 2;
	t151 = 0.1e1 / t155 ^ 2;
	t211 = t148 * t151;
	t122 = t123 * t211 + 0.1e1;
	t119 = 0.1e1 / t122;
	t147 = 0.1e1 / t154;
	t209 = t151 * t157;
	t187 = t147 * t209;
	t171 = t126 * t187 + t144;
	t105 = t171 * t119;
	t205 = -t105 + t144;
	t229 = t205 * t155 * t117 + t118 * t157;
	t150 = 0.1e1 / t155;
	t152 = t150 * t151;
	t227 = qJD(3) * (0.2e1 * t152 * t157 ^ 2 + t150);
	t199 = qJD(4) * t156;
	t179 = t145 * t199;
	t200 = qJD(4) * t154;
	t180 = t144 * t200;
	t202 = qJD(3) * t155;
	t183 = t154 * t202;
	t110 = t126 * qJD(1) + t145 * t183 - t157 * t179 - t180;
	t214 = t144 * t156;
	t129 = t145 * t208 - t214;
	t201 = qJD(3) * t157;
	t185 = t151 * t201;
	t212 = t147 * t150;
	t226 = (t148 * t150 * t199 + t147 * t185) * t129 + t110 * t212;
	t222 = t117 * t126;
	t109 = t118 * t207 - t222;
	t106 = 0.1e1 / t109;
	t141 = 0.1e1 / t145;
	t107 = 0.1e1 / t109 ^ 2;
	t142 = 0.1e1 / t145 ^ 2;
	t168 = t154 * t201 + t155 * t199;
	t189 = t126 * t211;
	t203 = qJD(1) * t145;
	t204 = qJD(1) * t144;
	t112 = -t144 * t183 - t145 * t200 - t156 * t204 + (t144 * t199 + t154 * t203) * t157;
	t191 = t112 * t212;
	t98 = (t168 * t189 - t191) * t119;
	t166 = -t126 * t98 + t168;
	t172 = -t98 * t207 - t112;
	t94 = t172 * t117 + t166 * t118;
	t225 = t106 * t107 * t94;
	t149 = t147 * t148;
	t181 = t151 * t199;
	t184 = t152 * t201;
	t224 = (t112 * t189 + (-t148 * t184 - t149 * t181) * t123) / t122 ^ 2;
	t223 = t107 * t129;
	t221 = t117 * t129;
	t219 = t118 * t126;
	t218 = t118 * t129;
	t216 = t142 * t144;
	t215 = t142 * t151;
	t213 = t145 * t106;
	t210 = t148 * t156;
	t206 = t156 * t157;
	t124 = t129 ^ 2;
	t104 = t107 * t124 + 0.1e1;
	t198 = 0.2e1 / t104 ^ 2 * (-t110 * t223 - t124 * t225);
	t197 = 0.2e1 * t225;
	t182 = t156 * t202;
	t111 = (-qJD(1) * t157 + qJD(4)) * t214 + (-t182 + (-qJD(4) * t157 + qJD(1)) * t154) * t145;
	t130 = t144 * t154 + t145 * t206;
	t125 = t130 ^ 2;
	t116 = t125 * t215 + 0.1e1;
	t143 = t141 * t142;
	t186 = t151 * t204;
	t196 = 0.2e1 / t116 ^ 2 * (t130 * t111 * t215 + (-t142 * t184 + t143 * t186) * t125);
	t195 = -0.2e1 * t224;
	t194 = t150 * t224;
	t193 = t107 * t221;
	t190 = t126 * t212;
	t188 = t141 * t209;
	t178 = t106 * t198;
	t177 = t107 * t198;
	t176 = t129 * t197;
	t175 = t150 * t196;
	t173 = t147 * t194;
	t128 = t144 * t206 - t145 * t154;
	t170 = t126 * t210 - t128 * t147;
	t169 = t128 * t141 - t130 * t216;
	t114 = 0.1e1 / t116;
	t113 = qJD(1) * t130 - t144 * t182 - t157 * t180 - t179;
	t102 = 0.1e1 / t104;
	t101 = t170 * t150 * t119;
	t97 = (-t117 + (t118 * t190 + t117) * t119) * t129;
	t96 = -t105 * t219 + t229 * t154;
	t95 = t118 * t155 * t156 - t117 * t128 + (-t117 * t207 - t219) * t101;
	t93 = t171 * t195 + (t112 * t187 + t203 + (-t148 * t157 * t181 - t147 * t227) * t126) * t119;
	t91 = -0.2e1 * t170 * t194 + (-t170 * t185 + (t112 * t210 - t113 * t147 + (t128 * t210 + (-0.2e1 * t149 * t156 ^ 2 - t147) * t126) * qJD(4)) * t150) * t119;
	t1 = [t226 * t119 + 0.2e1 * t129 * t173, 0, t93, t91, 0, 0; t126 * t178 + (-t112 * t106 + (t110 * t97 + t126 * t94) * t107) * t102 + (t97 * t177 + (t97 * t197 + (t110 * t119 - t110 - (-t119 * t98 * t190 + t195) * t129) * t107 * t117 + (-(-0.2e1 * t126 * t173 - t98) * t223 + (-(t98 + t191) * t129 + t226 * t126) * t107 * t119) * t118) * t102) * t129, 0, t96 * t129 * t177 + (t96 * t176 + (-(-t93 * t219 + (-t112 * t118 + t222 * t98) * t105) * t129 + t96 * t110) * t107 + (-t155 * t213 - t229 * t223) * t199) * t102 + (t145 * t178 * t155 + ((-qJD(3) * t213 - (t205 * qJD(3) - t98) * t193) * t157 + (t106 * t204 + (t145 * t94 - (-t93 + t203) * t221 - (t205 * t98 - qJD(3)) * t218) * t107) * t155) * t102) * t154, (-t106 * t130 + t95 * t223) * t198 + (t95 * t176 + t111 * t106 - (-t113 + (-t154 * t91 - t156 * t98) * t155 - t166 * t101) * t193 + (t95 * t110 - t130 * t94 - (t172 * t101 - t126 * t91 - t128 * t98 - t155 * t200 + t156 * t201) * t218) * t107) * t102, 0, 0; t169 * t175 + (t169 * t185 + (t111 * t216 - t113 * t141 + (-t128 * t216 + (0.2e1 * t143 * t144 ^ 2 + t141) * t130) * qJD(1)) * t150) * t114, 0, (t130 * t188 + t156) * t196 + (-t111 * t188 + t200 + (-t142 * t157 * t186 + t141 * t227) * t130) * t114, t129 * t141 * t175 + (t110 * t141 * t150 + (-t142 * t150 * t204 + t141 * t185) * t129) * t114, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:11
	% EndTime: 2019-10-10 01:13:12
	% DurationCPUTime: 1.61s
	% Computational Cost: add. (3939->115), mult. (5988->257), div. (1156->14), fcn. (7606->9), ass. (0->107)
	t141 = qJ(1) + pkin(9);
	t139 = sin(t141);
	t140 = cos(t141);
	t149 = sin(qJ(4));
	t151 = cos(qJ(4));
	t152 = cos(qJ(3));
	t203 = t151 * t152;
	t126 = t139 * t203 - t140 * t149;
	t150 = sin(qJ(3));
	t204 = t150 * t151;
	t120 = atan2(-t126, t204);
	t116 = sin(t120);
	t117 = cos(t120);
	t122 = t126 ^ 2;
	t143 = 0.1e1 / t150 ^ 2;
	t146 = 0.1e1 / t151 ^ 2;
	t208 = t143 * t146;
	t121 = t122 * t208 + 0.1e1;
	t118 = 0.1e1 / t121;
	t145 = 0.1e1 / t151;
	t207 = t143 * t152;
	t182 = t145 * t207;
	t166 = t126 * t182 + t139;
	t104 = t166 * t118;
	t226 = t104 - t139;
	t227 = t116 * t150 * t226 - t117 * t152;
	t194 = qJD(4) * t152;
	t176 = t149 * t194;
	t198 = qJD(3) * t150;
	t225 = t151 * t198 + t176;
	t142 = 0.1e1 / t150;
	t144 = t142 * t143;
	t224 = qJD(3) * (0.2e1 * t144 * t152 ^ 2 + t142);
	t199 = qJD(1) * t152;
	t181 = t139 * t199;
	t195 = qJD(4) * t151;
	t200 = qJD(1) * t140;
	t110 = -t139 * t195 + t225 * t140 - t149 * t200 + t151 * t181;
	t129 = t139 * t149 + t140 * t203;
	t197 = qJD(3) * t152;
	t180 = t143 * t197;
	t196 = qJD(4) * t149;
	t209 = t142 * t145;
	t223 = t129 * (-t142 * t146 * t196 + t145 * t180) + t110 * t209;
	t222 = t151 * (-qJD(1) + t194) - t149 * t198;
	t218 = t116 * t126;
	t108 = t117 * t204 - t218;
	t105 = 0.1e1 / t108;
	t136 = 0.1e1 / t140;
	t106 = 0.1e1 / t108 ^ 2;
	t137 = 0.1e1 / t140 ^ 2;
	t163 = -t150 * t196 + t151 * t197;
	t184 = t126 * t208;
	t112 = qJD(1) * t129 - t225 * t139 - t140 * t195;
	t186 = t112 * t209;
	t97 = (t163 * t184 - t186) * t118;
	t161 = -t126 * t97 + t163;
	t168 = -t204 * t97 - t112;
	t93 = t116 * t168 + t117 * t161;
	t221 = t105 * t106 * t93;
	t147 = t145 * t146;
	t179 = t144 * t197;
	t220 = 0.1e1 / t121 ^ 2 * (t112 * t184 + (t143 * t147 * t196 - t146 * t179) * t122);
	t219 = t106 * t129;
	t217 = t116 * t129;
	t215 = t117 * t126;
	t214 = t117 * t129;
	t212 = t137 * t139;
	t211 = t137 * t143;
	t210 = t140 * t105;
	t206 = t146 * t149;
	t205 = t149 * t152;
	t201 = qJD(1) * t139;
	t124 = t129 ^ 2;
	t103 = t124 * t106 + 0.1e1;
	t193 = 0.2e1 / t103 ^ 2 * (-t110 * t219 - t124 * t221);
	t192 = 0.2e1 * t221;
	t167 = (-qJD(4) + t199) * t149;
	t109 = t139 * t167 - t222 * t140;
	t128 = t139 * t151 - t140 * t205;
	t123 = t128 ^ 2;
	t115 = t123 * t211 + 0.1e1;
	t138 = t136 * t137;
	t191 = 0.2e1 / t115 ^ 2 * (t128 * t109 * t211 + (t138 * t143 * t201 - t137 * t179) * t123);
	t190 = -0.2e1 * t220;
	t189 = t142 * t220;
	t188 = t106 * t217;
	t185 = t126 * t209;
	t183 = t136 * t207;
	t175 = t105 * t193;
	t174 = t106 * t193;
	t173 = t142 * t191;
	t172 = t129 * t192;
	t170 = t145 * t189;
	t125 = t139 * t205 + t140 * t151;
	t165 = -t125 * t145 + t126 * t206;
	t164 = -t125 * t136 - t128 * t212;
	t113 = 0.1e1 / t115;
	t111 = t222 * t139 + t140 * t167;
	t101 = 0.1e1 / t103;
	t100 = t165 * t142 * t118;
	t96 = (-t116 + (t117 * t185 + t116) * t118) * t129;
	t95 = -t104 * t215 - t227 * t151;
	t94 = -t117 * t150 * t149 + t116 * t125 - (-t116 * t204 - t215) * t100;
	t92 = t166 * t190 + (t112 * t182 + t200 + (-t145 * t224 + t176 * t208) * t126) * t118;
	t90 = 0.2e1 * t165 * t189 + (t165 * t180 + (-t112 * t206 + t111 * t145 + (t125 * t206 + (-0.2e1 * t147 * t149 ^ 2 - t145) * t126) * qJD(4)) * t142) * t118;
	t1 = [t223 * t118 + 0.2e1 * t129 * t170, 0, t92, t90, 0, 0; t126 * t175 + (-t112 * t105 + (t110 * t96 + t126 * t93) * t106) * t101 + (t96 * t174 + (t96 * t192 + (t110 * t118 - t110 - (-t118 * t185 * t97 + t190) * t129) * t106 * t116 + (-(-0.2e1 * t126 * t170 - t97) * t219 + (-(t97 + t186) * t129 + t223 * t126) * t106 * t118) * t117) * t101) * t129, 0, t95 * t129 * t174 + (t95 * t172 + (-(-t215 * t92 + (-t112 * t117 + t218 * t97) * t104) * t129 + t95 * t110) * t106 + (t150 * t210 - t227 * t219) * t196) * t101 + (t140 * t175 * t150 + ((-qJD(3) * t210 - (-qJD(3) * t226 - t97) * t188) * t152 + (t105 * t201 + (t140 * t93 - (-t92 + t200) * t217 - (-t226 * t97 - qJD(3)) * t214) * t106) * t150) * t101) * t151, (-t105 * t128 + t219 * t94) * t193 + (t94 * t172 + t109 * t105 - (t111 + (t149 * t97 - t151 * t90) * t150 + t161 * t100) * t188 + (t94 * t110 - t128 * t93 - (-t100 * t168 + t125 * t97 - t126 * t90 - t149 * t197 - t150 * t195) * t214) * t106) * t101, 0, 0; t164 * t173 + (t164 * t180 + (t109 * t212 + t111 * t136 + (t125 * t212 + (0.2e1 * t138 * t139 ^ 2 + t136) * t128) * qJD(1)) * t142) * t113, 0, (t128 * t183 - t149) * t191 + (-t109 * t183 + t195 + (t136 * t224 - t181 * t211) * t128) * t113, t129 * t136 * t173 + (t110 * t136 * t142 + (-t137 * t142 * t201 + t136 * t180) * t129) * t113, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end