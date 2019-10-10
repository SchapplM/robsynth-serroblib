% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR10
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
%   Wie in S6RPRRPR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR10_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR10_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:31
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (813->90), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->96)
	t132 = cos(qJ(1));
	t195 = 0.2e1 * t132;
	t167 = qJD(3) * t132;
	t126 = t132 ^ 2;
	t128 = sin(qJ(3));
	t121 = 0.1e1 / t128 ^ 2;
	t131 = cos(qJ(3));
	t125 = t131 ^ 2;
	t179 = t121 * t125;
	t118 = t126 * t179 + 0.1e1;
	t116 = 0.1e1 / t118;
	t120 = 0.1e1 / t128;
	t158 = t121 * t167;
	t129 = sin(qJ(1));
	t171 = qJD(1) * t131;
	t160 = t129 * t171;
	t90 = ((t128 * t167 + t160) * t120 + t125 * t158) * t116;
	t152 = -t90 + t167;
	t173 = t132 * t131;
	t115 = atan2(-t173, t128);
	t113 = sin(t115);
	t114 = cos(t115);
	t99 = -t113 * t173 + t114 * t128;
	t96 = 0.1e1 / t99;
	t127 = sin(qJ(4));
	t175 = t132 * t127;
	t130 = cos(qJ(4));
	t177 = t129 * t130;
	t110 = t128 * t177 + t175;
	t106 = 0.1e1 / t110;
	t107 = 0.1e1 / t110 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t194 = t116 - 0.1e1;
	t123 = t129 ^ 2;
	t170 = qJD(1) * t132;
	t147 = t125 * t129 * t170;
	t169 = qJD(3) * t128;
	t163 = t97 * t169;
	t153 = -t132 * t90 + qJD(3);
	t181 = t114 * t131;
	t85 = t153 * t181 + (t152 * t128 + t160) * t113;
	t192 = t85 * t96 * t97;
	t95 = t123 * t125 * t97 + 0.1e1;
	t193 = (t97 * t147 + (-t125 * t192 - t131 * t163) * t123) / t95 ^ 2;
	t93 = 0.1e1 / t95;
	t191 = t93 * t97;
	t190 = t96 * t93;
	t174 = t132 * t130;
	t178 = t129 * t127;
	t109 = t128 * t178 - t174;
	t105 = t109 ^ 2;
	t104 = t105 * t107 + 0.1e1;
	t184 = t107 * t109;
	t150 = qJD(1) * t128 + qJD(4);
	t143 = t150 * t132;
	t151 = qJD(4) * t128 + qJD(1);
	t145 = t151 * t127;
	t168 = qJD(3) * t131;
	t92 = t130 * t143 + (t130 * t168 - t145) * t129;
	t188 = t106 * t107 * t92;
	t144 = t151 * t130;
	t91 = t129 * t144 + (t129 * t168 + t143) * t127;
	t189 = 0.1e1 / t104 ^ 2 * (-t105 * t188 + t91 * t184);
	t187 = t129 * t97;
	t124 = t131 * t125;
	t141 = qJD(3) * (-t121 * t124 - t131) * t120;
	t186 = (-t121 * t147 + t126 * t141) / t118 ^ 2;
	t185 = t106 * t127;
	t183 = t109 * t130;
	t182 = t113 * t132;
	t180 = t120 * t125;
	t176 = t129 * t131;
	t172 = qJD(1) * t129;
	t166 = -0.2e1 * t192;
	t165 = 0.2e1 * t189;
	t164 = t96 * t193;
	t162 = t131 * t186;
	t161 = t116 * t180;
	t159 = t131 * t170;
	t157 = -0.2e1 * t97 * t193;
	t156 = 0.1e1 + t179;
	t155 = 0.2e1 * t109 * t188;
	t154 = t186 * t195;
	t149 = t132 * t161;
	t148 = t194 * t131 * t113;
	t146 = t156 * t129;
	t142 = t107 * t183 - t185;
	t140 = -t150 * t129 + t131 * t167;
	t112 = t128 * t174 - t178;
	t111 = t128 * t175 + t177;
	t103 = t156 * t132 * t116;
	t101 = 0.1e1 / t104;
	t89 = (-t114 * t149 - t148) * t129;
	t88 = t128 * t182 + t181 + (-t113 * t128 - t114 * t173) * t103;
	t86 = -t156 * t154 + (-qJD(1) * t146 + t141 * t195) * t116;
	t1 = [-0.2e1 * t129 * t120 * t162 + (-qJD(3) * t146 + t120 * t159) * t116, 0, t86, 0, 0, 0; (t169 * t190 + (0.2e1 * t164 + (qJD(1) * t89 + t85) * t191) * t131) * t132 + (t89 * t157 * t131 + (-t89 * t163 + (t89 * t166 + ((t90 * t149 + t194 * t169 + 0.2e1 * t162) * t113 + (t154 * t180 + t131 * t90 + (t124 * t158 + (-t90 + 0.2e1 * t167) * t131) * t116) * t114) * t187) * t131 + (t96 + ((t123 - t126) * t114 * t161 - t132 * t148) * t97) * t171) * t93) * t129, 0, (t170 * t190 + (-0.2e1 * t164 + (-qJD(3) * t88 - t85) * t191) * t129) * t128 + (t88 * t129 * t157 + (t129 * qJD(3) * t96 + (-t114 * t132 * t86 + t152 * t113 + (-qJD(3) * t113 + t114 * t172 + t182 * t90) * t103) * t97 * t176 + (t129 * t166 + t97 * t170) * t88 + ((-t86 - t172) * t113 + (t152 * t103 - t153) * t114) * t128 * t187) * t93) * t131, 0, 0, 0; (-t106 * t111 + t112 * t184) * t165 + (t112 * t155 + t132 * t106 * t144 + t140 * t185 + (t132 * t109 * t145 - t111 * t92 - t112 * t91 - t140 * t183) * t107) * t101, 0, t142 * t165 * t176 + (-t142 * t159 + (t142 * t169 + ((qJD(4) * t106 + t155) * t130 + (-t130 * t91 + (qJD(4) * t109 - t92) * t127) * t107) * t131) * t129) * t101, -0.2e1 * t189 + 0.2e1 * (t91 * t107 * t101 + (-t101 * t188 - t107 * t189) * t109) * t109, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:32
	% DurationCPUTime: 1.38s
	% Computational Cost: add. (1824->120), mult. (6168->265), div. (1114->15), fcn. (7752->9), ass. (0->113)
	t153 = sin(qJ(1));
	t155 = cos(qJ(3));
	t208 = t153 * t155;
	t152 = sin(qJ(3));
	t151 = sin(qJ(4));
	t228 = cos(qJ(1));
	t190 = t228 * t151;
	t154 = cos(qJ(4));
	t209 = t153 * t154;
	t135 = t152 * t190 + t209;
	t207 = t155 * t151;
	t125 = atan2(t135, t207);
	t120 = cos(t125);
	t119 = sin(t125);
	t219 = t119 * t135;
	t114 = t120 * t207 + t219;
	t111 = 0.1e1 / t114;
	t134 = t152 * t209 + t190;
	t129 = 0.1e1 / t134;
	t143 = 0.1e1 / t151;
	t148 = 0.1e1 / t155;
	t112 = 0.1e1 / t114 ^ 2;
	t130 = 0.1e1 / t134 ^ 2;
	t144 = 0.1e1 / t151 ^ 2;
	t189 = t228 * t154;
	t210 = t153 * t151;
	t133 = t152 * t210 - t189;
	t128 = t133 ^ 2;
	t109 = t112 * t128 + 0.1e1;
	t182 = qJD(1) * t228;
	t172 = t152 * t182;
	t165 = t228 * qJD(4) + t172;
	t177 = qJD(4) * t152 + qJD(1);
	t203 = qJD(3) * t155;
	t186 = t153 * t203;
	t117 = t177 * t209 + (t165 + t186) * t151;
	t222 = t117 * t112;
	t132 = t135 ^ 2;
	t149 = 0.1e1 / t155 ^ 2;
	t213 = t144 * t149;
	t127 = t132 * t213 + 0.1e1;
	t123 = 0.1e1 / t127;
	t201 = qJD(4) * t155;
	t205 = qJD(3) * t152;
	t166 = -t151 * t205 + t154 * t201;
	t192 = t135 * t213;
	t188 = t228 * t155;
	t173 = qJD(3) * t188;
	t174 = t154 * t182;
	t175 = t152 * t189;
	t115 = -qJD(4) * t175 - t151 * t173 - t174 + (qJD(1) * t152 + qJD(4)) * t210;
	t214 = t143 * t148;
	t195 = t115 * t214;
	t103 = (-t166 * t192 - t195) * t123;
	t164 = -t103 * t135 - t166;
	t99 = (-t103 * t207 - t115) * t119 - t164 * t120;
	t226 = t111 * t112 * t99;
	t227 = 0.1e1 / t109 ^ 2 * (-t128 * t226 + t133 * t222);
	t146 = t153 ^ 2;
	t147 = t155 ^ 2;
	t211 = t146 * t147;
	t194 = t130 * t211;
	t126 = 0.1e1 + t194;
	t185 = t154 * t203;
	t118 = t165 * t154 + (-t177 * t151 + t185) * t153;
	t221 = t118 * t129 * t130;
	t176 = t211 * t221;
	t225 = (-t176 + (-t146 * t152 * t203 + t147 * t153 * t182) * t130) / t126 ^ 2;
	t145 = t143 * t144;
	t150 = t148 / t147;
	t202 = qJD(4) * t154;
	t184 = t149 * t202;
	t224 = (-t115 * t192 + (t144 * t150 * t205 - t145 * t184) * t132) / t127 ^ 2;
	t223 = t112 * t133;
	t220 = t119 * t133;
	t218 = t119 * t155;
	t217 = t120 * t133;
	t216 = t120 * t135;
	t215 = t120 * t152;
	t212 = t144 * t154;
	t206 = qJD(1) * t153;
	t204 = qJD(3) * t153;
	t200 = 0.2e1 * t227;
	t199 = 0.2e1 * t226;
	t198 = 0.2e1 * t225;
	t197 = -0.2e1 * t224;
	t196 = t112 * t220;
	t193 = t135 * t214;
	t191 = t143 * t149 * t152;
	t187 = t149 * t205;
	t168 = t135 * t191 + t228;
	t110 = t168 * t123;
	t183 = t228 - t110;
	t181 = -0.2e1 * t111 * t227;
	t180 = t112 * t200;
	t179 = 0.2e1 * t148 * t224;
	t178 = -0.2e1 * t133 * t208;
	t171 = t143 * t179;
	t170 = t133 * t199 - t222;
	t136 = t175 - t210;
	t169 = t135 * t212 - t136 * t143;
	t167 = t130 * t136 * t153 - t228 * t129;
	t163 = t117 * t214 - (t144 * t148 * t202 - t143 * t187) * t133;
	t121 = 0.1e1 / t126;
	t116 = t134 * qJD(1) + qJD(4) * t135 - t154 * t173;
	t107 = 0.1e1 / t109;
	t106 = t169 * t148 * t123;
	t102 = (-t119 + (-t120 * t193 + t119) * t123) * t133;
	t101 = t110 * t216 + (t183 * t218 - t215) * t151;
	t100 = t120 * t154 * t155 + t119 * t136 - (-t119 * t207 + t216) * t106;
	t98 = t168 * t197 + (-t115 * t191 - t206 + (-t144 * t152 * t184 + (0.2e1 * t150 * t152 ^ 2 + t148) * t143 * qJD(3)) * t135) * t123;
	t96 = t169 * t179 + (-t169 * t187 + (t115 * t212 - t116 * t143 + (-t136 * t212 + (0.2e1 * t145 * t154 ^ 2 + t143) * t135) * qJD(4)) * t148) * t123;
	t1 = [-t163 * t123 + t133 * t171, 0, t98, t96, 0, 0; t135 * t181 + (-t115 * t111 + (-t102 * t117 - t135 * t99) * t112) * t107 + (t102 * t180 + (t102 * t199 + (-t117 * t123 + t117 - (t103 * t123 * t193 + t197) * t133) * t112 * t119 + (-(t135 * t171 - t103) * t223 + (-(t103 + t195) * t133 + t163 * t135) * t112 * t123) * t120) * t107) * t133, 0, t101 * t133 * t180 + (-(t98 * t216 + (-t103 * t219 - t115 * t120) * t110) * t223 + t170 * t101 + (t111 * t208 - (-t110 * t218 + t119 * t188 - t215) * t223) * t202) * t107 + (t181 * t208 + ((-t111 * t204 - (-t183 * qJD(3) + t103) * t196) * t152 + (t111 * t182 + (-t153 * t99 - (-t98 - t206) * t220 - (t183 * t103 - qJD(3)) * t217) * t112) * t155) * t107) * t151, (t100 * t223 - t111 * t134) * t200 + (t118 * t111 + t170 * t100 - (-t116 + (-t103 * t154 - t151 * t96) * t155 - t164 * t106) * t196 + (-t134 * t99 - (-t154 * t205 - t151 * t201 + t106 * t115 + t135 * t96 + (t106 * t207 + t136) * t103) * t217) * t112) * t107, 0, 0; t167 * t155 * t198 + (t167 * t205 + ((-qJD(1) * t129 + 0.2e1 * t136 * t221) * t153 + (t116 * t153 - t228 * t118 - t136 * t182) * t130) * t155) * t121, 0, (t129 * t152 * t153 + t154 * t194) * t198 + (0.2e1 * t154 * t176 + (-t172 - t186) * t129 + ((t118 * t152 - 0.2e1 * t147 * t174) * t153 + (qJD(4) * t147 * t151 + 0.2e1 * t152 * t185) * t146) * t130) * t121, t130 * t178 * t225 + (t178 * t221 + (t117 * t208 + (-t152 * t204 + t155 * t182) * t133) * t130) * t121, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:32
	% DurationCPUTime: 1.49s
	% Computational Cost: add. (1454->115), mult. (4276->250), div. (493->12), fcn. (5073->11), ass. (0->111)
	t264 = -qJD(4) + qJD(6);
	t186 = sin(qJ(3));
	t178 = 0.1e1 / t186 ^ 2;
	t190 = cos(qJ(3));
	t182 = t190 ^ 2;
	t244 = t178 * t182;
	t263 = t190 * t244;
	t191 = cos(qJ(1));
	t217 = 0.1e1 + t244;
	t262 = t191 * t217;
	t183 = t191 ^ 2;
	t175 = t183 * t244 + 0.1e1;
	t173 = 0.1e1 / t175;
	t177 = 0.1e1 / t186;
	t187 = sin(qJ(1));
	t235 = qJD(1) * t190;
	t219 = t187 * t235;
	t231 = qJD(3) * t191;
	t139 = (-(-t186 * t231 - t219) * t177 + t231 * t244) * t173;
	t261 = t139 - t231;
	t189 = cos(qJ(4));
	t260 = t264 * t189;
	t185 = sin(qJ(4));
	t259 = t264 * t185;
	t213 = qJD(1) * t186 + qJD(4);
	t201 = t213 * t191;
	t214 = qJD(4) * t186 + qJD(1);
	t203 = t189 * t214;
	t232 = qJD(3) * t190;
	t142 = t187 * t203 + (t187 * t232 + t201) * t185;
	t202 = t214 * t185;
	t143 = t189 * t201 + (t189 * t232 - t202) * t187;
	t184 = sin(qJ(6));
	t188 = cos(qJ(6));
	t238 = t191 * t189;
	t242 = t187 * t185;
	t166 = t186 * t242 - t238;
	t239 = t191 * t185;
	t241 = t187 * t189;
	t167 = t186 * t241 + t239;
	t207 = t166 * t188 - t167 * t184;
	t133 = t207 * qJD(6) + t142 * t184 + t143 * t188;
	t154 = t166 * t184 + t167 * t188;
	t148 = 0.1e1 / t154;
	t204 = t184 * t185 + t188 * t189;
	t205 = t184 * t189 - t185 * t188;
	t149 = 0.1e1 / t154 ^ 2;
	t251 = t149 * t207;
	t258 = t205 * t148 + t204 * t251;
	t237 = t191 * t190;
	t172 = atan2(t237, -t186);
	t170 = sin(t172);
	t171 = cos(t172);
	t160 = t170 * t237 - t171 * t186;
	t157 = 0.1e1 / t160;
	t158 = 0.1e1 / t160 ^ 2;
	t257 = t173 - 0.1e1;
	t132 = t154 * qJD(6) - t142 * t188 + t143 * t184;
	t147 = t207 ^ 2;
	t138 = t147 * t149 + 0.1e1;
	t150 = t148 * t149;
	t253 = t133 * t150;
	t256 = (-t132 * t251 - t147 * t253) / t138 ^ 2;
	t180 = t187 ^ 2;
	t243 = t180 * t182;
	t146 = t158 * t243 + 0.1e1;
	t234 = qJD(1) * t191;
	t210 = t182 * t187 * t234;
	t246 = t171 * t190;
	t130 = (t139 * t191 - qJD(3)) * t246 + (t261 * t186 - t219) * t170;
	t254 = t130 * t157 * t158;
	t255 = (-t243 * t254 + (-t180 * t186 * t232 + t210) * t158) / t146 ^ 2;
	t252 = t139 * t190;
	t240 = t187 * t190;
	t164 = t204 * t240;
	t250 = t149 * t164;
	t249 = t158 * t187;
	t248 = t158 * t190;
	t247 = t170 * t191;
	t245 = t177 * t182;
	t236 = qJD(1) * t187;
	t233 = qJD(3) * t186;
	t226 = 0.2e1 * t256;
	t225 = 0.2e1 * t254;
	t224 = -0.2e1 * t150 * t207;
	t200 = qJD(3) * (-t190 - t263) * t177;
	t223 = -0.2e1 * (-t178 * t210 + t183 * t200) / t175 ^ 2;
	t222 = t190 * t255;
	t221 = t158 * t240;
	t220 = t191 * t245;
	t218 = t190 * t234;
	t216 = t133 * t224;
	t215 = t190 * t223;
	t212 = t173 * t220;
	t211 = t257 * t190 * t170;
	t209 = t217 * t187;
	t168 = t186 * t239 + t241;
	t169 = t186 * t238 - t242;
	t206 = t168 * t188 - t169 * t184;
	t156 = t168 * t184 + t169 * t188;
	t199 = -t213 * t187 + t190 * t231;
	t163 = t205 * t240;
	t162 = t173 * t262;
	t144 = 0.1e1 / t146;
	t141 = t199 * t189 - t191 * t202;
	t140 = t199 * t185 + t191 * t203;
	t136 = 0.1e1 / t138;
	t135 = (t171 * t212 + t211) * t187;
	t134 = -t186 * t247 - t246 + (t170 * t186 + t171 * t237) * t162;
	t131 = t223 * t262 + (-qJD(1) * t209 + 0.2e1 * t191 * t200) * t173;
	t1 = [t187 * t177 * t215 + (-qJD(3) * t209 + t177 * t218) * t173, 0, t131, 0, 0, 0; (-0.2e1 * t157 * t222 + (-t157 * t233 + (-qJD(1) * t135 - t130) * t248) * t144) * t191 + (0.2e1 * t158 * t222 * t135 + (-((-t139 * t212 - t257 * t233 + t215) * t170 + (t220 * t223 - t252 + (t252 + (-0.2e1 * t190 - t263) * t231) * t173) * t171) * t221 + (t158 * t233 + t190 * t225) * t135 + (-t157 + ((t180 - t183) * t173 * t171 * t245 - t191 * t211) * t158) * t235) * t144) * t187, 0, 0.2e1 * (t134 * t248 + t157 * t186) * t187 * t255 + ((-t157 * t234 + (qJD(3) * t134 + t130) * t249) * t186 + (-t187 * qJD(3) * t157 - (t131 * t171 * t191 + t261 * t170 + (qJD(3) * t170 - t139 * t247 - t171 * t236) * t162) * t221 + (-t158 * t234 + t187 * t225) * t134 - ((t131 + t236) * t170 + ((-t162 * t191 + 0.1e1) * qJD(3) + (t162 - t191) * t139) * t171) * t186 * t249) * t190) * t144, 0, 0, 0; (t148 * t206 - t156 * t251) * t226 + ((t156 * qJD(6) - t140 * t188 + t141 * t184) * t148 + t156 * t216 + (t206 * t133 + (t206 * qJD(6) + t140 * t184 + t141 * t188) * t207 - t156 * t132) * t149) * t136, 0, (-t148 * t163 - t207 * t250) * t226 + (-t132 * t250 + (-t163 * t149 + t164 * t224) * t133 + t258 * t218 + (-t258 * t233 + ((t260 * t148 + t259 * t251) * t188 + (t259 * t148 - t260 * t251) * t184) * t190) * t187) * t136, (t148 * t154 + t207 * t251) * t226 + (-t133 * t148 - t207 * t216 + (0.2e1 * t207 * t132 + t154 * t133) * t149) * t136, 0, -0.2e1 * t256 - 0.2e1 * (t132 * t149 * t136 - (-t136 * t253 - t149 * t256) * t207) * t207;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end