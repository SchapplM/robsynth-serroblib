% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP8
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
%   Wie in S6RPRRPP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:50
	% DurationCPUTime: 0.93s
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
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:50
	% DurationCPUTime: 1.44s
	% Computational Cost: add. (1786->114), mult. (5988->255), div. (1156->14), fcn. (7606->9), ass. (0->107)
	t150 = sin(qJ(3));
	t151 = sin(qJ(1));
	t152 = cos(qJ(4));
	t202 = t151 * t152;
	t149 = sin(qJ(4));
	t154 = cos(qJ(1));
	t205 = t149 * t154;
	t131 = t150 * t205 + t202;
	t153 = cos(qJ(3));
	t201 = t153 * t149;
	t123 = atan2(t131, t201);
	t117 = sin(t123);
	t118 = cos(t123);
	t128 = t131 ^ 2;
	t140 = 0.1e1 / t149 ^ 2;
	t147 = 0.1e1 / t153 ^ 2;
	t210 = t140 * t147;
	t125 = t128 * t210 + 0.1e1;
	t121 = 0.1e1 / t125;
	t139 = 0.1e1 / t149;
	t206 = t147 * t150;
	t184 = t139 * t206;
	t168 = t131 * t184 + t154;
	t108 = t168 * t121;
	t199 = t108 - t154;
	t226 = -t199 * t153 * t117 - t118 * t150;
	t146 = 0.1e1 / t153;
	t171 = qJD(1) * t150 + qJD(4);
	t169 = t171 * t154;
	t172 = qJD(4) * t150 + qJD(1);
	t195 = qJD(3) * t153;
	t116 = t152 * t169 + (-t149 * t172 + t152 * t195) * t151;
	t130 = t150 * t202 + t205;
	t127 = t130 ^ 2;
	t144 = 0.1e1 / t151 ^ 2;
	t208 = t144 * t147;
	t124 = t127 * t208 + 0.1e1;
	t143 = 0.1e1 / t151;
	t145 = t143 * t144;
	t148 = t146 * t147;
	t196 = qJD(3) * t150;
	t180 = t148 * t196;
	t197 = qJD(1) * t154;
	t182 = t147 * t197;
	t220 = (t130 * t116 * t208 + (t144 * t180 - t145 * t182) * t127) / t124 ^ 2;
	t225 = -0.2e1 * t146 * t220;
	t223 = qJD(3) * (0.2e1 * t148 * t150 ^ 2 + t146);
	t216 = t117 * t131;
	t112 = t118 * t201 + t216;
	t109 = 0.1e1 / t112;
	t110 = 0.1e1 / t112 ^ 2;
	t200 = t154 * t152;
	t203 = t151 * t149;
	t129 = t150 * t203 - t200;
	t126 = t129 ^ 2;
	t107 = t110 * t126 + 0.1e1;
	t115 = t172 * t202 + (t151 * t195 + t169) * t149;
	t218 = t110 * t129;
	t193 = qJD(4) * t153;
	t165 = -t149 * t196 + t152 * t193;
	t185 = t131 * t210;
	t194 = qJD(4) * t152;
	t178 = t150 * t194;
	t179 = t154 * t195;
	t113 = -t149 * t179 - t152 * t197 - t154 * t178 + t171 * t203;
	t211 = t139 * t146;
	t187 = t113 * t211;
	t101 = (-t165 * t185 - t187) * t121;
	t164 = -t101 * t131 - t165;
	t97 = (-t101 * t201 - t113) * t117 - t164 * t118;
	t221 = t109 * t110 * t97;
	t222 = 0.1e1 / t107 ^ 2 * (t115 * t218 - t126 * t221);
	t141 = t139 * t140;
	t219 = (-t113 * t185 + (-t141 * t147 * t194 + t140 * t180) * t128) / t125 ^ 2;
	t217 = t117 * t129;
	t214 = t118 * t129;
	t213 = t118 * t131;
	t209 = t140 * t152;
	t207 = t144 * t154;
	t204 = t151 * t109;
	t198 = qJD(1) * t151;
	t192 = 0.2e1 * t222;
	t191 = 0.2e1 * t221;
	t190 = -0.2e1 * t219;
	t188 = t110 * t217;
	t186 = t131 * t211;
	t183 = t143 * t206;
	t181 = t147 * t196;
	t177 = -0.2e1 * t109 * t222;
	t176 = t110 * t192;
	t175 = t129 * t191;
	t174 = 0.2e1 * t146 * t219;
	t170 = t139 * t174;
	t132 = t150 * t200 - t203;
	t167 = t130 * t207 - t132 * t143;
	t166 = t131 * t209 - t132 * t139;
	t163 = t115 * t211 - (t140 * t146 * t194 - t139 * t181) * t129;
	t119 = 0.1e1 / t124;
	t114 = qJD(1) * t130 + qJD(4) * t131 - t152 * t179;
	t105 = 0.1e1 / t107;
	t104 = t166 * t146 * t121;
	t100 = (-t117 + (-t118 * t186 + t117) * t121) * t129;
	t99 = t108 * t213 + t149 * t226;
	t98 = t118 * t153 * t152 + t117 * t132 - (-t117 * t201 + t213) * t104;
	t96 = t168 * t190 + (-t113 * t184 - t198 + (t139 * t223 - t178 * t210) * t131) * t121;
	t94 = t166 * t174 + (-t166 * t181 + (t113 * t209 - t114 * t139 + (-t132 * t209 + (0.2e1 * t141 * t152 ^ 2 + t139) * t131) * qJD(4)) * t146) * t121;
	t1 = [-t121 * t163 + t129 * t170, 0, t96, t94, 0, 0; t131 * t177 + (-t113 * t109 + (-t100 * t115 - t131 * t97) * t110) * t105 + (t100 * t176 + (t100 * t191 + (-t115 * t121 + t115 - (t101 * t121 * t186 + t190) * t129) * t110 * t117 + (-(t131 * t170 - t101) * t218 + (-(t101 + t187) * t129 + t163 * t131) * t110 * t121) * t118) * t105) * t129, 0, t99 * t129 * t176 + (t99 * t175 + (-(t96 * t213 + (-t101 * t216 - t113 * t118) * t108) * t129 - t99 * t115) * t110 + (t153 * t204 - t218 * t226) * t194) * t105 + (t151 * t177 * t153 + ((-qJD(3) * t204 - (qJD(3) * t199 + t101) * t188) * t150 + (t109 * t197 + (-t151 * t97 - (-t96 - t198) * t217 - (-t101 * t199 - qJD(3)) * t214) * t110) * t153) * t105) * t149, (-t109 * t130 + t218 * t98) * t192 + (t98 * t175 + t116 * t109 - (-t114 + (-t101 * t152 - t149 * t94) * t153 - t164 * t104) * t188 + (-t98 * t115 - t130 * t97 - (-t152 * t196 - t149 * t193 + t104 * t113 + t131 * t94 + (t104 * t201 + t132) * t101) * t214) * t110) * t105, 0, 0; t167 * t225 + (t167 * t181 + (t116 * t207 + t114 * t143 + (t132 * t207 + (-0.2e1 * t145 * t154 ^ 2 - t143) * t130) * qJD(1)) * t146) * t119, 0, 0.2e1 * (t130 * t183 + t152) * t220 + (-t116 * t183 + qJD(4) * t149 + (t144 * t150 * t182 - t143 * t223) * t130) * t119, t129 * t143 * t225 + (t115 * t143 * t146 + (-t144 * t146 * t197 + t143 * t181) * t129) * t119, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:50
	% DurationCPUTime: 1.42s
	% Computational Cost: add. (1786->113), mult. (5988->254), div. (1156->14), fcn. (7606->9), ass. (0->108)
	t148 = sin(qJ(3));
	t150 = cos(qJ(4));
	t152 = cos(qJ(1));
	t199 = t150 * t152;
	t147 = sin(qJ(4));
	t149 = sin(qJ(1));
	t201 = t149 * t147;
	t132 = t148 * t199 - t201;
	t151 = cos(qJ(3));
	t198 = t151 * t150;
	t123 = atan2(t132, t198);
	t117 = sin(t123);
	t118 = cos(t123);
	t128 = t132 ^ 2;
	t142 = 0.1e1 / t150 ^ 2;
	t145 = 0.1e1 / t151 ^ 2;
	t206 = t142 * t145;
	t125 = t128 * t206 + 0.1e1;
	t121 = 0.1e1 / t125;
	t141 = 0.1e1 / t150;
	t204 = t145 * t148;
	t180 = t141 * t204;
	t166 = t132 * t180 + t152;
	t108 = t166 * t121;
	t197 = t108 - t152;
	t226 = t197 * t117 * t151 + t118 * t148;
	t144 = 0.1e1 / t151;
	t191 = qJD(4) * t148;
	t167 = (-qJD(1) - t191) * t150;
	t169 = qJD(1) * t148 + qJD(4);
	t193 = qJD(3) * t151;
	t221 = t149 * t193 + t152 * t169;
	t115 = -t221 * t147 + t149 * t167;
	t129 = -t148 * t201 + t199;
	t126 = t129 ^ 2;
	t139 = 0.1e1 / t149 ^ 2;
	t209 = t139 * t145;
	t124 = t126 * t209 + 0.1e1;
	t138 = 0.1e1 / t149;
	t140 = t138 * t139;
	t146 = t144 * t145;
	t194 = qJD(3) * t148;
	t177 = t146 * t194;
	t195 = qJD(1) * t152;
	t179 = t145 * t195;
	t218 = (t129 * t115 * t209 + (t139 * t177 - t140 * t179) * t126) / t124 ^ 2;
	t225 = -0.2e1 * t144 * t218;
	t200 = t149 * t150;
	t203 = t147 * t152;
	t131 = t148 * t203 + t200;
	t205 = t142 * t147;
	t164 = -t131 * t141 + t132 * t205;
	t224 = t144 * t164;
	t222 = qJD(3) * (0.2e1 * t146 * t148 ^ 2 + t144);
	t214 = t117 * t132;
	t112 = t118 * t198 + t214;
	t109 = 0.1e1 / t112;
	t110 = 0.1e1 / t112 ^ 2;
	t130 = t148 * t200 + t203;
	t127 = t130 ^ 2;
	t107 = t110 * t127 + 0.1e1;
	t174 = t147 * t191;
	t196 = qJD(1) * t149;
	t116 = -t147 * t196 - t149 * t174 + t221 * t150;
	t216 = t110 * t130;
	t190 = qJD(4) * t151;
	t163 = -t147 * t190 - t150 * t194;
	t182 = t132 * t206;
	t175 = t152 * t193;
	t114 = qJD(1) * t130 + qJD(4) * t131 - t150 * t175;
	t207 = t141 * t144;
	t184 = t114 * t207;
	t101 = (-t163 * t182 - t184) * t121;
	t162 = -t101 * t132 - t163;
	t97 = (-t198 * t101 - t114) * t117 - t162 * t118;
	t219 = t109 * t110 * t97;
	t220 = 0.1e1 / t107 ^ 2 * (t116 * t216 - t127 * t219);
	t143 = t141 * t142;
	t192 = qJD(4) * t147;
	t217 = (-t114 * t182 + (t143 * t145 * t192 + t142 * t177) * t128) / t125 ^ 2;
	t215 = t117 * t130;
	t212 = t118 * t130;
	t211 = t118 * t132;
	t208 = t139 * t152;
	t202 = t149 * t109;
	t189 = 0.2e1 * t220;
	t188 = 0.2e1 * t219;
	t187 = -0.2e1 * t217;
	t185 = t110 * t215;
	t183 = t132 * t207;
	t181 = t138 * t204;
	t178 = t145 * t194;
	t173 = -0.2e1 * t109 * t220;
	t172 = t110 * t189;
	t171 = t130 * t188;
	t168 = 0.2e1 * t207 * t217;
	t165 = t129 * t208 + t131 * t138;
	t161 = t116 * t207 - (-t142 * t144 * t192 - t141 * t178) * t130;
	t119 = 0.1e1 / t124;
	t113 = t152 * t167 + (t149 * t169 - t175) * t147;
	t105 = 0.1e1 / t107;
	t104 = t121 * t224;
	t100 = (-t117 + (-t118 * t183 + t117) * t121) * t130;
	t99 = t108 * t211 - t226 * t150;
	t98 = -t118 * t151 * t147 - t117 * t131 + (-t117 * t198 + t211) * t104;
	t96 = t166 * t187 + (-t114 * t180 - t196 + (t141 * t222 + t174 * t206) * t132) * t121;
	t94 = t187 * t224 + (t164 * t178 + (-t114 * t205 + t113 * t141 + (-t131 * t205 + (0.2e1 * t143 * t147 ^ 2 + t141) * t132) * qJD(4)) * t144) * t121;
	t1 = [-t121 * t161 + t130 * t168, 0, t96, t94, 0, 0; t132 * t173 + (-t114 * t109 + (-t100 * t116 - t132 * t97) * t110) * t105 + (t100 * t172 + (t100 * t188 + (-t116 * t121 + t116 - (t101 * t121 * t183 + t187) * t130) * t110 * t117 + (-(t132 * t168 - t101) * t216 + (-(t101 + t184) * t130 + t161 * t132) * t110 * t121) * t118) * t105) * t130, 0, t99 * t130 * t172 + (t99 * t171 + (-(t211 * t96 + (-t101 * t214 - t114 * t118) * t108) * t130 - t99 * t116) * t110 + (-t151 * t202 - t226 * t216) * t192) * t105 + (t149 * t173 * t151 + ((-qJD(3) * t202 - (qJD(3) * t197 + t101) * t185) * t148 + (t109 * t195 + (-t149 * t97 - (-t96 - t196) * t215 - (-t101 * t197 - qJD(3)) * t212) * t110) * t151) * t105) * t150, (-t109 * t129 + t216 * t98) * t189 + (t98 * t171 + t115 * t109 - (t113 + (t101 * t147 - t150 * t94) * t151 + t162 * t104) * t185 + (-t98 * t116 - t129 * t97 - (t147 * t194 - t150 * t190 - t104 * t114 + t132 * t94 + (-t104 * t198 - t131) * t101) * t212) * t110) * t105, 0, 0; t165 * t225 + (t165 * t178 + (t115 * t208 - t113 * t138 + (-t131 * t208 + (-0.2e1 * t140 * t152 ^ 2 - t138) * t129) * qJD(1)) * t144) * t119, 0, 0.2e1 * (t129 * t181 - t147) * t218 + (-t115 * t181 + qJD(4) * t150 + (t139 * t148 * t179 - t138 * t222) * t129) * t119, t130 * t138 * t225 + (t116 * t138 * t144 + (-t139 * t144 * t195 + t138 * t178) * t130) * t119, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end