% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR1
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
%   Wie in S6RPPRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:34
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:18
	% EndTime: 2019-10-09 23:34:19
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (2921->87), mult. (2191->194), div. (456->12), fcn. (2616->9), ass. (0->91)
	t116 = qJ(1) + pkin(9);
	t112 = sin(t116);
	t106 = t112 ^ 2;
	t115 = pkin(10) + qJ(4);
	t111 = sin(t115);
	t105 = t111 ^ 2;
	t113 = cos(t115);
	t108 = 0.1e1 / t113 ^ 2;
	t156 = t105 * t108;
	t103 = t106 * t156 + 0.1e1;
	t104 = t111 * t105;
	t107 = 0.1e1 / t113;
	t155 = t107 * t111;
	t127 = qJD(4) * (t104 * t107 * t108 + t155);
	t114 = cos(t116);
	t147 = qJD(1) * t114;
	t135 = t112 * t147;
	t164 = 0.1e1 / t103 ^ 2 * (t106 * t127 + t135 * t156);
	t177 = -0.2e1 * t164;
	t101 = 0.1e1 / t103;
	t130 = 0.1e1 + t156;
	t173 = t112 * t130;
	t84 = t101 * t173;
	t176 = t112 * t84 - 0.1e1;
	t118 = cos(pkin(11));
	t146 = qJD(4) * t111;
	t132 = t118 * t146;
	t117 = sin(pkin(11));
	t150 = t114 * t117;
	t151 = t112 * t118;
	t95 = -t113 * t151 + t150;
	t89 = t95 * qJD(1) - t114 * t132;
	t149 = t114 * t118;
	t152 = t112 * t117;
	t97 = t113 * t149 + t152;
	t92 = 0.1e1 / t97 ^ 2;
	t175 = t89 * t92;
	t96 = t113 * t150 - t151;
	t165 = t92 * t96;
	t90 = t96 ^ 2;
	t87 = t90 * t92 + 0.1e1;
	t85 = 0.1e1 / t87;
	t91 = 0.1e1 / t97;
	t174 = (-t117 * t91 + t118 * t165) * t85;
	t153 = t112 * t111;
	t100 = atan2(-t153, -t113);
	t98 = sin(t100);
	t138 = t98 * t153;
	t99 = cos(t100);
	t82 = -t113 * t99 - t138;
	t79 = 0.1e1 / t82;
	t80 = 0.1e1 / t82 ^ 2;
	t110 = t114 ^ 2;
	t144 = qJD(4) * t113;
	t137 = t80 * t144;
	t145 = qJD(4) * t112;
	t160 = t113 * t98;
	t134 = t108 * t145;
	t75 = (-(-t111 * t147 - t112 * t144) * t107 + t105 * t134) * t101;
	t70 = (t75 - t145) * t160 + (-t98 * t147 + (-t112 * t75 + qJD(4)) * t99) * t111;
	t171 = t70 * t79 * t80;
	t78 = t105 * t110 * t80 + 0.1e1;
	t172 = (t110 * t111 * t137 + (-t110 * t171 - t80 * t135) * t105) / t78 ^ 2;
	t166 = t91 * t175;
	t133 = t117 * t146;
	t94 = -t113 * t152 - t149;
	t88 = t94 * qJD(1) - t114 * t133;
	t170 = (t88 * t165 - t90 * t166) / t87 ^ 2;
	t76 = 0.1e1 / t78;
	t168 = t76 * t80;
	t167 = t79 * t76;
	t163 = t111 * t98;
	t162 = t111 * t99;
	t159 = t114 * t80;
	t158 = qJD(4) * t84;
	t157 = t105 * t107;
	t154 = t111 * t114;
	t148 = qJD(1) * t112;
	t143 = qJD(4) * t114;
	t142 = 0.2e1 * t171;
	t141 = 0.2e1 * t170;
	t140 = t79 * t172;
	t139 = t96 * t166;
	t136 = t112 * t157;
	t131 = 0.2e1 * t80 * t172;
	t129 = t130 * t114;
	t126 = -t99 * t136 + t163;
	t74 = (t126 * t101 - t163) * t114;
	t72 = (-t112 + t84) * t160 - t176 * t162;
	t71 = t173 * t177 + (qJD(1) * t129 + 0.2e1 * t112 * t127) * t101;
	t1 = [t107 * t154 * t177 + (qJD(4) * t129 - t148 * t155) * t101, 0, 0, t71, 0, 0; (-t144 * t167 + (0.2e1 * t140 + (qJD(1) * t74 + t70) * t168) * t111) * t112 + (t74 * t131 * t111 + (-t74 * t137 + (t74 * t142 + (t98 * t144 + t75 * t162 + 0.2e1 * t126 * t164 + ((-t75 * t136 - t144) * t98 + (t104 * t134 - (t75 - 0.2e1 * t145) * t111) * t99) * t101) * t159) * t111 + (-t79 + (-t138 + (t138 - (t106 - t110) * t99 * t157) * t101) * t80) * t111 * qJD(1)) * t76) * t114, 0, 0, (-t148 * t167 + (-0.2e1 * t140 + (-qJD(4) * t72 - t70) * t168) * t114) * t113 + (t72 * t114 * t131 + (-t79 * t143 - ((-t112 * t71 - t147 * t84) * t99 + (t176 * t75 + t145 - t158) * t98) * t80 * t154 + (t114 * t142 + t80 * t148) * t72 - ((t71 - t147) * t98 + (t75 * t84 + qJD(4) + (-t75 - t158) * t112) * t99) * t113 * t159) * t76) * t111, 0, 0; (t95 * t165 - t91 * t94) * t141 + ((-t96 * qJD(1) + t112 * t133) * t91 + 0.2e1 * t95 * t139 + (-t94 * t89 - (-t97 * qJD(1) + t112 * t132) * t96 - t95 * t88) * t92) * t85, 0, 0, t113 * t143 * t174 + (-t148 * t174 + ((t91 * t141 + t175 * t85) * t117 + (-0.2e1 * t165 * t170 + (t88 * t92 - 0.2e1 * t139) * t85) * t118) * t114) * t111, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:18
	% EndTime: 2019-10-09 23:34:19
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (3592->98), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->94)
	t150 = pkin(10) + qJ(4);
	t144 = sin(t150);
	t137 = t144 ^ 2;
	t147 = cos(t150);
	t140 = 0.1e1 / t147 ^ 2;
	t196 = t137 * t140;
	t151 = qJ(1) + pkin(9);
	t145 = sin(t151);
	t213 = 0.2e1 * t145;
	t212 = t144 * t196;
	t149 = pkin(11) + qJ(6);
	t143 = sin(t149);
	t146 = cos(t149);
	t148 = cos(t151);
	t187 = t147 * t148;
	t126 = t143 * t145 + t146 * t187;
	t190 = t145 * t144;
	t129 = atan2(-t190, -t147);
	t128 = cos(t129);
	t127 = sin(t129);
	t176 = t127 * t190;
	t113 = -t128 * t147 - t176;
	t110 = 0.1e1 / t113;
	t120 = 0.1e1 / t126;
	t139 = 0.1e1 / t147;
	t111 = 0.1e1 / t113 ^ 2;
	t121 = 0.1e1 / t126 ^ 2;
	t211 = -0.2e1 * t144;
	t138 = t145 ^ 2;
	t132 = t138 * t196 + 0.1e1;
	t130 = 0.1e1 / t132;
	t210 = t130 - 0.1e1;
	t185 = qJD(1) * t148;
	t173 = t144 * t185;
	t183 = qJD(4) * t147;
	t184 = qJD(4) * t145;
	t104 = (-(-t145 * t183 - t173) * t139 + t184 * t196) * t130;
	t198 = t128 * t144;
	t99 = (-t104 * t145 + qJD(4)) * t198 + (-t173 + (t104 - t184) * t147) * t127;
	t209 = t110 * t111 * t99;
	t188 = t145 * t147;
	t160 = t143 * t188 + t146 * t148;
	t182 = qJD(4) * t148;
	t172 = t144 * t182;
	t105 = t160 * qJD(1) - qJD(6) * t126 + t143 * t172;
	t189 = t145 * t146;
	t125 = t143 * t187 - t189;
	t119 = t125 ^ 2;
	t118 = t119 * t121 + 0.1e1;
	t200 = t121 * t125;
	t166 = -qJD(1) * t147 + qJD(6);
	t167 = qJD(6) * t147 - qJD(1);
	t192 = t143 * t148;
	t106 = -t167 * t192 + (t166 * t145 - t172) * t146;
	t205 = t106 * t120 * t121;
	t208 = (-t105 * t200 - t119 * t205) / t118 ^ 2;
	t207 = t104 * t127;
	t206 = t104 * t144;
	t204 = t111 * t144;
	t194 = t139 * t144;
	t159 = qJD(4) * (t139 * t212 + t194);
	t164 = t137 * t145 * t185;
	t203 = (t138 * t159 + t140 * t164) / t132 ^ 2;
	t171 = 0.1e1 + t196;
	t117 = t171 * t145 * t130;
	t202 = t117 * t145;
	t201 = t120 * t143;
	t199 = t125 * t146;
	t197 = t137 * t139;
	t142 = t148 ^ 2;
	t195 = t137 * t142;
	t191 = t144 * t148;
	t186 = qJD(1) * t145;
	t109 = t111 * t195 + 0.1e1;
	t181 = 0.2e1 / t109 ^ 2 * (-t195 * t209 + (t142 * t144 * t183 - t164) * t111);
	t180 = 0.2e1 * t209;
	t179 = -0.2e1 * t208;
	t178 = t125 * t205;
	t177 = t111 * t191;
	t175 = t130 * t197;
	t170 = t144 * t181;
	t169 = t203 * t211;
	t168 = t203 * t213;
	t165 = t145 * t175;
	t163 = t171 * t148;
	t162 = t166 * t148;
	t161 = t121 * t199 - t201;
	t124 = -t146 * t188 + t192;
	t115 = 0.1e1 / t118;
	t107 = 0.1e1 / t109;
	t103 = (t210 * t144 * t127 - t128 * t165) * t148;
	t102 = -t127 * t188 + t198 + (t127 * t147 - t128 * t190) * t117;
	t100 = -t171 * t168 + (qJD(1) * t163 + t159 * t213) * t130;
	t1 = [t139 * t148 * t169 + (qJD(4) * t163 - t186 * t194) * t130, 0, 0, t100, 0, 0; (t110 * t170 + (-t110 * t183 + (qJD(1) * t103 + t99) * t204) * t107) * t145 + (t111 * t170 * t103 + (-((t104 * t165 + t210 * t183 + t169) * t127 + (t168 * t197 - t206 + (t206 + (t211 - t212) * t184) * t130) * t128) * t177 + (-t111 * t183 + t144 * t180) * t103 + (-t110 + ((-t138 + t142) * t128 * t175 + t210 * t176) * t111) * t144 * qJD(1)) * t107) * t148, 0, 0, (t102 * t204 - t110 * t147) * t148 * t181 + ((-t110 * t186 + (-qJD(4) * t102 - t99) * t148 * t111) * t147 + (-t110 * t182 - (-t100 * t128 * t145 + t127 * t184 + t202 * t207 - t207 + (-qJD(4) * t127 - t128 * t185) * t117) * t177 + (t111 * t186 + t148 * t180) * t102 - ((t100 - t185) * t127 + ((0.1e1 - t202) * qJD(4) + (t117 - t145) * t104) * t128) * t111 * t187) * t144) * t107, 0, 0; 0.2e1 * (t120 * t160 + t124 * t200) * t208 + (0.2e1 * t124 * t178 - t167 * t120 * t189 + (t144 * t184 + t162) * t201 + (t124 * t105 + t160 * t106 - t162 * t199 - (qJD(4) * t144 * t146 + t167 * t143) * t125 * t145) * t121) * t115, 0, 0, t161 * t179 * t191 + (t161 * t147 * t182 + (-t161 * t186 + ((-qJD(6) * t120 - 0.2e1 * t178) * t146 + (-t105 * t146 + (-qJD(6) * t125 + t106) * t143) * t121) * t148) * t144) * t115, 0, t179 + 0.2e1 * (-t105 * t115 * t121 + (-t115 * t205 - t121 * t208) * t125) * t125;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end