% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4RRRP7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RRRP7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:21:24
	% EndTime: 2019-12-31 17:21:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:21:24
	% EndTime: 2019-12-31 17:21:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:21:24
	% EndTime: 2019-12-31 17:21:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:21:24
	% EndTime: 2019-12-31 17:21:25
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (1002->94), mult. (2519->211), div. (480->12), fcn. (2968->9), ass. (0->92)
	t126 = sin(qJ(1));
	t119 = t126 ^ 2;
	t125 = sin(qJ(2));
	t118 = t125 ^ 2;
	t128 = cos(qJ(2));
	t121 = 0.1e1 / t128 ^ 2;
	t173 = t118 * t121;
	t113 = t119 * t173 + 0.1e1;
	t117 = t125 * t118;
	t120 = 0.1e1 / t128;
	t172 = t120 * t125;
	t137 = qJD(2) * (t117 * t120 * t121 + t172);
	t129 = cos(qJ(1));
	t163 = qJD(1) * t129;
	t151 = t126 * t163;
	t181 = 0.1e1 / t113 ^ 2 * (t119 * t137 + t151 * t173);
	t193 = -0.2e1 * t181;
	t111 = 0.1e1 / t113;
	t146 = 0.1e1 + t173;
	t190 = t126 * t146;
	t98 = t111 * t190;
	t192 = t126 * t98 - 0.1e1;
	t124 = sin(qJ(3));
	t127 = cos(qJ(3));
	t165 = t129 * t127;
	t107 = t126 * t124 + t128 * t165;
	t102 = 0.1e1 / t107 ^ 2;
	t166 = t129 * t124;
	t168 = t126 * t127;
	t106 = t128 * t166 - t168;
	t175 = t106 * t127;
	t101 = 0.1e1 / t107;
	t177 = t101 * t124;
	t139 = t102 * t175 - t177;
	t100 = t106 ^ 2;
	t99 = t100 * t102 + 0.1e1;
	t96 = 0.1e1 / t99;
	t191 = t139 * t96;
	t169 = t126 * t125;
	t110 = atan2(-t169, -t128);
	t109 = cos(t110);
	t108 = sin(t110);
	t154 = t108 * t169;
	t94 = -t109 * t128 - t154;
	t91 = 0.1e1 / t94;
	t92 = 0.1e1 / t94 ^ 2;
	t189 = t111 - 0.1e1;
	t123 = t129 ^ 2;
	t161 = qJD(2) * t128;
	t155 = t92 * t161;
	t150 = t125 * t163;
	t162 = qJD(2) * t126;
	t174 = t109 * t125;
	t149 = t121 * t162;
	t85 = (-(-t126 * t161 - t150) * t120 + t118 * t149) * t111;
	t80 = (-t126 * t85 + qJD(2)) * t174 + (-t150 + (t85 - t162) * t128) * t108;
	t187 = t80 * t91 * t92;
	t90 = t123 * t118 * t92 + 0.1e1;
	t188 = (t123 * t125 * t155 + (-t123 * t187 - t92 * t151) * t118) / t90 ^ 2;
	t176 = t102 * t106;
	t143 = -qJD(1) * t128 + qJD(3);
	t144 = qJD(3) * t128 - qJD(1);
	t160 = qJD(2) * t129;
	t148 = t125 * t160;
	t87 = -t144 * t166 + (t143 * t126 - t148) * t127;
	t183 = t101 * t102 * t87;
	t167 = t126 * t128;
	t138 = t124 * t167 + t165;
	t86 = t138 * qJD(1) - t107 * qJD(3) + t124 * t148;
	t186 = (-t100 * t183 - t86 * t176) / t99 ^ 2;
	t88 = 0.1e1 / t90;
	t185 = t88 * t91;
	t184 = t88 * t92;
	t179 = t129 * t92;
	t178 = qJD(2) * t98;
	t171 = t125 * t129;
	t164 = qJD(1) * t126;
	t159 = 0.2e1 * t187;
	t158 = -0.2e1 * t186;
	t157 = t91 * t188;
	t156 = t106 * t183;
	t153 = t111 * t118 * t120;
	t147 = 0.2e1 * t92 * t188;
	t145 = t120 * t193;
	t142 = t126 * t153;
	t141 = t146 * t129;
	t140 = t143 * t129;
	t105 = -t127 * t167 + t166;
	t84 = (t189 * t125 * t108 - t109 * t142) * t129;
	t83 = -t192 * t174 + (-t126 + t98) * t128 * t108;
	t81 = t190 * t193 + (qJD(1) * t141 + 0.2e1 * t126 * t137) * t111;
	t1 = [t145 * t171 + (qJD(2) * t141 - t164 * t172) * t111, t81, 0, 0; (-t161 * t185 + (0.2e1 * t157 + (qJD(1) * t84 + t80) * t184) * t125) * t126 + (t84 * t147 * t125 + (-t84 * t155 + (t84 * t159 + ((0.2e1 * t125 * t181 - t85 * t142 - t189 * t161) * t108 + (t118 * t126 * t145 + t125 * t85 + (t117 * t149 - (t85 - 0.2e1 * t162) * t125) * t111) * t109) * t179) * t125 + (-t91 + (-(t119 - t123) * t109 * t153 + t189 * t154) * t92) * t125 * qJD(1)) * t88) * t129, (-t164 * t185 + (-0.2e1 * t157 + (-qJD(2) * t83 - t80) * t184) * t129) * t128 + (t83 * t129 * t147 + (-t91 * t160 - ((-t126 * t81 - t163 * t98) * t109 + (t192 * t85 + t162 - t178) * t108) * t92 * t171 + (t129 * t159 + t92 * t164) * t83 - ((t81 - t163) * t108 + (t85 * t98 + qJD(2) + (-t85 - t178) * t126) * t109) * t128 * t179) * t88) * t125, 0, 0; 0.2e1 * (t101 * t138 + t105 * t176) * t186 + (0.2e1 * t105 * t156 - t144 * t101 * t168 + (t125 * t162 + t140) * t177 + (t138 * t87 + t105 * t86 - t140 * t175 - (qJD(2) * t125 * t127 + t144 * t124) * t106 * t126) * t102) * t96, t128 * t160 * t191 + (-t164 * t191 + (t139 * t158 + ((-qJD(3) * t101 - 0.2e1 * t156) * t127 + (-t127 * t86 + (-qJD(3) * t106 + t87) * t124) * t102) * t96) * t129) * t125, t158 + 0.2e1 * (-t86 * t102 * t96 + (-t102 * t186 - t96 * t183) * t106) * t106, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:21:24
	% EndTime: 2019-12-31 17:21:25
	% DurationCPUTime: 0.90s
	% Computational Cost: add. (1824->123), mult. (6168->270), div. (1114->15), fcn. (7752->9), ass. (0->112)
	t155 = sin(qJ(3));
	t158 = cos(qJ(2));
	t157 = cos(qJ(3));
	t159 = cos(qJ(1));
	t211 = t159 * t157;
	t230 = sin(qJ(1));
	t136 = t230 * t155 + t158 * t211;
	t130 = 0.1e1 / t136 ^ 2;
	t156 = sin(qJ(2));
	t150 = t156 ^ 2;
	t154 = t159 ^ 2;
	t215 = t150 * t154;
	t193 = t130 * t215;
	t125 = 0.1e1 + t193;
	t184 = qJD(1) * t230;
	t208 = qJD(2) * t159;
	t188 = t156 * t208;
	t169 = t158 * t184 + t188;
	t183 = t230 * qJD(3);
	t212 = t159 * t155;
	t115 = (-qJD(3) * t158 + qJD(1)) * t212 + (t183 - t169) * t157;
	t129 = 0.1e1 / t136;
	t225 = t115 * t129 * t130;
	t177 = t215 * t225;
	t189 = qJD(2) * t154 * t156;
	t233 = (-t177 + (-t150 * t159 * t184 + t158 * t189) * t130) / t125 ^ 2;
	t213 = t156 * t159;
	t191 = t230 * t158;
	t132 = t155 * t191 + t211;
	t174 = t155 * t183;
	t205 = qJD(3) * t159;
	t186 = t157 * t205;
	t114 = t132 * qJD(1) + t155 * t188 - t158 * t186 - t174;
	t135 = -t230 * t157 + t158 * t212;
	t147 = 0.1e1 / t155;
	t148 = 0.1e1 / t155 ^ 2;
	t151 = 0.1e1 / t156;
	t152 = 0.1e1 / t156 ^ 2;
	t209 = qJD(2) * t158;
	t190 = t152 * t209;
	t206 = qJD(3) * t157;
	t218 = t147 * t151;
	t232 = (t148 * t151 * t206 + t147 * t190) * t135 + t114 * t218;
	t214 = t156 * t155;
	t124 = atan2(-t132, t214);
	t119 = cos(t124);
	t118 = sin(t124);
	t224 = t118 * t132;
	t113 = t119 * t214 - t224;
	t110 = 0.1e1 / t113;
	t111 = 0.1e1 / t113 ^ 2;
	t231 = 0.2e1 * t135;
	t127 = t132 ^ 2;
	t217 = t148 * t152;
	t126 = t127 * t217 + 0.1e1;
	t122 = 0.1e1 / t126;
	t170 = t155 * t209 + t156 * t206;
	t195 = t132 * t217;
	t192 = t230 * t156;
	t175 = qJD(2) * t192;
	t176 = t157 * t184;
	t210 = qJD(1) * t159;
	t116 = t157 * t183 * t158 - t176 + (t210 * t158 - t175 - t205) * t155;
	t197 = t116 * t218;
	t102 = (t170 * t195 - t197) * t122;
	t167 = -t102 * t132 + t170;
	t98 = (-t102 * t214 - t116) * t118 + t167 * t119;
	t229 = t110 * t111 * t98;
	t149 = t147 * t148;
	t153 = t151 / t150;
	t187 = t152 * t206;
	t228 = (t116 * t195 + (-t148 * t153 * t209 - t149 * t187) * t127) / t126 ^ 2;
	t227 = t111 * t135;
	t226 = t114 * t111;
	t223 = t118 * t135;
	t222 = t118 * t156;
	t221 = t119 * t132;
	t220 = t119 * t135;
	t219 = t119 * t158;
	t216 = t148 * t157;
	t207 = qJD(3) * t155;
	t128 = t135 ^ 2;
	t108 = t128 * t111 + 0.1e1;
	t204 = 0.2e1 / t108 ^ 2 * (-t128 * t229 - t135 * t226);
	t203 = 0.2e1 * t229;
	t202 = 0.2e1 * t233;
	t201 = -0.2e1 * t228;
	t200 = t151 * t228;
	t199 = t111 * t223;
	t196 = t132 * t218;
	t194 = t147 * t152 * t158;
	t172 = t132 * t194 + t230;
	t109 = t172 * t122;
	t185 = t230 - t109;
	t182 = t110 * t204;
	t181 = t111 * t204;
	t180 = t135 * t203;
	t179 = t213 * t231;
	t178 = t147 * t200;
	t134 = t157 * t191 - t212;
	t173 = t132 * t216 - t134 * t147;
	t171 = t130 * t134 * t159 - t230 * t129;
	t120 = 0.1e1 / t125;
	t117 = t136 * qJD(1) - t157 * t175 - t158 * t174 - t186;
	t106 = 0.1e1 / t108;
	t105 = t173 * t151 * t122;
	t101 = (-t118 + (t119 * t196 + t118) * t122) * t135;
	t100 = -t109 * t221 + (t185 * t222 + t219) * t155;
	t99 = t119 * t156 * t157 - t118 * t134 + (-t118 * t214 - t221) * t105;
	t97 = t172 * t201 + (t116 * t194 + t210 + (-t148 * t158 * t187 + (-0.2e1 * t153 * t158 ^ 2 - t151) * t147 * qJD(2)) * t132) * t122;
	t95 = -0.2e1 * t173 * t200 + (-t173 * t190 + (t116 * t216 - t117 * t147 + (t134 * t216 + (-0.2e1 * t149 * t157 ^ 2 - t147) * t132) * qJD(3)) * t151) * t122;
	t1 = [t232 * t122 + t178 * t231, t97, t95, 0; t132 * t182 + (-t116 * t110 + (t101 * t114 + t132 * t98) * t111) * t106 + (t101 * t181 + (t101 * t203 + (t114 * t122 - t114 - (-t102 * t122 * t196 + t201) * t135) * t111 * t118 + (-(-0.2e1 * t132 * t178 - t102) * t227 + (-(t102 + t197) * t135 + t232 * t132) * t111 * t122) * t119) * t106) * t135, t100 * t135 * t181 + (-(-t97 * t221 + (t102 * t224 - t116 * t119) * t109) * t227 + (t180 + t226) * t100 + (-t110 * t213 - (-t109 * t222 + t118 * t192 + t219) * t227) * t206) * t106 + (t182 * t213 + ((-t110 * t208 - (t185 * qJD(2) - t102) * t199) * t158 + (t110 * t184 + (t159 * t98 - (-t97 + t210) * t223 - (t185 * t102 - qJD(2)) * t220) * t111) * t156) * t106) * t155, (-t110 * t136 + t99 * t227) * t204 + (t99 * t180 + t115 * t110 - (-t117 + (-t102 * t157 - t155 * t95) * t156 - t167 * t105) * t199 + (t99 * t114 - t136 * t98 - (t157 * t209 - t156 * t207 - t105 * t116 - t132 * t95 + (-t105 * t214 - t134) * t102) * t220) * t111) * t106, 0; t171 * t156 * t202 + (-t171 * t209 + ((qJD(1) * t129 + 0.2e1 * t134 * t225) * t159 + (-t230 * t115 - t117 * t159 + t134 * t184) * t130) * t156) * t120, (t129 * t158 * t159 + t157 * t193) * t202 + (0.2e1 * t157 * t177 + t169 * t129 + ((t115 * t159 - 0.2e1 * t157 * t189) * t158 + (t154 * t207 + 0.2e1 * t159 * t176) * t150) * t130) * t120, t130 * t179 * t233 + (t179 * t225 + (t114 * t213 + (t156 * t184 - t158 * t208) * t135) * t130) * t120, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end