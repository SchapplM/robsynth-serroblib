% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP2
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
%   Wie in S6RPRRPP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:11
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:27
	% DurationCPUTime: 1.10s
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
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:28
	% DurationCPUTime: 1.46s
	% Computational Cost: add. (4034->125), mult. (6168->273), div. (1114->15), fcn. (7752->9), ass. (0->111)
	t157 = cos(qJ(4));
	t155 = sin(qJ(4));
	t206 = qJ(1) + pkin(9);
	t185 = sin(t206);
	t180 = t185 * t155;
	t147 = cos(t206);
	t158 = cos(qJ(3));
	t216 = t147 * t158;
	t135 = t157 * t216 + t180;
	t129 = 0.1e1 / t135 ^ 2;
	t146 = t147 ^ 2;
	t156 = sin(qJ(3));
	t151 = t156 ^ 2;
	t218 = t146 * t151;
	t195 = t129 * t218;
	t119 = 0.1e1 + t195;
	t177 = qJD(1) * t185;
	t173 = t158 * t177;
	t176 = t185 * qJD(4);
	t210 = qJD(3) * t156;
	t114 = (-t173 + t176) * t157 + (-t157 * t210 + (-qJD(4) * t158 + qJD(1)) * t155) * t147;
	t128 = 0.1e1 / t135;
	t225 = t114 * t128 * t129;
	t182 = t218 * t225;
	t209 = qJD(3) * t158;
	t231 = (-t182 + (t146 * t156 * t209 - t147 * t151 * t177) * t129) / t119 ^ 2;
	t217 = t147 * t156;
	t131 = t147 * t157 + t158 * t180;
	t172 = t155 * t176;
	t207 = qJD(4) * t157;
	t189 = t147 * t207;
	t193 = t147 * t210;
	t113 = t131 * qJD(1) + t155 * t193 - t158 * t189 - t172;
	t179 = t185 * t157;
	t134 = t155 * t216 - t179;
	t148 = 0.1e1 / t155;
	t149 = 0.1e1 / t155 ^ 2;
	t152 = 0.1e1 / t156;
	t153 = 0.1e1 / t156 ^ 2;
	t192 = t153 * t209;
	t215 = t148 * t152;
	t230 = (t149 * t152 * t207 + t148 * t192) * t134 + t113 * t215;
	t212 = t156 * t155;
	t124 = atan2(-t131, t212);
	t121 = cos(t124);
	t120 = sin(t124);
	t224 = t120 * t131;
	t112 = t121 * t212 - t224;
	t109 = 0.1e1 / t112;
	t110 = 0.1e1 / t112 ^ 2;
	t229 = 0.2e1 * t134;
	t126 = t131 ^ 2;
	t214 = t149 * t153;
	t125 = t126 * t214 + 0.1e1;
	t122 = 0.1e1 / t125;
	t170 = t155 * t209 + t156 * t207;
	t196 = t131 * t214;
	t181 = t156 * t185;
	t174 = qJD(3) * t181;
	t175 = t157 * t177;
	t208 = qJD(4) * t155;
	t211 = qJD(1) * t147;
	t115 = -t155 * t174 - t147 * t208 - t175 + (t155 * t211 + t157 * t176) * t158;
	t198 = t115 * t215;
	t101 = (t170 * t196 - t198) * t122;
	t166 = -t101 * t131 + t170;
	t97 = (-t101 * t212 - t115) * t120 + t166 * t121;
	t228 = t109 * t110 * t97;
	t150 = t148 * t149;
	t154 = t152 / t151;
	t190 = t153 * t207;
	t227 = (t115 * t196 + (-t149 * t154 * t209 - t150 * t190) * t126) / t125 ^ 2;
	t226 = t110 * t134;
	t223 = t120 * t134;
	t222 = t120 * t156;
	t221 = t121 * t131;
	t220 = t121 * t134;
	t219 = t121 * t158;
	t213 = t149 * t157;
	t127 = t134 ^ 2;
	t107 = t127 * t110 + 0.1e1;
	t205 = 0.2e1 / t107 ^ 2 * (-t113 * t226 - t127 * t228);
	t204 = 0.2e1 * t228;
	t203 = 0.2e1 * t231;
	t202 = -0.2e1 * t227;
	t201 = t152 * t227;
	t200 = t110 * t223;
	t197 = t131 * t215;
	t194 = t148 * t153 * t158;
	t191 = t157 * t209;
	t188 = t109 * t205;
	t187 = t110 * t205;
	t186 = t134 * t204;
	t184 = t217 * t229;
	t183 = t148 * t201;
	t169 = t131 * t194 + t185;
	t108 = t169 * t122;
	t178 = t185 - t108;
	t133 = -t147 * t155 + t158 * t179;
	t171 = t131 * t213 - t133 * t148;
	t168 = t129 * t133 * t147 - t185 * t128;
	t117 = 0.1e1 / t119;
	t116 = t135 * qJD(1) - t157 * t174 - t158 * t172 - t189;
	t105 = 0.1e1 / t107;
	t104 = t171 * t152 * t122;
	t100 = (-t120 + (t121 * t197 + t120) * t122) * t134;
	t99 = -t108 * t221 + (t178 * t222 + t219) * t155;
	t98 = t121 * t156 * t157 - t120 * t133 + (-t120 * t212 - t221) * t104;
	t96 = t169 * t202 + (t115 * t194 + t211 + (-t149 * t158 * t190 + (-0.2e1 * t154 * t158 ^ 2 - t152) * t148 * qJD(3)) * t131) * t122;
	t94 = -0.2e1 * t171 * t201 + (-t171 * t192 + (t115 * t213 - t116 * t148 + (t133 * t213 + (-0.2e1 * t150 * t157 ^ 2 - t148) * t131) * qJD(4)) * t152) * t122;
	t1 = [t230 * t122 + t183 * t229, 0, t96, t94, 0, 0; t131 * t188 + (-t115 * t109 + (t100 * t113 + t131 * t97) * t110) * t105 + (t100 * t187 + (t100 * t204 + (t113 * t122 - t113 - (-t101 * t122 * t197 + t202) * t134) * t110 * t120 + (-(-0.2e1 * t131 * t183 - t101) * t226 + (-(t101 + t198) * t134 + t230 * t131) * t110 * t122) * t121) * t105) * t134, 0, t99 * t134 * t187 + (t99 * t186 + (-(-t96 * t221 + (t101 * t224 - t115 * t121) * t108) * t134 + t99 * t113) * t110 + (-t109 * t217 - (-t108 * t222 + t120 * t181 + t219) * t226) * t207) * t105 + (t188 * t217 + ((-t147 * qJD(3) * t109 - (t178 * qJD(3) - t101) * t200) * t158 + (t109 * t177 + (t147 * t97 - (-t96 + t211) * t223 - (t178 * t101 - qJD(3)) * t220) * t110) * t156) * t105) * t155, (-t109 * t135 + t98 * t226) * t205 + (t98 * t186 + t114 * t109 - (-t116 + (-t101 * t157 - t155 * t94) * t156 - t166 * t104) * t200 + (t98 * t113 - t135 * t97 - (t191 - t156 * t208 - t104 * t115 - t131 * t94 + (-t104 * t212 - t133) * t101) * t220) * t110) * t105, 0, 0; t168 * t156 * t203 + (-t168 * t209 + ((qJD(1) * t128 + 0.2e1 * t133 * t225) * t147 + (-t185 * t114 - t116 * t147 + t133 * t177) * t129) * t156) * t117, 0, (t128 * t216 + t157 * t195) * t203 + (0.2e1 * t157 * t182 + (t173 + t193) * t128 + ((t114 * t158 + 0.2e1 * t151 * t175) * t147 + (t151 * t208 - 0.2e1 * t156 * t191) * t146) * t129) * t117, t129 * t184 * t231 + (t184 * t225 + (t113 * t217 + (-t147 * t209 + t156 * t177) * t134) * t129) * t117, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:11:26
	% EndTime: 2019-10-10 01:11:27
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (1618->95), mult. (2545->213), div. (484->12), fcn. (2996->9), ass. (0->99)
	t127 = qJ(1) + pkin(9);
	t125 = sin(t127);
	t123 = t125 ^ 2;
	t134 = sin(qJ(3));
	t129 = t134 ^ 2;
	t136 = cos(qJ(3));
	t131 = 0.1e1 / t136 ^ 2;
	t180 = t129 * t131;
	t120 = t123 * t180 + 0.1e1;
	t128 = t134 * t129;
	t130 = 0.1e1 / t136;
	t144 = qJD(3) * (t128 * t131 + t134) * t130;
	t126 = cos(t127);
	t176 = qJD(1) * t126;
	t185 = t125 * t129;
	t148 = t176 * t185;
	t193 = 0.1e1 / t120 ^ 2 * (t123 * t144 + t131 * t148);
	t202 = -0.2e1 * t193;
	t133 = sin(qJ(4));
	t179 = t133 * t136;
	t135 = cos(qJ(4));
	t183 = t125 * t135;
	t113 = -t126 * t179 + t183;
	t107 = t113 ^ 2;
	t178 = t135 * t136;
	t114 = t125 * t133 + t126 * t178;
	t109 = 0.1e1 / t114 ^ 2;
	t191 = t107 * t109;
	t102 = 0.1e1 + t191;
	t189 = t109 * t113;
	t150 = qJD(1) * t136 - qJD(4);
	t146 = t133 * t150;
	t171 = qJD(4) * t136;
	t151 = -qJD(1) + t171;
	t173 = qJD(3) * t134;
	t198 = -t133 * t173 + t151 * t135;
	t93 = t125 * t146 - t198 * t126;
	t165 = t93 * t189;
	t108 = 0.1e1 / t114;
	t110 = t108 * t109;
	t190 = t107 * t110;
	t156 = t135 * t173;
	t161 = t125 * t178;
	t94 = qJD(1) * t161 - qJD(4) * t183 - t133 * t176 + (t133 * t171 + t156) * t126;
	t201 = 0.1e1 / t102 ^ 2 * (t94 * t190 + t165);
	t145 = t108 * t133 + t135 * t189;
	t99 = 0.1e1 / t102;
	t200 = t145 * t99;
	t155 = 0.1e1 + t180;
	t199 = t125 * t155;
	t184 = t125 * t134;
	t119 = atan2(t184, t136);
	t116 = cos(t119);
	t115 = sin(t119);
	t163 = t115 * t184;
	t106 = t116 * t136 + t163;
	t103 = 0.1e1 / t106;
	t104 = 0.1e1 / t106 ^ 2;
	t117 = 0.1e1 / t120;
	t197 = t117 - 0.1e1;
	t124 = t126 ^ 2;
	t172 = qJD(3) * t136;
	t186 = t124 * t129;
	t174 = qJD(3) * t125;
	t158 = t131 * t174;
	t175 = qJD(1) * t134;
	t159 = t126 * t175;
	t92 = ((t125 * t172 + t159) * t130 + t129 * t158) * t117;
	t152 = -t92 + t174;
	t153 = t125 * t92 - qJD(3);
	t187 = t116 * t134;
	t87 = t153 * t187 + (t152 * t136 + t159) * t115;
	t194 = t103 * t104 * t87;
	t97 = t104 * t186 + 0.1e1;
	t196 = (-t186 * t194 + (t124 * t134 * t172 - t148) * t104) / t97 ^ 2;
	t95 = 0.1e1 / t97;
	t195 = t104 * t95;
	t188 = t115 * t136;
	t182 = t126 * t134;
	t181 = t126 * t135;
	t177 = qJD(1) * t125;
	t170 = -0.2e1 * t196;
	t169 = 0.2e1 * t201;
	t168 = -0.2e1 * t194;
	t167 = t103 * t196;
	t166 = t110 * t113 * t94;
	t164 = t95 * t172;
	t162 = t117 * t129 * t130;
	t160 = t125 * t175;
	t154 = t130 * t202;
	t149 = t125 * t162;
	t147 = t155 * t126;
	t112 = t126 * t133 - t161;
	t111 = t125 * t179 + t181;
	t101 = t117 * t199;
	t91 = (-t197 * t134 * t115 + t116 * t149) * t126;
	t90 = t125 * t188 - t187 + (t116 * t184 - t188) * t101;
	t88 = t199 * t202 + (qJD(1) * t147 + 0.2e1 * t125 * t144) * t117;
	t1 = [t154 * t182 + (qJD(3) * t147 - t130 * t160) * t117, 0, t88, 0, 0, 0; (t103 * t164 + (-0.2e1 * t167 + (-qJD(1) * t91 - t87) * t195) * t134) * t125 + (t91 * t134 * t95 * t168 + (t91 * t164 + (t91 * t170 + ((0.2e1 * t134 * t193 - t92 * t149 - t197 * t172) * t115 + (t154 * t185 + t134 * t92 + (t128 * t158 + (-t92 + 0.2e1 * t174) * t134) * t117) * t116) * t95 * t126) * t134) * t104 + (t103 + ((-t123 + t124) * t116 * t162 + t197 * t163) * t104) * t95 * t175) * t126, 0, (t103 * t95 * t177 + (0.2e1 * t167 + (qJD(3) * t90 + t87) * t195) * t126) * t136 + (((qJD(3) * t103 + t90 * t168) * t126 + (-t90 * t177 + ((t101 * t176 + t125 * t88) * t116 + ((-t101 * t125 + 0.1e1) * t92 + (t101 - t125) * qJD(3)) * t115) * t182) * t104) * t95 + (t90 * t170 + ((-t88 + t176) * t115 + (t152 * t101 + t153) * t116) * t95 * t136) * t104 * t126) * t134, 0, 0, 0; (-t108 * t111 + t112 * t189) * t169 + (-0.2e1 * t112 * t166 + (t111 * t94 - t112 * t93 + (t150 * t181 - (t151 * t133 + t156) * t125) * t113) * t109 + (t198 * t125 + t126 * t146) * t108) * t99, 0, -t160 * t200 + (t172 * t200 + (-0.2e1 * t145 * t201 + ((qJD(4) * t108 + 0.2e1 * t166) * t135 + (t135 * t93 + (-qJD(4) * t113 + t94) * t133) * t109) * t99) * t134) * t126, (t108 * t114 + t191) * t169 + (-0.2e1 * t165 + (-t109 * t114 + t108 - 0.2e1 * t190) * t94) * t99, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end