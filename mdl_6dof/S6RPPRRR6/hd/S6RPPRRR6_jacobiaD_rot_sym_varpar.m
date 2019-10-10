% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR6
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
%   Wie in S6RPPRRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:03
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (624->91), mult. (2519->209), div. (480->12), fcn. (2968->9), ass. (0->95)
	t131 = sin(qJ(1));
	t171 = qJD(4) * t131;
	t125 = t131 ^ 2;
	t130 = sin(qJ(4));
	t123 = 0.1e1 / t130 ^ 2;
	t133 = cos(qJ(4));
	t127 = t133 ^ 2;
	t183 = t123 * t127;
	t118 = t125 * t183 + 0.1e1;
	t115 = 0.1e1 / t118;
	t122 = 0.1e1 / t130;
	t158 = t123 * t171;
	t134 = cos(qJ(1));
	t173 = qJD(1) * t134;
	t159 = t133 * t173;
	t90 = ((-t130 * t171 + t159) * t122 - t127 * t158) * t115;
	t150 = -t90 - t171;
	t154 = 0.1e1 + t183;
	t196 = t131 * t154;
	t178 = t131 * t133;
	t117 = atan2(t178, t130);
	t114 = cos(t117);
	t113 = sin(t117);
	t163 = t113 * t178;
	t99 = t114 * t130 + t163;
	t96 = 0.1e1 / t99;
	t132 = cos(qJ(5));
	t177 = t132 * t134;
	t161 = t130 * t177;
	t129 = sin(qJ(5));
	t180 = t131 * t129;
	t112 = t161 - t180;
	t106 = 0.1e1 / t112;
	t107 = 0.1e1 / t112 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t195 = t115 - 0.1e1;
	t151 = t131 * t90 + qJD(4);
	t184 = t114 * t133;
	t85 = t151 * t184 + (t150 * t130 + t159) * t113;
	t194 = t85 * t96 * t97;
	t128 = t134 ^ 2;
	t95 = t127 * t128 * t97 + 0.1e1;
	t93 = 0.1e1 / t95;
	t193 = t93 * t97;
	t192 = t96 * t93;
	t179 = t131 * t132;
	t181 = t130 * t134;
	t111 = t129 * t181 + t179;
	t105 = t111 ^ 2;
	t104 = t105 * t107 + 0.1e1;
	t187 = t107 * t111;
	t148 = qJD(1) * t130 + qJD(5);
	t149 = qJD(5) * t130 + qJD(1);
	t169 = qJD(4) * t134;
	t157 = t133 * t169;
	t182 = t129 * t134;
	t92 = -t149 * t182 + (-t148 * t131 + t157) * t132;
	t190 = t106 * t107 * t92;
	t91 = -qJD(5) * t161 - t129 * t157 - t132 * t173 + t148 * t180;
	t191 = 0.1e1 / t104 ^ 2 * (-t105 * t190 - t91 * t187);
	t126 = t133 * t127;
	t145 = (t123 * t126 + t133) * t122;
	t160 = t131 * t173;
	t189 = (-t145 * t125 * qJD(4) + t160 * t183) / t118 ^ 2;
	t188 = t106 * t129;
	t186 = t111 * t132;
	t185 = t113 * t131;
	t176 = t133 * t134;
	t175 = qJD(1) * t131;
	t174 = qJD(1) * t133;
	t172 = qJD(4) * t130;
	t170 = qJD(4) * t133;
	t164 = t97 * t172;
	t168 = -0.2e1 * (-t128 * t133 * t164 + (-t128 * t194 - t97 * t160) * t127) / t95 ^ 2;
	t167 = -0.2e1 * t194;
	t166 = 0.2e1 * t191;
	t165 = 0.2e1 * t189;
	t162 = t115 * t122 * t127;
	t156 = t96 * t168;
	t155 = t97 * t168;
	t153 = 0.2e1 * t111 * t190;
	t152 = -0.2e1 * t122 * t189;
	t147 = t131 * t162;
	t146 = t154 * t134;
	t144 = t148 * t134;
	t143 = t107 * t186 - t188;
	t142 = t143 * t134;
	t110 = -t130 * t179 - t182;
	t109 = -t130 * t180 + t177;
	t103 = t115 * t196;
	t101 = 0.1e1 / t104;
	t89 = (-t195 * t133 * t113 + t114 * t147) * t134;
	t88 = -t130 * t185 + t184 - (-t113 * t130 + t114 * t178) * t103;
	t86 = t165 * t196 + (-qJD(1) * t146 + 0.2e1 * t145 * t171) * t115;
	t1 = [t152 * t176 + (-t122 * t131 * t174 - qJD(4) * t146) * t115, 0, 0, t86, 0, 0; (-t172 * t192 + (t156 + (-qJD(1) * t89 - t85) * t193) * t133) * t131 + (t89 * t155 * t133 + (-t89 * t164 + (t89 * t167 + ((t133 * t165 - t90 * t147 + t195 * t172) * t113 + (t127 * t131 * t152 + t133 * t90 + (-t126 * t158 + (-t90 - 0.2e1 * t171) * t133) * t115) * t114) * t97 * t134) * t133 + (t96 + ((-t125 + t128) * t114 * t162 + t195 * t163) * t97) * t174) * t93) * t134, 0, 0, (-t175 * t192 + (t156 + (-qJD(4) * t88 - t85) * t193) * t134) * t130 + (t88 * t134 * t155 + (t96 * t169 + (t114 * t131 * t86 + t150 * t113 - (-qJD(4) * t113 + t114 * t173 - t185 * t90) * t103) * t97 * t176 + (t134 * t167 - t97 * t175) * t88) * t93 + ((-t86 - t173) * t113 + (-t150 * t103 - t151) * t114) * t181 * t193) * t133, 0, 0; (-t106 * t109 + t110 * t187) * t166 + (t110 * t153 - t149 * t106 * t179 + (-t131 * t170 - t144) * t188 + (-t109 * t92 + t110 * t91 + t144 * t186 - (t149 * t129 - t132 * t170) * t111 * t131) * t107) * t101, 0, 0, t133 * t142 * t166 + (t142 * t172 + (t143 * t175 + ((qJD(5) * t106 + t153) * t132 + (t132 * t91 + (qJD(5) * t111 - t92) * t129) * t107) * t134) * t133) * t101, -0.2e1 * t191 + 0.2e1 * (-t101 * t107 * t91 + (-t101 * t190 - t107 * t191) * t111) * t111, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:03
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (1192->94), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->97)
	t154 = sin(qJ(4));
	t147 = 0.1e1 / t154 ^ 2;
	t156 = cos(qJ(4));
	t151 = t156 ^ 2;
	t201 = t147 * t151;
	t219 = t156 * t201;
	t155 = sin(qJ(1));
	t177 = 0.1e1 + t201;
	t218 = t155 * t177;
	t149 = t155 ^ 2;
	t141 = t149 * t201 + 0.1e1;
	t138 = 0.1e1 / t141;
	t146 = 0.1e1 / t154;
	t157 = cos(qJ(1));
	t191 = qJD(1) * t157;
	t179 = t156 * t191;
	t189 = qJD(4) * t155;
	t113 = ((-t154 * t189 + t179) * t146 - t189 * t201) * t138;
	t217 = -t113 - t189;
	t195 = t155 * t156;
	t140 = atan2(t195, t154);
	t137 = cos(t140);
	t136 = sin(t140);
	t182 = t136 * t195;
	t123 = t137 * t154 + t182;
	t120 = 0.1e1 / t123;
	t153 = qJ(5) + qJ(6);
	t144 = cos(t153);
	t198 = t154 * t157;
	t180 = t144 * t198;
	t143 = sin(t153);
	t197 = t155 * t143;
	t133 = t180 - t197;
	t127 = 0.1e1 / t133;
	t121 = 0.1e1 / t123 ^ 2;
	t128 = 0.1e1 / t133 ^ 2;
	t216 = -0.2e1 * t156;
	t215 = t138 - 0.1e1;
	t152 = t157 ^ 2;
	t200 = t151 * t152;
	t118 = t121 * t200 + 0.1e1;
	t199 = t151 * t155;
	t170 = t191 * t199;
	t188 = qJD(4) * t156;
	t203 = t137 * t156;
	t107 = (t113 * t155 + qJD(4)) * t203 + (t217 * t154 + t179) * t136;
	t212 = t107 * t120 * t121;
	t214 = (-t200 * t212 + (-t152 * t154 * t188 - t170) * t121) / t118 ^ 2;
	t145 = qJD(5) + qJD(6);
	t172 = qJD(1) * t154 + t145;
	t187 = qJD(4) * t157;
	t178 = t156 * t187;
	t111 = -t143 * t178 - t144 * t191 - t145 * t180 + t172 * t197;
	t196 = t155 * t144;
	t132 = t143 * t198 + t196;
	t126 = t132 ^ 2;
	t119 = t126 * t128 + 0.1e1;
	t206 = t128 * t132;
	t173 = t145 * t154 + qJD(1);
	t202 = t143 * t157;
	t112 = -t173 * t202 + (-t172 * t155 + t178) * t144;
	t211 = t112 * t127 * t128;
	t213 = (-t111 * t206 - t126 * t211) / t119 ^ 2;
	t210 = t113 * t156;
	t209 = t121 * t156;
	t168 = (t156 + t219) * t146;
	t208 = (-t168 * t149 * qJD(4) + t147 * t170) / t141 ^ 2;
	t207 = t127 * t143;
	t205 = t132 * t144;
	t204 = t136 * t155;
	t194 = t156 * t157;
	t193 = qJD(1) * t155;
	t192 = qJD(1) * t156;
	t190 = qJD(4) * t154;
	t186 = 0.2e1 * t213;
	t185 = -0.2e1 * t212;
	t184 = 0.2e1 * t208;
	t183 = t121 * t194;
	t181 = t138 * t146 * t151;
	t176 = t214 * t216;
	t175 = 0.2e1 * t132 * t211;
	t174 = -0.2e1 * t146 * t208;
	t171 = t155 * t181;
	t169 = t177 * t157;
	t167 = t128 * t205 - t207;
	t166 = t167 * t157;
	t165 = -t155 * t188 - t172 * t157;
	t131 = -t154 * t196 - t202;
	t130 = t144 * t157 - t154 * t197;
	t125 = t138 * t218;
	t116 = 0.1e1 / t119;
	t114 = 0.1e1 / t118;
	t110 = (-t215 * t156 * t136 + t137 * t171) * t157;
	t109 = -t154 * t204 + t203 - (-t136 * t154 + t137 * t195) * t125;
	t108 = t184 * t218 + (-qJD(1) * t169 + 0.2e1 * t168 * t189) * t138;
	t104 = -0.2e1 * t213 + 0.2e1 * (-t111 * t116 * t128 + (-t116 * t211 - t128 * t213) * t132) * t132;
	t1 = [t174 * t194 + (-t146 * t155 * t192 - qJD(4) * t169) * t138, 0, 0, t108, 0, 0; (t120 * t176 + (-t120 * t190 + (-qJD(1) * t110 - t107) * t209) * t114) * t155 + (t121 * t176 * t110 + (((-t113 * t171 + t156 * t184 + t215 * t190) * t136 + (t174 * t199 + t210 + (-t210 + (t216 - t219) * t189) * t138) * t137) * t183 + (-t121 * t190 + t156 * t185) * t110 + (t120 + ((-t149 + t152) * t137 * t181 + t215 * t182) * t121) * t192) * t114) * t157, 0, 0, 0.2e1 * (-t109 * t209 - t120 * t154) * t157 * t214 + ((-t120 * t193 + (-qJD(4) * t109 - t107) * t157 * t121) * t154 + (t120 * t187 + (t108 * t137 * t155 + t217 * t136 - (-qJD(4) * t136 - t113 * t204 + t137 * t191) * t125) * t183 + (-t121 * t193 + t157 * t185) * t109 + ((-t108 - t191) * t136 + ((t125 * t155 - 0.1e1) * qJD(4) + (t125 - t155) * t113) * t137) * t121 * t198) * t156) * t114, 0, 0; (-t127 * t130 + t131 * t206) * t186 + (t131 * t175 - t173 * t127 * t196 + t165 * t207 + (-t173 * t132 * t197 + t131 * t111 - t130 * t112 - t165 * t205) * t128) * t116, 0, 0, t156 * t166 * t186 + (t166 * t190 + (t167 * t193 + ((t127 * t145 + t175) * t144 + (t111 * t144 + (t132 * t145 - t112) * t143) * t128) * t157) * t156) * t116, t104, t104;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end