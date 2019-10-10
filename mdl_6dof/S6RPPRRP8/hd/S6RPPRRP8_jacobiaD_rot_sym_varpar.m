% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP8
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
%   Wie in S6RPPRRP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:59
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:42
	% DurationCPUTime: 0.99s
	% Computational Cost: add. (2081->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->95)
	t136 = cos(qJ(1));
	t200 = 0.2e1 * t136;
	t130 = pkin(9) + qJ(4);
	t128 = sin(t130);
	t129 = cos(t130);
	t177 = t136 * t129;
	t118 = atan2(-t177, t128);
	t116 = sin(t118);
	t117 = cos(t118);
	t102 = -t116 * t177 + t117 * t128;
	t99 = 0.1e1 / t102;
	t134 = sin(qJ(1));
	t135 = cos(qJ(5));
	t178 = t134 * t135;
	t133 = sin(qJ(5));
	t180 = t133 * t136;
	t113 = t128 * t178 + t180;
	t109 = 0.1e1 / t113;
	t123 = 0.1e1 / t128;
	t100 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t113 ^ 2;
	t124 = 0.1e1 / t128 ^ 2;
	t132 = t136 ^ 2;
	t127 = t129 ^ 2;
	t183 = t124 * t127;
	t121 = t132 * t183 + 0.1e1;
	t119 = 0.1e1 / t121;
	t199 = t119 - 0.1e1;
	t131 = t134 ^ 2;
	t174 = qJD(1) * t136;
	t152 = t127 * t134 * t174;
	t172 = qJD(4) * t129;
	t182 = t127 * t131;
	t171 = qJD(4) * t136;
	t162 = t124 * t171;
	t175 = qJD(1) * t134;
	t163 = t129 * t175;
	t93 = ((t128 * t171 + t163) * t123 + t127 * t162) * t119;
	t157 = -t93 + t171;
	t158 = -t136 * t93 + qJD(4);
	t186 = t117 * t129;
	t88 = t158 * t186 + (t157 * t128 + t163) * t116;
	t196 = t99 * t100 * t88;
	t96 = t100 * t182 + 0.1e1;
	t198 = (-t182 * t196 + (-t128 * t131 * t172 + t152) * t100) / t96 ^ 2;
	t94 = 0.1e1 / t96;
	t197 = t100 * t94;
	t176 = t136 * t135;
	t179 = t134 * t133;
	t112 = t128 * t179 - t176;
	t108 = t112 ^ 2;
	t107 = t108 * t110 + 0.1e1;
	t189 = t110 * t112;
	t155 = qJD(1) * t128 + qJD(5);
	t148 = t155 * t136;
	t156 = qJD(5) * t128 + qJD(1);
	t150 = t156 * t133;
	t98 = t135 * t148 + (t135 * t172 - t150) * t134;
	t194 = t109 * t110 * t98;
	t149 = t156 * t135;
	t97 = t134 * t149 + (t134 * t172 + t148) * t133;
	t195 = 0.1e1 / t107 ^ 2 * (-t108 * t194 + t97 * t189);
	t126 = t129 * t127;
	t184 = t123 * t129;
	t146 = qJD(4) * (-t123 * t124 * t126 - t184);
	t191 = (-t124 * t152 + t132 * t146) / t121 ^ 2;
	t190 = t109 * t133;
	t188 = t112 * t135;
	t187 = t116 * t128;
	t185 = t123 * t127;
	t173 = qJD(4) * t128;
	t170 = -0.2e1 * t198;
	t169 = -0.2e1 * t196;
	t168 = 0.2e1 * t195;
	t167 = t99 * t198;
	t166 = t94 * t173;
	t165 = t129 * t191;
	t164 = t119 * t185;
	t161 = 0.1e1 + t183;
	t160 = 0.2e1 * t112 * t194;
	t159 = t191 * t200;
	t154 = t136 * t164;
	t153 = t199 * t129 * t116;
	t151 = t161 * t134;
	t147 = t110 * t188 - t190;
	t145 = t147 * t134;
	t144 = t129 * t171 - t155 * t134;
	t115 = t128 * t176 - t179;
	t114 = t128 * t180 + t178;
	t105 = 0.1e1 / t107;
	t104 = t161 * t136 * t119;
	t92 = (-t117 * t154 - t153) * t134;
	t90 = t136 * t187 + t186 + (-t117 * t177 - t187) * t104;
	t89 = -t161 * t159 + (-qJD(1) * t151 + t146 * t200) * t119;
	t1 = [-0.2e1 * t134 * t123 * t165 + (-qJD(4) * t151 + t174 * t184) * t119, 0, 0, t89, 0, 0; (t99 * t166 + (0.2e1 * t167 + (qJD(1) * t92 + t88) * t197) * t129) * t136 + ((-t92 * t166 + (t92 * t170 + ((t93 * t154 + t199 * t173 + 0.2e1 * t165) * t116 + (t159 * t185 + t93 * t129 + (t126 * t162 + (-t93 + 0.2e1 * t171) * t129) * t119) * t117) * t94 * t134) * t129) * t100 + (t92 * t169 + (t99 + ((t131 - t132) * t117 * t164 - t136 * t153) * t100) * qJD(1)) * t129 * t94) * t134, 0, 0, (t99 * t94 * t174 + (-0.2e1 * t167 + (-qJD(4) * t90 - t88) * t197) * t134) * t128 + (((qJD(4) * t99 + t90 * t169) * t134 + (t90 * t174 + ((t104 * t175 - t136 * t89) * t117 + ((t104 * t136 - 0.1e1) * t93 + (-t104 + t136) * qJD(4)) * t116) * t129 * t134) * t100) * t94 + (t90 * t170 + ((-t89 - t175) * t116 + (t157 * t104 - t158) * t117) * t94 * t128) * t100 * t134) * t129, 0, 0; (-t109 * t114 + t115 * t189) * t168 + (t115 * t160 + t136 * t109 * t149 + t144 * t190 + (t136 * t112 * t150 - t114 * t98 - t115 * t97 - t144 * t188) * t110) * t105, 0, 0, t129 * t145 * t168 + (t145 * t173 + (-t147 * t174 + ((qJD(5) * t109 + t160) * t135 + (-t135 * t97 + (qJD(5) * t112 - t98) * t133) * t110) * t134) * t129) * t105, -0.2e1 * t195 + 0.2e1 * (t105 * t110 * t97 + (-t105 * t194 - t110 * t195) * t112) * t112, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:42
	% DurationCPUTime: 1.38s
	% Computational Cost: add. (3877->122), mult. (6168->268), div. (1114->15), fcn. (7752->9), ass. (0->115)
	t153 = pkin(9) + qJ(4);
	t152 = cos(t153);
	t159 = sin(qJ(1));
	t214 = t152 * t159;
	t151 = sin(t153);
	t158 = sin(qJ(5));
	t233 = cos(qJ(1));
	t193 = t233 * t158;
	t160 = cos(qJ(5));
	t211 = t159 * t160;
	t139 = t151 * t193 + t211;
	t215 = t152 * t158;
	t129 = atan2(t139, t215);
	t124 = cos(t129);
	t123 = sin(t129);
	t224 = t123 * t139;
	t118 = t124 * t215 + t224;
	t115 = 0.1e1 / t118;
	t138 = t151 * t211 + t193;
	t133 = 0.1e1 / t138;
	t148 = 0.1e1 / t152;
	t154 = 0.1e1 / t158;
	t116 = 0.1e1 / t118 ^ 2;
	t134 = 0.1e1 / t138 ^ 2;
	t155 = 0.1e1 / t158 ^ 2;
	t192 = t233 * t160;
	t212 = t159 * t158;
	t137 = t151 * t212 - t192;
	t132 = t137 ^ 2;
	t113 = t116 * t132 + 0.1e1;
	t186 = qJD(1) * t233;
	t176 = t151 * t186;
	t170 = t233 * qJD(5) + t176;
	t181 = qJD(5) * t151 + qJD(1);
	t208 = qJD(4) * t159;
	t190 = t152 * t208;
	t121 = t181 * t211 + (t170 + t190) * t158;
	t227 = t121 * t116;
	t136 = t139 ^ 2;
	t149 = 0.1e1 / t152 ^ 2;
	t216 = t149 * t155;
	t131 = t136 * t216 + 0.1e1;
	t127 = 0.1e1 / t131;
	t205 = qJD(5) * t160;
	t209 = qJD(4) * t151;
	t171 = t152 * t205 - t158 * t209;
	t196 = t139 * t216;
	t194 = t233 * t152;
	t177 = qJD(4) * t194;
	t178 = t160 * t186;
	t179 = t151 * t192;
	t119 = -qJD(5) * t179 - t158 * t177 - t178 + (qJD(1) * t151 + qJD(5)) * t212;
	t218 = t148 * t154;
	t199 = t119 * t218;
	t107 = (-t171 * t196 - t199) * t127;
	t169 = -t107 * t139 - t171;
	t103 = (-t107 * t215 - t119) * t123 - t169 * t124;
	t117 = t115 * t116;
	t231 = t103 * t117;
	t232 = (-t132 * t231 + t137 * t227) / t113 ^ 2;
	t147 = t152 ^ 2;
	t157 = t159 ^ 2;
	t219 = t147 * t157;
	t198 = t134 * t219;
	t130 = 0.1e1 + t198;
	t207 = qJD(4) * t160;
	t189 = t152 * t207;
	t122 = t170 * t160 + (-t181 * t158 + t189) * t159;
	t226 = t122 * t133 * t134;
	t180 = t219 * t226;
	t230 = (-t180 + (t147 * t159 * t186 - t152 * t157 * t209) * t134) / t130 ^ 2;
	t150 = t148 / t147;
	t156 = t154 * t155;
	t229 = (-t119 * t196 + (-t149 * t156 * t205 + t150 * t155 * t209) * t136) / t131 ^ 2;
	t228 = t116 * t137;
	t225 = t123 * t137;
	t223 = t123 * t152;
	t222 = t124 * t137;
	t221 = t124 * t139;
	t220 = t124 * t151;
	t217 = t149 * t151;
	t213 = t155 * t160;
	t210 = qJD(1) * t159;
	t206 = qJD(5) * t158;
	t204 = 0.2e1 * t232;
	t203 = 0.2e1 * t230;
	t202 = -0.2e1 * t229;
	t201 = 0.2e1 * t117 * t137;
	t200 = t116 * t225;
	t197 = t139 * t218;
	t195 = t154 * t217;
	t191 = t149 * t209;
	t188 = t155 * t205;
	t173 = t139 * t195 + t233;
	t114 = t173 * t127;
	t187 = t233 - t114;
	t185 = -0.2e1 * t115 * t232;
	t184 = t116 * t204;
	t183 = 0.2e1 * t148 * t229;
	t182 = -0.2e1 * t137 * t214;
	t175 = t154 * t183;
	t140 = t179 - t212;
	t174 = t139 * t213 - t140 * t154;
	t172 = t134 * t140 * t159 - t233 * t133;
	t168 = t121 * t218 - (t148 * t188 - t154 * t191) * t137;
	t125 = 0.1e1 / t130;
	t120 = t138 * qJD(1) + qJD(5) * t139 - t160 * t177;
	t111 = 0.1e1 / t113;
	t110 = t174 * t148 * t127;
	t106 = (-t123 + (-t124 * t197 + t123) * t127) * t137;
	t105 = t114 * t221 + (t187 * t223 - t220) * t158;
	t104 = t124 * t152 * t160 + t123 * t140 - (-t123 * t215 + t221) * t110;
	t102 = t173 * t202 + (-t119 * t195 - t210 + (-t188 * t217 + (0.2e1 * t150 * t151 ^ 2 + t148) * t154 * qJD(4)) * t139) * t127;
	t100 = t174 * t183 + (-t174 * t191 + (t119 * t213 - t120 * t154 + (-t140 * t213 + (0.2e1 * t156 * t160 ^ 2 + t154) * t139) * qJD(5)) * t148) * t127;
	t1 = [-t168 * t127 + t137 * t175, 0, 0, t102, t100, 0; t139 * t185 + (-t119 * t115 + (-t103 * t139 - t106 * t121) * t116) * t111 + (t106 * t184 + (0.2e1 * t106 * t231 + (-t121 * t127 + t121 - (t107 * t127 * t197 + t202) * t137) * t116 * t123 + (-(t139 * t175 - t107) * t228 + (-(t107 + t199) * t137 + t168 * t139) * t116 * t127) * t124) * t111) * t137, 0, 0, t105 * t137 * t184 + (-(t102 * t221 + (-t107 * t224 - t119 * t124) * t114) * t228 + (t103 * t201 - t227) * t105 + (t115 * t214 - (-t114 * t223 + t123 * t194 - t220) * t228) * t205) * t111 + (t185 * t214 + ((-t115 * t208 - (-t187 * qJD(4) + t107) * t200) * t151 + (t115 * t186 + (-t159 * t103 - (-t102 - t210) * t225 - (t187 * t107 - qJD(4)) * t222) * t116) * t152) * t111) * t158, (t104 * t228 - t115 * t138) * t204 + (-t104 * t227 + t122 * t115 + (t104 * t201 - t116 * t138) * t103 - (-t151 * t207 - t152 * t206 + t100 * t139 + t110 * t119 + (t110 * t215 + t140) * t107) * t116 * t222 - (-t120 + (-t100 * t158 - t107 * t160) * t152 - t169 * t110) * t200) * t111, 0; t172 * t152 * t203 + (t172 * t209 + ((-qJD(1) * t133 + 0.2e1 * t140 * t226) * t159 + (t120 * t159 - t233 * t122 - t140 * t186) * t134) * t152) * t125, 0, 0, (t133 * t151 * t159 + t160 * t198) * t203 + (0.2e1 * t160 * t180 + (-t176 - t190) * t133 + ((t122 * t159 + 0.2e1 * t157 * t189) * t151 + (t157 * t206 - 0.2e1 * t159 * t178) * t147) * t134) * t125, t134 * t182 * t230 + (t182 * t226 + (t121 * t214 + (-t151 * t208 + t152 * t186) * t137) * t134) * t125, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end