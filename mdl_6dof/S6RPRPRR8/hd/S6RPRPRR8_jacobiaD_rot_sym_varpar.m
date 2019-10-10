% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR8
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
%   Wie in S6RPRPRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:23
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (2081->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->95)
	t139 = cos(qJ(1));
	t203 = 0.2e1 * t139;
	t133 = qJ(3) + pkin(10);
	t131 = sin(t133);
	t132 = cos(t133);
	t181 = t139 * t132;
	t121 = atan2(-t181, t131);
	t119 = sin(t121);
	t120 = cos(t121);
	t105 = -t119 * t181 + t120 * t131;
	t102 = 0.1e1 / t105;
	t136 = sin(qJ(5));
	t180 = t139 * t136;
	t137 = sin(qJ(1));
	t138 = cos(qJ(5));
	t182 = t137 * t138;
	t116 = t131 * t182 + t180;
	t112 = 0.1e1 / t116;
	t126 = 0.1e1 / t131;
	t103 = 0.1e1 / t105 ^ 2;
	t113 = 0.1e1 / t116 ^ 2;
	t127 = 0.1e1 / t131 ^ 2;
	t135 = t139 ^ 2;
	t130 = t132 ^ 2;
	t186 = t127 * t130;
	t124 = t135 * t186 + 0.1e1;
	t122 = 0.1e1 / t124;
	t202 = t122 - 0.1e1;
	t134 = t137 ^ 2;
	t177 = qJD(1) * t139;
	t155 = t130 * t137 * t177;
	t175 = qJD(3) * t132;
	t185 = t130 * t134;
	t174 = qJD(3) * t139;
	t165 = t127 * t174;
	t178 = qJD(1) * t137;
	t166 = t132 * t178;
	t96 = ((t131 * t174 + t166) * t126 + t130 * t165) * t122;
	t160 = -t96 + t174;
	t161 = -t139 * t96 + qJD(3);
	t189 = t120 * t132;
	t91 = t161 * t189 + (t160 * t131 + t166) * t119;
	t199 = t102 * t103 * t91;
	t99 = t103 * t185 + 0.1e1;
	t201 = (-t185 * t199 + (-t131 * t134 * t175 + t155) * t103) / t99 ^ 2;
	t97 = 0.1e1 / t99;
	t200 = t103 * t97;
	t158 = qJD(1) * t131 + qJD(5);
	t151 = t158 * t139;
	t159 = qJD(5) * t131 + qJD(1);
	t152 = t159 * t138;
	t100 = t137 * t152 + (t137 * t175 + t151) * t136;
	t179 = t139 * t138;
	t183 = t137 * t136;
	t115 = t131 * t183 - t179;
	t111 = t115 ^ 2;
	t110 = t111 * t113 + 0.1e1;
	t192 = t113 * t115;
	t153 = t159 * t136;
	t101 = t138 * t151 + (t138 * t175 - t153) * t137;
	t196 = t101 * t112 * t113;
	t198 = 0.1e1 / t110 ^ 2 * (t100 * t192 - t111 * t196);
	t129 = t132 * t130;
	t187 = t126 * t132;
	t149 = qJD(3) * (-t126 * t127 * t129 - t187);
	t194 = (-t127 * t155 + t135 * t149) / t124 ^ 2;
	t193 = t112 * t136;
	t191 = t115 * t138;
	t190 = t119 * t131;
	t188 = t126 * t130;
	t176 = qJD(3) * t131;
	t173 = -0.2e1 * t201;
	t172 = -0.2e1 * t199;
	t171 = 0.2e1 * t198;
	t170 = t102 * t201;
	t169 = t97 * t176;
	t168 = t132 * t194;
	t167 = t122 * t188;
	t164 = 0.1e1 + t186;
	t163 = t194 * t203;
	t162 = 0.2e1 * t115 * t196;
	t157 = t139 * t167;
	t156 = t202 * t132 * t119;
	t154 = t164 * t137;
	t150 = t113 * t191 - t193;
	t148 = t150 * t137;
	t147 = t132 * t174 - t158 * t137;
	t118 = t131 * t179 - t183;
	t117 = t131 * t180 + t182;
	t108 = 0.1e1 / t110;
	t107 = t164 * t139 * t122;
	t95 = (-t120 * t157 - t156) * t137;
	t93 = t139 * t190 + t189 + (-t120 * t181 - t190) * t107;
	t92 = -t164 * t163 + (-qJD(1) * t154 + t149 * t203) * t122;
	t1 = [-0.2e1 * t137 * t126 * t168 + (-qJD(3) * t154 + t177 * t187) * t122, 0, t92, 0, 0, 0; (t102 * t169 + (0.2e1 * t170 + (qJD(1) * t95 + t91) * t200) * t132) * t139 + ((-t95 * t169 + (t95 * t173 + ((t96 * t157 + t202 * t176 + 0.2e1 * t168) * t119 + (t163 * t188 + t132 * t96 + (t129 * t165 + (-t96 + 0.2e1 * t174) * t132) * t122) * t120) * t97 * t137) * t132) * t103 + (t95 * t172 + (t102 + ((t134 - t135) * t120 * t167 - t139 * t156) * t103) * qJD(1)) * t132 * t97) * t137, 0, (t102 * t97 * t177 + (-0.2e1 * t170 + (-qJD(3) * t93 - t91) * t200) * t137) * t131 + (((qJD(3) * t102 + t93 * t172) * t137 + (t93 * t177 + ((t107 * t178 - t139 * t92) * t120 + ((t107 * t139 - 0.1e1) * t96 + (-t107 + t139) * qJD(3)) * t119) * t132 * t137) * t103) * t97 + (t93 * t173 + ((-t92 - t178) * t119 + (t160 * t107 - t161) * t120) * t97 * t131) * t103 * t137) * t132, 0, 0, 0; (-t112 * t117 + t118 * t192) * t171 + (t118 * t162 + t139 * t112 * t152 + t147 * t193 + (t139 * t115 * t153 - t118 * t100 - t117 * t101 - t147 * t191) * t113) * t108, 0, t132 * t148 * t171 + (t148 * t176 + (-t150 * t177 + ((qJD(5) * t112 + t162) * t138 + (-t100 * t138 + (qJD(5) * t115 - t101) * t136) * t113) * t137) * t132) * t108, 0, -0.2e1 * t198 + 0.2e1 * (t100 * t113 * t108 + (-t108 * t196 - t113 * t198) * t115) * t115, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:23
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (2698->94), mult. (2734->203), div. (498->12), fcn. (3199->9), ass. (0->97)
	t162 = qJ(3) + pkin(10);
	t157 = sin(t162);
	t153 = 0.1e1 / t157 ^ 2;
	t158 = cos(t162);
	t156 = t158 ^ 2;
	t209 = t153 * t156;
	t167 = cos(qJ(1));
	t228 = 0.2e1 * t167;
	t227 = t158 * t209;
	t164 = t167 ^ 2;
	t150 = t164 * t209 + 0.1e1;
	t148 = 0.1e1 / t150;
	t152 = 0.1e1 / t157;
	t166 = sin(qJ(1));
	t202 = qJD(1) * t166;
	t191 = t158 * t202;
	t198 = qJD(3) * t167;
	t122 = ((t157 * t198 + t191) * t152 + t198 * t209) * t148;
	t226 = -t122 + t198;
	t204 = t167 * t158;
	t147 = atan2(-t204, t157);
	t145 = sin(t147);
	t146 = cos(t147);
	t131 = -t145 * t204 + t146 * t157;
	t128 = 0.1e1 / t131;
	t165 = qJ(5) + qJ(6);
	t160 = cos(t165);
	t205 = t166 * t160;
	t159 = sin(t165);
	t207 = t159 * t167;
	t142 = t157 * t205 + t207;
	t138 = 0.1e1 / t142;
	t129 = 0.1e1 / t131 ^ 2;
	t139 = 0.1e1 / t142 ^ 2;
	t225 = t148 - 0.1e1;
	t163 = t166 ^ 2;
	t208 = t156 * t163;
	t127 = t129 * t208 + 0.1e1;
	t201 = qJD(1) * t167;
	t183 = t156 * t166 * t201;
	t200 = qJD(3) * t157;
	t212 = t146 * t158;
	t117 = (-t122 * t167 + qJD(3)) * t212 + (t226 * t157 + t191) * t145;
	t223 = t117 * t128 * t129;
	t224 = (-t208 * t223 + (-t158 * t163 * t200 + t183) * t129) / t127 ^ 2;
	t161 = qJD(5) + qJD(6);
	t186 = qJD(1) * t157 + t161;
	t199 = qJD(3) * t166;
	t175 = t158 * t199 + t186 * t167;
	t187 = t157 * t161 + qJD(1);
	t180 = t160 * t187;
	t123 = t175 * t159 + t166 * t180;
	t203 = t167 * t160;
	t206 = t166 * t159;
	t141 = t157 * t206 - t203;
	t137 = t141 ^ 2;
	t136 = t137 * t139 + 0.1e1;
	t215 = t139 * t141;
	t181 = t159 * t187;
	t124 = t175 * t160 - t166 * t181;
	t220 = t124 * t138 * t139;
	t222 = (t123 * t215 - t137 * t220) / t136 ^ 2;
	t221 = t122 * t158;
	t219 = t129 * t158;
	t218 = t129 * t166;
	t210 = t152 * t158;
	t178 = qJD(3) * (-t152 * t227 - t210);
	t217 = (-t153 * t183 + t164 * t178) / t150 ^ 2;
	t216 = t138 * t159;
	t214 = t141 * t160;
	t213 = t145 * t167;
	t211 = t152 * t156;
	t197 = -0.2e1 * t223;
	t196 = 0.2e1 * t222;
	t195 = t158 * t224;
	t194 = t158 * t218;
	t193 = t158 * t217;
	t192 = t148 * t211;
	t190 = 0.1e1 + t209;
	t189 = 0.2e1 * t141 * t220;
	t188 = t217 * t228;
	t185 = t167 * t192;
	t184 = t225 * t158 * t145;
	t182 = t190 * t166;
	t179 = t139 * t214 - t216;
	t177 = t179 * t166;
	t176 = t158 * t198 - t186 * t166;
	t144 = t157 * t203 - t206;
	t143 = t157 * t207 + t205;
	t134 = 0.1e1 / t136;
	t133 = t190 * t167 * t148;
	t125 = 0.1e1 / t127;
	t121 = (-t146 * t185 - t184) * t166;
	t120 = t157 * t213 + t212 + (-t145 * t157 - t146 * t204) * t133;
	t118 = -t190 * t188 + (-qJD(1) * t182 + t178 * t228) * t148;
	t115 = -0.2e1 * t222 + 0.2e1 * (t123 * t134 * t139 + (-t134 * t220 - t139 * t222) * t141) * t141;
	t1 = [-0.2e1 * t166 * t152 * t193 + (-qJD(3) * t182 + t201 * t210) * t148, 0, t118, 0, 0, 0; (0.2e1 * t128 * t195 + (t128 * t200 + (qJD(1) * t121 + t117) * t219) * t125) * t167 + (-0.2e1 * t129 * t195 * t121 + (((t122 * t185 + t225 * t200 + 0.2e1 * t193) * t145 + (t188 * t211 + t221 + (-t221 + (0.2e1 * t158 + t227) * t198) * t148) * t146) * t194 + (-t129 * t200 + t158 * t197) * t121 + (t128 + ((t163 - t164) * t146 * t192 - t167 * t184) * t129) * t158 * qJD(1)) * t125) * t166, 0, 0.2e1 * (-t120 * t219 - t128 * t157) * t166 * t224 + ((t128 * t201 + (-qJD(3) * t120 - t117) * t218) * t157 + (t128 * t199 + (-t118 * t146 * t167 + t226 * t145 + (-qJD(3) * t145 + t122 * t213 + t146 * t202) * t133) * t194 + (t129 * t201 + t166 * t197) * t120 + ((-t118 - t202) * t145 + ((t133 * t167 - 0.1e1) * qJD(3) + (-t133 + t167) * t122) * t146) * t157 * t218) * t158) * t125, 0, 0, 0; (-t138 * t143 + t144 * t215) * t196 + (t144 * t189 + t167 * t138 * t180 + t176 * t216 + (t167 * t141 * t181 - t144 * t123 - t143 * t124 - t176 * t214) * t139) * t134, 0, t158 * t177 * t196 + (t177 * t200 + (-t179 * t201 + ((t138 * t161 + t189) * t160 + (-t123 * t160 + (t141 * t161 - t124) * t159) * t139) * t166) * t158) * t134, 0, t115, t115;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end