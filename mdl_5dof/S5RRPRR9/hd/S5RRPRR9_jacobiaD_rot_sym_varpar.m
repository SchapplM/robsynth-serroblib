% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRPRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPRR9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:32
	% EndTime: 2019-12-29 19:10:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:32
	% EndTime: 2019-12-29 19:10:32
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:32
	% EndTime: 2019-12-29 19:10:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:27
	% EndTime: 2019-12-29 19:10:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:27
	% EndTime: 2019-12-29 19:10:29
	% DurationCPUTime: 1.55s
	% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
	t139 = sin(qJ(1));
	t136 = t139 ^ 2;
	t135 = qJ(2) + pkin(9);
	t133 = sin(t135);
	t129 = t133 ^ 2;
	t134 = cos(t135);
	t131 = 0.1e1 / t134 ^ 2;
	t188 = t129 * t131;
	t124 = t136 * t188 + 0.1e1;
	t128 = t133 * t129;
	t130 = 0.1e1 / t134;
	t185 = t130 * t133;
	t149 = qJD(2) * (t128 * t130 * t131 + t185);
	t141 = cos(qJ(1));
	t177 = qJD(1) * t141;
	t186 = t129 * t139;
	t154 = t177 * t186;
	t194 = (t131 * t154 + t136 * t149) / t124 ^ 2;
	t204 = -0.2e1 * t194;
	t161 = 0.1e1 + t188;
	t203 = t139 * t161;
	t140 = cos(qJ(4));
	t179 = t141 * t140;
	t138 = sin(qJ(4));
	t182 = t139 * t138;
	t120 = t134 * t179 + t182;
	t183 = t139 * t133;
	t121 = atan2(-t183, -t134);
	t116 = cos(t121);
	t115 = sin(t121);
	t168 = t115 * t183;
	t105 = -t116 * t134 - t168;
	t102 = 0.1e1 / t105;
	t112 = 0.1e1 / t120;
	t103 = 0.1e1 / t105 ^ 2;
	t113 = 0.1e1 / t120 ^ 2;
	t122 = 0.1e1 / t124;
	t202 = t122 - 0.1e1;
	t137 = t141 ^ 2;
	t176 = qJD(2) * t134;
	t187 = t129 * t137;
	t175 = qJD(2) * t139;
	t163 = t131 * t175;
	t164 = t133 * t177;
	t96 = (-(-t134 * t175 - t164) * t130 + t129 * t163) * t122;
	t158 = t96 - t175;
	t159 = -t139 * t96 + qJD(2);
	t190 = t116 * t133;
	t91 = t159 * t190 + (t158 * t134 - t164) * t115;
	t199 = t102 * t103 * t91;
	t99 = t103 * t187 + 0.1e1;
	t201 = (-t187 * t199 + (t133 * t137 * t176 - t154) * t103) / t99 ^ 2;
	t97 = 0.1e1 / t99;
	t200 = t103 * t97;
	t150 = t134 * t182 + t179;
	t174 = qJD(2) * t141;
	t162 = t133 * t174;
	t100 = t150 * qJD(1) - t120 * qJD(4) + t138 * t162;
	t180 = t141 * t138;
	t181 = t139 * t140;
	t119 = t134 * t180 - t181;
	t111 = t119 ^ 2;
	t110 = t111 * t113 + 0.1e1;
	t192 = t113 * t119;
	t156 = -qJD(1) * t134 + qJD(4);
	t157 = qJD(4) * t134 - qJD(1);
	t101 = -t157 * t180 + (t156 * t139 - t162) * t140;
	t196 = t101 * t112 * t113;
	t198 = 0.1e1 / t110 ^ 2 * (-t100 * t192 - t111 * t196);
	t193 = t112 * t138;
	t191 = t115 * t134;
	t189 = t119 * t140;
	t184 = t133 * t141;
	t178 = qJD(1) * t139;
	t173 = 0.2e1 * t201;
	t172 = 0.2e1 * t199;
	t171 = -0.2e1 * t198;
	t170 = t102 * t201;
	t169 = t97 * t176;
	t167 = t119 * t196;
	t166 = t122 * t129 * t130;
	t160 = t130 * t204;
	t155 = t139 * t166;
	t153 = t161 * t141;
	t152 = t156 * t141;
	t151 = t113 * t189 - t193;
	t118 = -t134 * t181 + t180;
	t108 = 0.1e1 / t110;
	t107 = t122 * t203;
	t95 = (t202 * t133 * t115 - t116 * t155) * t141;
	t93 = -t139 * t191 + t190 + (-t116 * t183 + t191) * t107;
	t92 = t203 * t204 + (qJD(1) * t153 + 0.2e1 * t139 * t149) * t122;
	t1 = [t160 * t184 + (qJD(2) * t153 - t178 * t185) * t122, t92, 0, 0, 0; (-t102 * t169 + (0.2e1 * t170 + (qJD(1) * t95 + t91) * t200) * t133) * t139 + ((-t95 * t169 + (t95 * t173 + ((0.2e1 * t133 * t194 - t96 * t155 - t202 * t176) * t115 + (t160 * t186 + t133 * t96 + (t128 * t163 - (t96 - 0.2e1 * t175) * t133) * t122) * t116) * t97 * t141) * t133) * t103 + (t95 * t172 + (-t102 + ((-t136 + t137) * t116 * t166 + t202 * t168) * t103) * qJD(1)) * t133 * t97) * t141, (-t102 * t97 * t178 + (-0.2e1 * t170 + (-qJD(2) * t93 - t91) * t200) * t141) * t134 + (((-qJD(2) * t102 + t93 * t172) * t141 + (t93 * t178 + (-(-t107 * t177 - t139 * t92) * t116 - ((t107 * t139 - 0.1e1) * t96 + (-t107 + t139) * qJD(2)) * t115) * t184) * t103) * t97 + (t93 * t173 - ((t92 - t177) * t115 + (t158 * t107 + t159) * t116) * t97 * t134) * t103 * t141) * t133, 0, 0, 0; 0.2e1 * (t112 * t150 + t118 * t192) * t198 + (0.2e1 * t118 * t167 - t157 * t112 * t181 + (t133 * t175 + t152) * t193 + (t118 * t100 + t150 * t101 - t152 * t189 - (qJD(2) * t133 * t140 + t157 * t138) * t119 * t139) * t113) * t108, t151 * t171 * t184 + (t151 * t134 * t174 + (-t151 * t178 + ((-qJD(4) * t112 - 0.2e1 * t167) * t140 + (-t100 * t140 + (-qJD(4) * t119 + t101) * t138) * t113) * t141) * t133) * t108, 0, t171 + 0.2e1 * (-t100 * t113 * t108 + (-t108 * t196 - t113 * t198) * t119) * t119, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:10:27
	% EndTime: 2019-12-29 19:10:29
	% DurationCPUTime: 1.61s
	% Computational Cost: add. (2887->96), mult. (2734->202), div. (498->12), fcn. (3199->9), ass. (0->96)
	t165 = qJ(2) + pkin(9);
	t160 = sin(t165);
	t156 = t160 ^ 2;
	t161 = cos(t165);
	t158 = 0.1e1 / t161 ^ 2;
	t214 = t156 * t158;
	t169 = sin(qJ(1));
	t232 = 0.2e1 * t169;
	t231 = t160 * t214;
	t166 = t169 ^ 2;
	t152 = t166 * t214 + 0.1e1;
	t150 = 0.1e1 / t152;
	t157 = 0.1e1 / t161;
	t170 = cos(qJ(1));
	t204 = qJD(1) * t170;
	t192 = t160 * t204;
	t202 = qJD(2) * t169;
	t123 = (-(-t161 * t202 - t192) * t157 + t202 * t214) * t150;
	t230 = t123 - t202;
	t168 = qJ(4) + qJ(5);
	t163 = cos(t168);
	t206 = t170 * t163;
	t162 = sin(t168);
	t209 = t169 * t162;
	t145 = t161 * t206 + t209;
	t210 = t169 * t160;
	t148 = atan2(-t210, -t161);
	t147 = cos(t148);
	t146 = sin(t148);
	t195 = t146 * t210;
	t132 = -t147 * t161 - t195;
	t129 = 0.1e1 / t132;
	t139 = 0.1e1 / t145;
	t130 = 0.1e1 / t132 ^ 2;
	t140 = 0.1e1 / t145 ^ 2;
	t229 = -0.2e1 * t160;
	t228 = t150 - 0.1e1;
	t216 = t147 * t160;
	t118 = (-t123 * t169 + qJD(2)) * t216 + (t230 * t161 - t192) * t146;
	t227 = t118 * t129 * t130;
	t164 = qJD(4) + qJD(5);
	t180 = t161 * t209 + t206;
	t201 = qJD(2) * t170;
	t191 = t160 * t201;
	t124 = qJD(1) * t180 - t145 * t164 + t162 * t191;
	t207 = t170 * t162;
	t208 = t169 * t163;
	t144 = t161 * t207 - t208;
	t138 = t144 ^ 2;
	t137 = t138 * t140 + 0.1e1;
	t219 = t140 * t144;
	t185 = -qJD(1) * t161 + t164;
	t186 = t161 * t164 - qJD(1);
	t125 = -t186 * t207 + (t169 * t185 - t191) * t163;
	t224 = t125 * t139 * t140;
	t226 = (-t124 * t219 - t138 * t224) / t137 ^ 2;
	t225 = t123 * t160;
	t223 = t130 * t160;
	t222 = t130 * t170;
	t212 = t157 * t160;
	t179 = qJD(2) * (t157 * t231 + t212);
	t183 = t156 * t169 * t204;
	t221 = (t158 * t183 + t166 * t179) / t152 ^ 2;
	t220 = t139 * t162;
	t218 = t144 * t163;
	t217 = t146 * t169;
	t215 = t156 * t157;
	t167 = t170 ^ 2;
	t213 = t156 * t167;
	t211 = t160 * t170;
	t205 = qJD(1) * t169;
	t203 = qJD(2) * t161;
	t128 = t130 * t213 + 0.1e1;
	t200 = 0.2e1 * (-t213 * t227 + (t160 * t167 * t203 - t183) * t130) / t128 ^ 2;
	t199 = 0.2e1 * t227;
	t198 = -0.2e1 * t226;
	t197 = t130 * t211;
	t196 = t144 * t224;
	t194 = t150 * t215;
	t190 = 0.1e1 + t214;
	t189 = t160 * t200;
	t188 = t221 * t229;
	t187 = t221 * t232;
	t184 = t169 * t194;
	t182 = t190 * t170;
	t181 = t218 * t140 - t220;
	t178 = t160 * t202 + t170 * t185;
	t143 = -t161 * t208 + t207;
	t135 = 0.1e1 / t137;
	t134 = t190 * t169 * t150;
	t126 = 0.1e1 / t128;
	t122 = (t146 * t160 * t228 - t147 * t184) * t170;
	t121 = -t161 * t217 + t216 + (t146 * t161 - t147 * t210) * t134;
	t119 = -t190 * t187 + (qJD(1) * t182 + t179 * t232) * t150;
	t116 = t198 + 0.2e1 * (-t124 * t140 * t135 + (-t135 * t224 - t140 * t226) * t144) * t144;
	t1 = [t170 * t157 * t188 + (qJD(2) * t182 - t205 * t212) * t150, t119, 0, 0, 0; (t129 * t189 + (-t129 * t203 + (qJD(1) * t122 + t118) * t223) * t126) * t169 + (t130 * t189 * t122 + (-((t123 * t184 + t203 * t228 + t188) * t146 + (t187 * t215 - t225 + (t225 + (t229 - t231) * t202) * t150) * t147) * t197 + (-t130 * t203 + t160 * t199) * t122 + (-t129 + ((-t166 + t167) * t147 * t194 + t228 * t195) * t130) * t160 * qJD(1)) * t126) * t170, (t121 * t223 - t129 * t161) * t170 * t200 + ((-t129 * t205 + (-qJD(2) * t121 - t118) * t222) * t161 + (-t129 * t201 - (-t119 * t147 * t169 - t230 * t146 + (-qJD(2) * t146 + t123 * t217 - t147 * t204) * t134) * t197 + (t130 * t205 + t170 * t199) * t121 - ((t119 - t204) * t146 + ((-t134 * t169 + 0.1e1) * qJD(2) + (t134 - t169) * t123) * t147) * t161 * t222) * t160) * t126, 0, 0, 0; 0.2e1 * (t139 * t180 + t143 * t219) * t226 + (0.2e1 * t143 * t196 - t186 * t139 * t208 + t178 * t220 + (-t144 * t186 * t209 + t143 * t124 + t125 * t180 - t178 * t218) * t140) * t135, t181 * t198 * t211 + (t181 * t161 * t201 + (-t181 * t205 + ((-t139 * t164 - 0.2e1 * t196) * t163 + (-t124 * t163 + (-t144 * t164 + t125) * t162) * t140) * t170) * t160) * t135, 0, t116, t116;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end