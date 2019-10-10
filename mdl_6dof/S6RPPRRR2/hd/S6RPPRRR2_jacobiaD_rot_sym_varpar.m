% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR2
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
%   Wie in S6RPPRRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:03
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:06
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:07
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (3244->94), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->93)
	t144 = qJ(1) + pkin(10);
	t140 = sin(t144);
	t134 = t140 ^ 2;
	t143 = pkin(11) + qJ(4);
	t139 = sin(t143);
	t133 = t139 ^ 2;
	t141 = cos(t143);
	t136 = 0.1e1 / t141 ^ 2;
	t193 = t133 * t136;
	t128 = t134 * t193 + 0.1e1;
	t132 = t139 * t133;
	t135 = 0.1e1 / t141;
	t190 = t135 * t139;
	t154 = qJD(4) * (t132 * t135 * t136 + t190);
	t142 = cos(t144);
	t181 = qJD(1) * t142;
	t191 = t133 * t140;
	t159 = t181 * t191;
	t197 = (t134 * t154 + t136 * t159) / t128 ^ 2;
	t207 = -0.2e1 * t197;
	t165 = 0.1e1 + t193;
	t206 = t140 * t165;
	t146 = cos(qJ(5));
	t184 = t142 * t146;
	t145 = sin(qJ(5));
	t187 = t140 * t145;
	t122 = t141 * t184 + t187;
	t162 = qJD(5) * t141 - qJD(1);
	t180 = qJD(4) * t139;
	t205 = t162 * t145 + t146 * t180;
	t188 = t140 * t139;
	t125 = atan2(-t188, -t141);
	t124 = cos(t125);
	t123 = sin(t125);
	t173 = t123 * t188;
	t109 = -t124 * t141 - t173;
	t106 = 0.1e1 / t109;
	t116 = 0.1e1 / t122;
	t107 = 0.1e1 / t109 ^ 2;
	t117 = 0.1e1 / t122 ^ 2;
	t126 = 0.1e1 / t128;
	t204 = t126 - 0.1e1;
	t172 = t126 * t133 * t135;
	t160 = t140 * t172;
	t99 = (t204 * t139 * t123 - t124 * t160) * t142;
	t203 = t107 * t99;
	t179 = qJD(4) * t140;
	t169 = t136 * t179;
	t170 = t139 * t181;
	t178 = qJD(4) * t141;
	t100 = (-(-t140 * t178 - t170) * t135 + t133 * t169) * t126;
	t194 = t124 * t139;
	t95 = (-t100 * t140 + qJD(4)) * t194 + (-t170 + (t100 - t179) * t141) * t123;
	t202 = t106 * t107 * t95;
	t155 = t141 * t187 + t184;
	t168 = t145 * t180;
	t104 = t155 * qJD(1) - t122 * qJD(5) + t142 * t168;
	t185 = t142 * t145;
	t186 = t140 * t146;
	t121 = t141 * t185 - t186;
	t115 = t121 ^ 2;
	t114 = t115 * t117 + 0.1e1;
	t196 = t117 * t121;
	t161 = -qJD(1) * t141 + qJD(5);
	t157 = t161 * t146;
	t105 = t140 * t157 - t205 * t142;
	t200 = t105 * t116 * t117;
	t201 = 0.1e1 / t114 ^ 2 * (-t104 * t196 - t115 * t200);
	t199 = t107 * t139;
	t198 = t107 * t142;
	t195 = t123 * t141;
	t138 = t142 ^ 2;
	t192 = t133 * t138;
	t189 = t139 * t142;
	t111 = t126 * t206;
	t183 = t111 - t140;
	t182 = qJD(1) * t140;
	t103 = t107 * t192 + 0.1e1;
	t177 = 0.2e1 / t103 ^ 2 * (-t192 * t202 + (t138 * t139 * t178 - t159) * t107);
	t176 = 0.2e1 * t202;
	t175 = -0.2e1 * t201;
	t174 = t121 * t200;
	t166 = t111 * t140 - 0.1e1;
	t164 = t139 * t177;
	t163 = t135 * t207;
	t158 = t165 * t142;
	t156 = -t116 * t145 + t146 * t196;
	t120 = -t141 * t186 + t185;
	t112 = 0.1e1 / t114;
	t101 = 0.1e1 / t103;
	t97 = -t140 * t195 + t194 + (-t124 * t188 + t195) * t111;
	t96 = t206 * t207 + (qJD(1) * t158 + 0.2e1 * t140 * t154) * t126;
	t1 = [t163 * t189 + (qJD(4) * t158 - t182 * t190) * t126, 0, 0, t96, 0, 0; (t106 * t164 + (-t106 * t178 + (qJD(1) * t99 + t95) * t199) * t101) * t140 + (t164 * t203 + (-t178 * t203 + (t99 * t176 + ((-t100 * t160 + 0.2e1 * t139 * t197 - t204 * t178) * t123 + (t163 * t191 + t100 * t139 + (t132 * t169 - (t100 - 0.2e1 * t179) * t139) * t126) * t124) * t198) * t139 + (-t106 + ((-t134 + t138) * t124 * t172 + t204 * t173) * t107) * t139 * qJD(1)) * t101) * t142, 0, 0, (-t106 * t141 + t97 * t199) * t142 * t177 + ((-t106 * t182 + (-qJD(4) * t97 - t95) * t198) * t141 + ((-qJD(4) * t106 + t97 * t176) * t142 + (t97 * t182 + (-(-t111 * t181 - t140 * t96) * t124 - (-t183 * qJD(4) + t166 * t100) * t123) * t189) * t107 - ((t96 - t181) * t123 + (-t166 * qJD(4) + t183 * t100) * t124) * t141 * t198) * t139) * t101, 0, 0; 0.2e1 * (t116 * t155 + t120 * t196) * t201 + (0.2e1 * t120 * t174 + (t120 * t104 + t155 * t105 + (-t205 * t140 - t142 * t157) * t121) * t117 + (t161 * t185 + (-t162 * t146 + t168) * t140) * t116) * t112, 0, 0, t156 * t175 * t189 + (t156 * t142 * t178 + (-t156 * t182 + ((-qJD(5) * t116 - 0.2e1 * t174) * t146 + (-t104 * t146 + (-qJD(5) * t121 + t105) * t145) * t117) * t142) * t139) * t112, t175 + 0.2e1 * (-t104 * t112 * t117 + (-t112 * t200 - t117 * t201) * t121) * t121, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:03:06
	% EndTime: 2019-10-10 00:03:07
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (3952->98), mult. (2734->203), div. (498->12), fcn. (3199->9), ass. (0->98)
	t172 = pkin(11) + qJ(4);
	t166 = sin(t172);
	t160 = t166 ^ 2;
	t168 = cos(t172);
	t163 = 0.1e1 / t168 ^ 2;
	t221 = t160 * t163;
	t174 = qJ(1) + pkin(10);
	t167 = sin(t174);
	t239 = 0.2e1 * t167;
	t238 = t166 * t221;
	t169 = cos(t174);
	t175 = qJ(5) + qJ(6);
	t171 = cos(t175);
	t213 = t169 * t171;
	t170 = sin(t175);
	t216 = t167 * t170;
	t149 = t168 * t213 + t216;
	t173 = qJD(5) + qJD(6);
	t191 = t168 * t173 - qJD(1);
	t210 = qJD(4) * t166;
	t237 = t191 * t170 + t171 * t210;
	t217 = t167 * t166;
	t152 = atan2(-t217, -t168);
	t151 = cos(t152);
	t150 = sin(t152);
	t201 = t150 * t217;
	t136 = -t151 * t168 - t201;
	t133 = 0.1e1 / t136;
	t143 = 0.1e1 / t149;
	t162 = 0.1e1 / t168;
	t134 = 0.1e1 / t136 ^ 2;
	t144 = 0.1e1 / t149 ^ 2;
	t236 = -0.2e1 * t166;
	t161 = t167 ^ 2;
	t156 = t161 * t221 + 0.1e1;
	t154 = 0.1e1 / t156;
	t235 = t154 - 0.1e1;
	t211 = qJD(1) * t169;
	t198 = t166 * t211;
	t208 = qJD(4) * t168;
	t209 = qJD(4) * t167;
	t127 = (-(-t167 * t208 - t198) * t162 + t209 * t221) * t154;
	t223 = t151 * t166;
	t122 = (-t127 * t167 + qJD(4)) * t223 + (-t198 + (t127 - t209) * t168) * t150;
	t234 = t122 * t133 * t134;
	t184 = t168 * t216 + t213;
	t197 = t170 * t210;
	t128 = t184 * qJD(1) - t149 * t173 + t169 * t197;
	t214 = t169 * t170;
	t215 = t167 * t171;
	t148 = t168 * t214 - t215;
	t142 = t148 ^ 2;
	t141 = t142 * t144 + 0.1e1;
	t225 = t144 * t148;
	t190 = -qJD(1) * t168 + t173;
	t186 = t190 * t171;
	t129 = t167 * t186 - t169 * t237;
	t230 = t129 * t143 * t144;
	t233 = (-t128 * t225 - t142 * t230) / t141 ^ 2;
	t232 = t127 * t150;
	t231 = t127 * t166;
	t229 = t134 * t166;
	t228 = t134 * t169;
	t219 = t162 * t166;
	t183 = qJD(4) * (t162 * t238 + t219);
	t188 = t160 * t167 * t211;
	t227 = (t161 * t183 + t163 * t188) / t156 ^ 2;
	t195 = 0.1e1 + t221;
	t138 = t195 * t167 * t154;
	t226 = t138 * t167;
	t224 = t150 * t168;
	t222 = t160 * t162;
	t165 = t169 ^ 2;
	t220 = t160 * t165;
	t218 = t166 * t169;
	t212 = qJD(1) * t167;
	t207 = qJD(4) * t169;
	t132 = t134 * t220 + 0.1e1;
	t206 = 0.2e1 * (-t220 * t234 + (t165 * t166 * t208 - t188) * t134) / t132 ^ 2;
	t205 = 0.2e1 * t234;
	t204 = -0.2e1 * t233;
	t203 = t148 * t230;
	t202 = t134 * t218;
	t200 = t154 * t222;
	t194 = t166 * t206;
	t193 = t227 * t236;
	t192 = t227 * t239;
	t189 = t167 * t200;
	t187 = t195 * t169;
	t185 = -t143 * t170 + t171 * t225;
	t147 = -t168 * t215 + t214;
	t139 = 0.1e1 / t141;
	t130 = 0.1e1 / t132;
	t126 = (t235 * t166 * t150 - t151 * t189) * t169;
	t125 = -t167 * t224 + t223 + (-t151 * t217 + t224) * t138;
	t123 = -t195 * t192 + (qJD(1) * t187 + t183 * t239) * t154;
	t120 = t204 + 0.2e1 * (-t128 * t139 * t144 + (-t139 * t230 - t144 * t233) * t148) * t148;
	t1 = [t162 * t169 * t193 + (qJD(4) * t187 - t212 * t219) * t154, 0, 0, t123, 0, 0; (t133 * t194 + (-t133 * t208 + (qJD(1) * t126 + t122) * t229) * t130) * t167 + (t134 * t194 * t126 + (-((t127 * t189 + t235 * t208 + t193) * t150 + (t192 * t222 - t231 + (t231 + (t236 - t238) * t209) * t154) * t151) * t202 + (-t134 * t208 + t166 * t205) * t126 + (-t133 + ((-t161 + t165) * t151 * t200 + t235 * t201) * t134) * t166 * qJD(1)) * t130) * t169, 0, 0, (t125 * t229 - t133 * t168) * t169 * t206 + ((-t133 * t212 + (-qJD(4) * t125 - t122) * t228) * t168 + (-t133 * t207 - (-t123 * t151 * t167 + t150 * t209 + t226 * t232 - t232 + (-qJD(4) * t150 - t151 * t211) * t138) * t202 + (t134 * t212 + t169 * t205) * t125 - ((t123 - t211) * t150 + ((0.1e1 - t226) * qJD(4) + (t138 - t167) * t127) * t151) * t168 * t228) * t166) * t130, 0, 0; 0.2e1 * (t143 * t184 + t147 * t225) * t233 + (0.2e1 * t147 * t203 + (t147 * t128 + t184 * t129 + (-t167 * t237 - t169 * t186) * t148) * t144 + (t190 * t214 + (-t191 * t171 + t197) * t167) * t143) * t139, 0, 0, t185 * t204 * t218 + (t185 * t168 * t207 + (-t185 * t212 + ((-t143 * t173 - 0.2e1 * t203) * t171 + (-t128 * t171 + (-t148 * t173 + t129) * t170) * t144) * t169) * t166) * t139, t120, t120;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end