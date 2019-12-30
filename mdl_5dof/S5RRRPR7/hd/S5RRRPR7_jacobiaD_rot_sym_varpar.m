% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR7
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
%   Wie in S5RRRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:18
	% EndTime: 2019-12-29 20:05:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:15
	% DurationCPUTime: 1.54s
	% Computational Cost: add. (3413->87), mult. (3414->191), div. (723->12), fcn. (4027->9), ass. (0->91)
	t151 = sin(qJ(1));
	t206 = 0.2e1 * t151;
	t148 = qJ(2) + qJ(3);
	t144 = cos(t148);
	t150 = cos(pkin(9));
	t180 = t151 * t150;
	t149 = sin(pkin(9));
	t152 = cos(qJ(1));
	t184 = t149 * t152;
	t130 = t144 * t184 - t180;
	t124 = t130 ^ 2;
	t181 = t151 * t149;
	t183 = t150 * t152;
	t131 = t144 * t183 + t181;
	t126 = 0.1e1 / t131 ^ 2;
	t121 = t124 * t126 + 0.1e1;
	t128 = -t144 * t181 - t183;
	t143 = sin(t148);
	t145 = qJD(2) + qJD(3);
	t185 = t145 * t152;
	t170 = t143 * t185;
	t122 = qJD(1) * t128 - t149 * t170;
	t196 = t126 * t130;
	t129 = -t144 * t180 + t184;
	t123 = qJD(1) * t129 - t150 * t170;
	t125 = 0.1e1 / t131;
	t198 = t123 * t125 * t126;
	t205 = (t122 * t196 - t124 * t198) / t121 ^ 2;
	t146 = t151 ^ 2;
	t139 = t143 ^ 2;
	t141 = 0.1e1 / t144 ^ 2;
	t192 = t139 * t141;
	t137 = t146 * t192 + 0.1e1;
	t135 = 0.1e1 / t137;
	t140 = 0.1e1 / t144;
	t178 = qJD(1) * t152;
	t169 = t143 * t178;
	t186 = t145 * t151;
	t172 = t141 * t186;
	t109 = (-(-t144 * t186 - t169) * t140 + t139 * t172) * t135;
	t204 = t109 - t186;
	t182 = t151 * t143;
	t134 = atan2(-t182, -t144);
	t133 = cos(t134);
	t132 = sin(t134);
	t173 = t132 * t182;
	t117 = -t133 * t144 - t173;
	t114 = 0.1e1 / t117;
	t115 = 0.1e1 / t117 ^ 2;
	t203 = t135 - 0.1e1;
	t194 = t133 * t143;
	t104 = (-t109 * t151 + t145) * t194 + (t204 * t144 - t169) * t132;
	t202 = t104 * t114 * t115;
	t138 = t143 * t139;
	t189 = t140 * t143;
	t160 = t145 * (t138 * t140 * t141 + t189);
	t190 = t139 * t151;
	t163 = t178 * t190;
	t201 = (t141 * t163 + t146 * t160) / t137 ^ 2;
	t200 = t115 * t143;
	t199 = t115 * t152;
	t197 = t125 * t149;
	t195 = t132 * t151;
	t193 = t139 * t140;
	t147 = t152 ^ 2;
	t191 = t139 * t147;
	t188 = t143 * t152;
	t187 = t144 * t145;
	t179 = qJD(1) * t151;
	t112 = t115 * t191 + 0.1e1;
	t177 = 0.2e1 * (-t191 * t202 + (t143 * t147 * t187 - t163) * t115) / t112 ^ 2;
	t176 = 0.2e1 * t202;
	t175 = t115 * t188;
	t174 = t130 * t198;
	t171 = t145 * t182;
	t168 = 0.1e1 + t192;
	t167 = t143 * t177;
	t166 = -0.2e1 * t143 * t201;
	t165 = t201 * t206;
	t164 = t133 * t135 * t193;
	t162 = t168 * t152;
	t161 = t150 * t196 - t197;
	t119 = 0.1e1 / t121;
	t118 = t168 * t151 * t135;
	t110 = 0.1e1 / t112;
	t108 = (t203 * t143 * t132 - t151 * t164) * t152;
	t106 = -t144 * t195 + t194 + (t132 * t144 - t133 * t182) * t118;
	t105 = -t168 * t165 + (qJD(1) * t162 + t160 * t206) * t135;
	t102 = -0.2e1 * t161 * t188 * t205 + (t161 * t144 * t185 + (-0.2e1 * t174 * t183 + t179 * t197 + (t123 * t184 + (t122 * t152 - t130 * t179) * t150) * t126) * t143) * t119;
	t101 = (t106 * t200 - t114 * t144) * t152 * t177 + ((-t114 * t179 + (-t106 * t145 - t104) * t199) * t144 + (-t114 * t185 - (-t105 * t133 * t151 - t204 * t132 + (t109 * t195 - t132 * t145 - t133 * t178) * t118) * t175 + (t115 * t179 + t152 * t176) * t106 - ((t105 - t178) * t132 + ((-t118 * t151 + 0.1e1) * t145 + (t118 - t151) * t109) * t133) * t144 * t199) * t143) * t110;
	t1 = [t140 * t152 * t166 + (t145 * t162 - t179 * t189) * t135, t105, t105, 0, 0; (t114 * t167 + (-t114 * t187 + (qJD(1) * t108 + t104) * t200) * t110) * t151 + (t115 * t167 * t108 + (-((t166 - t187 + (t109 * t140 * t190 + t187) * t135) * t132 + (t165 * t193 - t109 * t143 + (-t138 * t172 + (t109 - 0.2e1 * t186) * t143) * t135) * t133) * t175 + (-t115 * t187 + t143 * t176) * t108 + (-t114 + ((-t146 + t147) * t164 + t203 * t173) * t115) * t143 * qJD(1)) * t110) * t152, t101, t101, 0, 0; 0.2e1 * (-t125 * t128 + t129 * t196) * t205 + ((-qJD(1) * t130 + t149 * t171) * t125 + 0.2e1 * t129 * t174 + (-t128 * t123 - (-qJD(1) * t131 + t150 * t171) * t130 - t129 * t122) * t126) * t119, t102, t102, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:15
	% DurationCPUTime: 1.65s
	% Computational Cost: add. (4131->97), mult. (3810->205), div. (753->12), fcn. (4455->9), ass. (0->96)
	t172 = sin(qJ(1));
	t233 = 0.2e1 * t172;
	t169 = t172 ^ 2;
	t171 = qJ(2) + qJ(3);
	t165 = sin(t171);
	t159 = t165 ^ 2;
	t166 = cos(t171);
	t161 = 0.1e1 / t166 ^ 2;
	t218 = t159 * t161;
	t155 = t169 * t218 + 0.1e1;
	t152 = 0.1e1 / t155;
	t160 = 0.1e1 / t166;
	t173 = cos(qJ(1));
	t204 = qJD(1) * t173;
	t194 = t165 * t204;
	t168 = qJD(2) + qJD(3);
	t212 = t168 * t172;
	t197 = t161 * t212;
	t126 = (-(-t166 * t212 - t194) * t160 + t159 * t197) * t152;
	t232 = t126 - t212;
	t167 = pkin(9) + qJ(5);
	t164 = cos(t167);
	t206 = t173 * t164;
	t163 = sin(t167);
	t210 = t172 * t163;
	t148 = t166 * t206 + t210;
	t208 = t172 * t165;
	t151 = atan2(-t208, -t166);
	t150 = cos(t151);
	t149 = sin(t151);
	t198 = t149 * t208;
	t136 = -t150 * t166 - t198;
	t133 = 0.1e1 / t136;
	t142 = 0.1e1 / t148;
	t134 = 0.1e1 / t136 ^ 2;
	t143 = 0.1e1 / t148 ^ 2;
	t231 = t152 - 0.1e1;
	t220 = t150 * t165;
	t121 = (-t126 * t172 + t168) * t220 + (t232 * t166 - t194) * t149;
	t230 = t121 * t133 * t134;
	t183 = t166 * t210 + t206;
	t211 = t168 * t173;
	t195 = t165 * t211;
	t127 = t183 * qJD(1) - t148 * qJD(5) + t163 * t195;
	t207 = t173 * t163;
	t209 = t172 * t164;
	t147 = t166 * t207 - t209;
	t141 = t147 ^ 2;
	t139 = t141 * t143 + 0.1e1;
	t223 = t143 * t147;
	t188 = -qJD(1) * t166 + qJD(5);
	t189 = qJD(5) * t166 - qJD(1);
	t128 = -t189 * t207 + (t188 * t172 - t195) * t164;
	t228 = t128 * t142 * t143;
	t229 = (-t127 * t223 - t141 * t228) / t139 ^ 2;
	t158 = t165 * t159;
	t215 = t160 * t165;
	t182 = t168 * (t158 * t160 * t161 + t215);
	t216 = t159 * t172;
	t186 = t204 * t216;
	t227 = (t161 * t186 + t169 * t182) / t155 ^ 2;
	t226 = t134 * t165;
	t225 = t134 * t173;
	t224 = t142 * t163;
	t222 = t147 * t164;
	t221 = t149 * t172;
	t219 = t159 * t160;
	t170 = t173 ^ 2;
	t217 = t159 * t170;
	t214 = t165 * t173;
	t213 = t166 * t168;
	t205 = qJD(1) * t172;
	t131 = t134 * t217 + 0.1e1;
	t203 = 0.2e1 * (-t217 * t230 + (t165 * t170 * t213 - t186) * t134) / t131 ^ 2;
	t202 = 0.2e1 * t230;
	t201 = -0.2e1 * t229;
	t200 = t134 * t214;
	t199 = t147 * t228;
	t193 = 0.1e1 + t218;
	t192 = t165 * t203;
	t191 = -0.2e1 * t165 * t227;
	t190 = t227 * t233;
	t187 = t150 * t152 * t219;
	t185 = t193 * t173;
	t184 = t143 * t222 - t224;
	t181 = t168 * t208 + t188 * t173;
	t146 = -t166 * t209 + t207;
	t140 = t193 * t172 * t152;
	t137 = 0.1e1 / t139;
	t129 = 0.1e1 / t131;
	t125 = (t231 * t165 * t149 - t172 * t187) * t173;
	t124 = -t166 * t221 + t220 + (t149 * t166 - t150 * t208) * t140;
	t122 = -t193 * t190 + (qJD(1) * t185 + t182 * t233) * t152;
	t119 = t184 * t201 * t214 + (t184 * t166 * t211 + (-t184 * t205 + ((-qJD(5) * t142 - 0.2e1 * t199) * t164 + (-t127 * t164 + (-qJD(5) * t147 + t128) * t163) * t143) * t173) * t165) * t137;
	t118 = (t124 * t226 - t133 * t166) * t173 * t203 + ((-t133 * t205 + (-t124 * t168 - t121) * t225) * t166 + (-t133 * t211 - (-t122 * t150 * t172 - t232 * t149 + (t126 * t221 - t149 * t168 - t150 * t204) * t140) * t200 + (t134 * t205 + t173 * t202) * t124 - ((t122 - t204) * t149 + ((-t140 * t172 + 0.1e1) * t168 + (t140 - t172) * t126) * t150) * t166 * t225) * t165) * t129;
	t1 = [t173 * t160 * t191 + (t168 * t185 - t205 * t215) * t152, t122, t122, 0, 0; (t133 * t192 + (-t133 * t213 + (qJD(1) * t125 + t121) * t226) * t129) * t172 + (t134 * t192 * t125 + (-((t191 - t213 + (t126 * t160 * t216 + t213) * t152) * t149 + (t190 * t219 - t126 * t165 + (-t158 * t197 + (t126 - 0.2e1 * t212) * t165) * t152) * t150) * t200 + (-t134 * t213 + t165 * t202) * t125 + (-t133 + ((-t169 + t170) * t187 + t231 * t198) * t134) * t165 * qJD(1)) * t129) * t173, t118, t118, 0, 0; 0.2e1 * (t142 * t183 + t146 * t223) * t229 + (0.2e1 * t146 * t199 - t189 * t142 * t209 + t181 * t224 + (-t189 * t147 * t210 + t146 * t127 + t128 * t183 - t181 * t222) * t143) * t137, t119, t119, 0, t201 + 0.2e1 * (-t127 * t143 * t137 + (-t137 * t228 - t143 * t229) * t147) * t147;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end