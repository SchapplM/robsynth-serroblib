% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR1
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
%   Wie in S6RRPRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:27
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (5265->87), mult. (3414->191), div. (723->12), fcn. (4027->9), ass. (0->91)
	t154 = sin(qJ(1));
	t209 = 0.2e1 * t154;
	t148 = qJ(2) + pkin(10) + qJ(4);
	t147 = cos(t148);
	t152 = sin(pkin(11));
	t155 = cos(qJ(1));
	t184 = t155 * t152;
	t153 = cos(pkin(11));
	t185 = t154 * t153;
	t135 = t147 * t184 - t185;
	t129 = t135 ^ 2;
	t183 = t155 * t153;
	t186 = t154 * t152;
	t136 = t147 * t183 + t186;
	t131 = 0.1e1 / t136 ^ 2;
	t124 = t129 * t131 + 0.1e1;
	t133 = -t147 * t186 - t183;
	t146 = sin(t148);
	t149 = qJD(2) + qJD(4);
	t188 = t149 * t155;
	t173 = t146 * t188;
	t125 = t133 * qJD(1) - t152 * t173;
	t197 = t131 * t135;
	t134 = -t147 * t185 + t184;
	t126 = t134 * qJD(1) - t153 * t173;
	t130 = 0.1e1 / t136;
	t201 = t126 * t130 * t131;
	t208 = (t125 * t197 - t129 * t201) / t124 ^ 2;
	t150 = t154 ^ 2;
	t142 = t146 ^ 2;
	t144 = 0.1e1 / t147 ^ 2;
	t195 = t142 * t144;
	t140 = t150 * t195 + 0.1e1;
	t138 = 0.1e1 / t140;
	t143 = 0.1e1 / t147;
	t181 = qJD(1) * t155;
	t172 = t146 * t181;
	t189 = t149 * t154;
	t175 = t144 * t189;
	t112 = (-(-t147 * t189 - t172) * t143 + t142 * t175) * t138;
	t207 = t112 - t189;
	t187 = t154 * t146;
	t137 = atan2(-t187, -t147);
	t128 = cos(t137);
	t127 = sin(t137);
	t177 = t127 * t187;
	t120 = -t128 * t147 - t177;
	t117 = 0.1e1 / t120;
	t118 = 0.1e1 / t120 ^ 2;
	t206 = t138 - 0.1e1;
	t199 = t128 * t146;
	t107 = (-t112 * t154 + t149) * t199 + (t207 * t147 - t172) * t127;
	t205 = t107 * t117 * t118;
	t141 = t146 * t142;
	t192 = t143 * t146;
	t163 = t149 * (t141 * t143 * t144 + t192);
	t193 = t142 * t154;
	t166 = t181 * t193;
	t204 = (t144 * t166 + t150 * t163) / t140 ^ 2;
	t203 = t118 * t146;
	t202 = t118 * t155;
	t200 = t127 * t154;
	t198 = t130 * t152;
	t196 = t142 * t143;
	t151 = t155 ^ 2;
	t194 = t142 * t151;
	t191 = t146 * t155;
	t190 = t147 * t149;
	t182 = qJD(1) * t154;
	t115 = t118 * t194 + 0.1e1;
	t180 = 0.2e1 * (-t194 * t205 + (t146 * t151 * t190 - t166) * t118) / t115 ^ 2;
	t179 = 0.2e1 * t205;
	t178 = t118 * t191;
	t176 = t135 * t201;
	t174 = t149 * t187;
	t171 = 0.1e1 + t195;
	t170 = t146 * t180;
	t169 = -0.2e1 * t146 * t204;
	t168 = t204 * t209;
	t167 = t128 * t138 * t196;
	t165 = t171 * t155;
	t164 = t153 * t197 - t198;
	t122 = 0.1e1 / t124;
	t121 = t171 * t154 * t138;
	t113 = 0.1e1 / t115;
	t111 = (t206 * t146 * t127 - t154 * t167) * t155;
	t109 = -t147 * t200 + t199 + (t127 * t147 - t128 * t187) * t121;
	t108 = -t171 * t168 + (qJD(1) * t165 + t163 * t209) * t138;
	t105 = -0.2e1 * t164 * t191 * t208 + (t164 * t147 * t188 + (-0.2e1 * t176 * t183 + t182 * t198 + (t126 * t184 + (t125 * t155 - t135 * t182) * t153) * t131) * t146) * t122;
	t104 = (t109 * t203 - t117 * t147) * t155 * t180 + ((-t117 * t182 + (-t109 * t149 - t107) * t202) * t147 + (-t117 * t188 - (-t108 * t128 * t154 - t207 * t127 + (t112 * t200 - t127 * t149 - t128 * t181) * t121) * t178 + (t118 * t182 + t155 * t179) * t109 - ((t108 - t181) * t127 + ((-t121 * t154 + 0.1e1) * t149 + (t121 - t154) * t112) * t128) * t147 * t202) * t146) * t113;
	t1 = [t155 * t143 * t169 + (t149 * t165 - t182 * t192) * t138, t108, 0, t108, 0, 0; (t117 * t170 + (-t117 * t190 + (qJD(1) * t111 + t107) * t203) * t113) * t154 + (t118 * t170 * t111 + (-((t169 - t190 + (t112 * t143 * t193 + t190) * t138) * t127 + (t168 * t196 - t112 * t146 + (-t141 * t175 + (t112 - 0.2e1 * t189) * t146) * t138) * t128) * t178 + (-t118 * t190 + t146 * t179) * t111 + (-t117 + ((-t150 + t151) * t167 + t206 * t177) * t118) * t146 * qJD(1)) * t113) * t155, t104, 0, t104, 0, 0; 0.2e1 * (-t130 * t133 + t134 * t197) * t208 + ((-t135 * qJD(1) + t152 * t174) * t130 + 0.2e1 * t134 * t176 + (-t133 * t126 - (-t136 * qJD(1) + t153 * t174) * t135 - t134 * t125) * t131) * t122, t105, 0, t105, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:27
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (6073->97), mult. (3810->205), div. (753->12), fcn. (4455->9), ass. (0->96)
	t179 = sin(qJ(1));
	t240 = 0.2e1 * t179;
	t177 = t179 ^ 2;
	t174 = qJ(2) + pkin(10) + qJ(4);
	t170 = sin(t174);
	t166 = t170 ^ 2;
	t171 = cos(t174);
	t168 = 0.1e1 / t171 ^ 2;
	t225 = t166 * t168;
	t161 = t177 * t225 + 0.1e1;
	t159 = 0.1e1 / t161;
	t167 = 0.1e1 / t171;
	t180 = cos(qJ(1));
	t211 = qJD(1) * t180;
	t201 = t170 * t211;
	t176 = qJD(2) + qJD(4);
	t217 = t176 * t179;
	t204 = t168 * t217;
	t133 = (-(-t171 * t217 - t201) * t167 + t166 * t204) * t159;
	t239 = t133 - t217;
	t175 = pkin(11) + qJ(6);
	t172 = sin(t175);
	t214 = t179 * t172;
	t173 = cos(t175);
	t218 = t173 * t180;
	t155 = t171 * t218 + t214;
	t215 = t179 * t170;
	t158 = atan2(-t215, -t171);
	t157 = cos(t158);
	t156 = sin(t158);
	t205 = t156 * t215;
	t143 = -t157 * t171 - t205;
	t140 = 0.1e1 / t143;
	t149 = 0.1e1 / t155;
	t141 = 0.1e1 / t143 ^ 2;
	t150 = 0.1e1 / t155 ^ 2;
	t238 = t159 - 0.1e1;
	t227 = t157 * t170;
	t128 = (-t133 * t179 + t176) * t227 + (t239 * t171 - t201) * t156;
	t237 = t128 * t140 * t141;
	t190 = t171 * t214 + t218;
	t216 = t176 * t180;
	t202 = t170 * t216;
	t134 = t190 * qJD(1) - t155 * qJD(6) + t172 * t202;
	t213 = t179 * t173;
	t219 = t172 * t180;
	t154 = t171 * t219 - t213;
	t148 = t154 ^ 2;
	t147 = t148 * t150 + 0.1e1;
	t230 = t150 * t154;
	t195 = -qJD(1) * t171 + qJD(6);
	t196 = qJD(6) * t171 - qJD(1);
	t135 = -t196 * t219 + (t195 * t179 - t202) * t173;
	t235 = t135 * t149 * t150;
	t236 = (-t134 * t230 - t148 * t235) / t147 ^ 2;
	t165 = t170 * t166;
	t222 = t167 * t170;
	t189 = t176 * (t165 * t167 * t168 + t222);
	t223 = t166 * t179;
	t193 = t211 * t223;
	t234 = (t168 * t193 + t177 * t189) / t161 ^ 2;
	t233 = t141 * t170;
	t232 = t141 * t180;
	t231 = t149 * t172;
	t229 = t154 * t173;
	t228 = t156 * t179;
	t226 = t166 * t167;
	t178 = t180 ^ 2;
	t224 = t166 * t178;
	t221 = t170 * t180;
	t220 = t171 * t176;
	t212 = qJD(1) * t179;
	t138 = t141 * t224 + 0.1e1;
	t210 = 0.2e1 * (-t224 * t237 + (t170 * t178 * t220 - t193) * t141) / t138 ^ 2;
	t209 = 0.2e1 * t237;
	t208 = -0.2e1 * t236;
	t207 = t141 * t221;
	t206 = t154 * t235;
	t200 = 0.1e1 + t225;
	t199 = t170 * t210;
	t198 = -0.2e1 * t170 * t234;
	t197 = t234 * t240;
	t194 = t157 * t159 * t226;
	t192 = t200 * t180;
	t191 = t150 * t229 - t231;
	t188 = t176 * t215 + t195 * t180;
	t153 = -t171 * t213 + t219;
	t145 = 0.1e1 / t147;
	t144 = t200 * t179 * t159;
	t136 = 0.1e1 / t138;
	t132 = (t238 * t170 * t156 - t179 * t194) * t180;
	t131 = -t171 * t228 + t227 + (t156 * t171 - t157 * t215) * t144;
	t129 = -t200 * t197 + (qJD(1) * t192 + t189 * t240) * t159;
	t126 = t191 * t208 * t221 + (t191 * t171 * t216 + (-t191 * t212 + ((-qJD(6) * t149 - 0.2e1 * t206) * t173 + (-t134 * t173 + (-qJD(6) * t154 + t135) * t172) * t150) * t180) * t170) * t145;
	t125 = (t131 * t233 - t140 * t171) * t180 * t210 + ((-t140 * t212 + (-t131 * t176 - t128) * t232) * t171 + (-t140 * t216 - (-t129 * t157 * t179 - t239 * t156 + (t133 * t228 - t156 * t176 - t157 * t211) * t144) * t207 + (t141 * t212 + t180 * t209) * t131 - ((t129 - t211) * t156 + ((-t144 * t179 + 0.1e1) * t176 + (t144 - t179) * t133) * t157) * t171 * t232) * t170) * t136;
	t1 = [t167 * t180 * t198 + (t176 * t192 - t212 * t222) * t159, t129, 0, t129, 0, 0; (t140 * t199 + (-t140 * t220 + (qJD(1) * t132 + t128) * t233) * t136) * t179 + (t141 * t199 * t132 + (-((t198 - t220 + (t133 * t167 * t223 + t220) * t159) * t156 + (t197 * t226 - t133 * t170 + (-t165 * t204 + (t133 - 0.2e1 * t217) * t170) * t159) * t157) * t207 + (-t141 * t220 + t170 * t209) * t132 + (-t140 + ((-t177 + t178) * t194 + t238 * t205) * t141) * t170 * qJD(1)) * t136) * t180, t125, 0, t125, 0, 0; 0.2e1 * (t149 * t190 + t153 * t230) * t236 + (0.2e1 * t153 * t206 - t196 * t149 * t213 + t188 * t231 + (-t196 * t154 * t214 + t153 * t134 + t135 * t190 - t188 * t229) * t150) * t145, t126, 0, t126, 0, t208 + 0.2e1 * (-t134 * t150 * t145 + (-t145 * t235 - t150 * t236) * t154) * t154;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end