% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR1
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
%   Wie in S6RRRPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:17
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:05
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (5265->87), mult. (3414->191), div. (723->12), fcn. (4027->9), ass. (0->91)
	t156 = sin(qJ(1));
	t211 = 0.2e1 * t156;
	t150 = qJ(2) + qJ(3) + pkin(10);
	t149 = cos(t150);
	t154 = sin(pkin(11));
	t157 = cos(qJ(1));
	t186 = t157 * t154;
	t155 = cos(pkin(11));
	t187 = t156 * t155;
	t137 = t149 * t186 - t187;
	t131 = t137 ^ 2;
	t185 = t157 * t155;
	t188 = t156 * t154;
	t138 = t149 * t185 + t188;
	t133 = 0.1e1 / t138 ^ 2;
	t126 = t131 * t133 + 0.1e1;
	t135 = -t149 * t188 - t185;
	t148 = sin(t150);
	t151 = qJD(2) + qJD(3);
	t190 = t151 * t157;
	t175 = t148 * t190;
	t127 = t135 * qJD(1) - t154 * t175;
	t199 = t133 * t137;
	t136 = -t149 * t187 + t186;
	t128 = t136 * qJD(1) - t155 * t175;
	t132 = 0.1e1 / t138;
	t203 = t128 * t132 * t133;
	t210 = (t127 * t199 - t131 * t203) / t126 ^ 2;
	t152 = t156 ^ 2;
	t144 = t148 ^ 2;
	t146 = 0.1e1 / t149 ^ 2;
	t197 = t144 * t146;
	t142 = t152 * t197 + 0.1e1;
	t140 = 0.1e1 / t142;
	t145 = 0.1e1 / t149;
	t183 = qJD(1) * t157;
	t174 = t148 * t183;
	t191 = t151 * t156;
	t177 = t146 * t191;
	t114 = (-(-t149 * t191 - t174) * t145 + t144 * t177) * t140;
	t209 = t114 - t191;
	t189 = t156 * t148;
	t139 = atan2(-t189, -t149);
	t130 = cos(t139);
	t129 = sin(t139);
	t178 = t129 * t189;
	t122 = -t130 * t149 - t178;
	t119 = 0.1e1 / t122;
	t120 = 0.1e1 / t122 ^ 2;
	t208 = t140 - 0.1e1;
	t201 = t130 * t148;
	t109 = (-t114 * t156 + t151) * t201 + (t209 * t149 - t174) * t129;
	t207 = t109 * t119 * t120;
	t143 = t148 * t144;
	t194 = t145 * t148;
	t165 = t151 * (t143 * t145 * t146 + t194);
	t195 = t144 * t156;
	t168 = t183 * t195;
	t206 = (t146 * t168 + t152 * t165) / t142 ^ 2;
	t205 = t120 * t148;
	t204 = t120 * t157;
	t202 = t129 * t156;
	t200 = t132 * t154;
	t198 = t144 * t145;
	t153 = t157 ^ 2;
	t196 = t144 * t153;
	t193 = t148 * t157;
	t192 = t149 * t151;
	t184 = qJD(1) * t156;
	t117 = t120 * t196 + 0.1e1;
	t182 = 0.2e1 * (-t196 * t207 + (t148 * t153 * t192 - t168) * t120) / t117 ^ 2;
	t181 = 0.2e1 * t207;
	t180 = t120 * t193;
	t179 = t137 * t203;
	t176 = t151 * t189;
	t173 = 0.1e1 + t197;
	t172 = t148 * t182;
	t171 = -0.2e1 * t148 * t206;
	t170 = t206 * t211;
	t169 = t130 * t140 * t198;
	t167 = t173 * t157;
	t166 = t155 * t199 - t200;
	t124 = 0.1e1 / t126;
	t123 = t173 * t156 * t140;
	t115 = 0.1e1 / t117;
	t113 = (t208 * t148 * t129 - t156 * t169) * t157;
	t111 = -t149 * t202 + t201 + (t129 * t149 - t130 * t189) * t123;
	t110 = -t173 * t170 + (qJD(1) * t167 + t165 * t211) * t140;
	t107 = -0.2e1 * t166 * t193 * t210 + (t166 * t149 * t190 + (-0.2e1 * t179 * t185 + t184 * t200 + (t128 * t186 + (t127 * t157 - t137 * t184) * t155) * t133) * t148) * t124;
	t106 = (t111 * t205 - t119 * t149) * t157 * t182 + ((-t119 * t184 + (-t111 * t151 - t109) * t204) * t149 + (-t119 * t190 - (-t110 * t130 * t156 - t209 * t129 + (t114 * t202 - t129 * t151 - t130 * t183) * t123) * t180 + (t120 * t184 + t157 * t181) * t111 - ((t110 - t183) * t129 + ((-t123 * t156 + 0.1e1) * t151 + (t123 - t156) * t114) * t130) * t149 * t204) * t148) * t115;
	t1 = [t157 * t145 * t171 + (t151 * t167 - t184 * t194) * t140, t110, t110, 0, 0, 0; (t119 * t172 + (-t119 * t192 + (qJD(1) * t113 + t109) * t205) * t115) * t156 + (t120 * t172 * t113 + (-((t171 - t192 + (t114 * t145 * t195 + t192) * t140) * t129 + (t170 * t198 - t114 * t148 + (-t143 * t177 + (t114 - 0.2e1 * t191) * t148) * t140) * t130) * t180 + (-t120 * t192 + t148 * t181) * t113 + (-t119 + ((-t152 + t153) * t169 + t208 * t178) * t120) * t148 * qJD(1)) * t115) * t157, t106, t106, 0, 0, 0; 0.2e1 * (-t132 * t135 + t136 * t199) * t210 + ((-t137 * qJD(1) + t154 * t176) * t132 + 0.2e1 * t136 * t179 + (-t135 * t128 - (-t138 * qJD(1) + t155 * t176) * t137 - t136 * t127) * t133) * t124, t107, t107, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:05
	% DurationCPUTime: 1.09s
	% Computational Cost: add. (6073->97), mult. (3810->205), div. (753->12), fcn. (4455->9), ass. (0->96)
	t182 = sin(qJ(1));
	t243 = 0.2e1 * t182;
	t180 = t182 ^ 2;
	t177 = qJ(2) + qJ(3) + pkin(10);
	t173 = sin(t177);
	t169 = t173 ^ 2;
	t174 = cos(t177);
	t171 = 0.1e1 / t174 ^ 2;
	t228 = t169 * t171;
	t164 = t180 * t228 + 0.1e1;
	t162 = 0.1e1 / t164;
	t170 = 0.1e1 / t174;
	t183 = cos(qJ(1));
	t214 = qJD(1) * t183;
	t204 = t173 * t214;
	t179 = qJD(2) + qJD(3);
	t222 = t179 * t182;
	t207 = t171 * t222;
	t136 = (-(-t174 * t222 - t204) * t170 + t169 * t207) * t162;
	t242 = t136 - t222;
	t178 = pkin(11) + qJ(6);
	t176 = cos(t178);
	t216 = t183 * t176;
	t175 = sin(t178);
	t219 = t182 * t175;
	t158 = t174 * t216 + t219;
	t220 = t182 * t173;
	t161 = atan2(-t220, -t174);
	t160 = cos(t161);
	t159 = sin(t161);
	t208 = t159 * t220;
	t146 = -t160 * t174 - t208;
	t143 = 0.1e1 / t146;
	t152 = 0.1e1 / t158;
	t144 = 0.1e1 / t146 ^ 2;
	t153 = 0.1e1 / t158 ^ 2;
	t241 = t162 - 0.1e1;
	t230 = t160 * t173;
	t131 = (-t136 * t182 + t179) * t230 + (t242 * t174 - t204) * t159;
	t240 = t131 * t143 * t144;
	t193 = t174 * t219 + t216;
	t221 = t179 * t183;
	t205 = t173 * t221;
	t137 = t193 * qJD(1) - t158 * qJD(6) + t175 * t205;
	t217 = t183 * t175;
	t218 = t182 * t176;
	t157 = t174 * t217 - t218;
	t151 = t157 ^ 2;
	t150 = t151 * t153 + 0.1e1;
	t233 = t153 * t157;
	t198 = -qJD(1) * t174 + qJD(6);
	t199 = qJD(6) * t174 - qJD(1);
	t138 = -t199 * t217 + (t198 * t182 - t205) * t176;
	t238 = t138 * t152 * t153;
	t239 = (-t137 * t233 - t151 * t238) / t150 ^ 2;
	t168 = t173 * t169;
	t225 = t170 * t173;
	t192 = t179 * (t168 * t170 * t171 + t225);
	t226 = t169 * t182;
	t196 = t214 * t226;
	t237 = (t171 * t196 + t180 * t192) / t164 ^ 2;
	t236 = t144 * t173;
	t235 = t144 * t183;
	t234 = t152 * t175;
	t232 = t157 * t176;
	t231 = t159 * t182;
	t229 = t169 * t170;
	t181 = t183 ^ 2;
	t227 = t169 * t181;
	t224 = t173 * t183;
	t223 = t174 * t179;
	t215 = qJD(1) * t182;
	t141 = t144 * t227 + 0.1e1;
	t213 = 0.2e1 * (-t227 * t240 + (t173 * t181 * t223 - t196) * t144) / t141 ^ 2;
	t212 = 0.2e1 * t240;
	t211 = -0.2e1 * t239;
	t210 = t144 * t224;
	t209 = t157 * t238;
	t203 = 0.1e1 + t228;
	t202 = t173 * t213;
	t201 = -0.2e1 * t173 * t237;
	t200 = t237 * t243;
	t197 = t160 * t162 * t229;
	t195 = t203 * t183;
	t194 = t153 * t232 - t234;
	t191 = t179 * t220 + t198 * t183;
	t156 = -t174 * t218 + t217;
	t148 = 0.1e1 / t150;
	t147 = t203 * t182 * t162;
	t139 = 0.1e1 / t141;
	t135 = (t241 * t173 * t159 - t182 * t197) * t183;
	t134 = -t174 * t231 + t230 + (t159 * t174 - t160 * t220) * t147;
	t132 = -t203 * t200 + (qJD(1) * t195 + t192 * t243) * t162;
	t129 = t194 * t211 * t224 + (t194 * t174 * t221 + (-t194 * t215 + ((-qJD(6) * t152 - 0.2e1 * t209) * t176 + (-t137 * t176 + (-qJD(6) * t157 + t138) * t175) * t153) * t183) * t173) * t148;
	t128 = (t134 * t236 - t143 * t174) * t183 * t213 + ((-t143 * t215 + (-t134 * t179 - t131) * t235) * t174 + (-t143 * t221 - (-t132 * t160 * t182 - t242 * t159 + (t136 * t231 - t159 * t179 - t160 * t214) * t147) * t210 + (t144 * t215 + t183 * t212) * t134 - ((t132 - t214) * t159 + ((-t147 * t182 + 0.1e1) * t179 + (t147 - t182) * t136) * t160) * t174 * t235) * t173) * t139;
	t1 = [t183 * t170 * t201 + (t179 * t195 - t215 * t225) * t162, t132, t132, 0, 0, 0; (t143 * t202 + (-t143 * t223 + (qJD(1) * t135 + t131) * t236) * t139) * t182 + (t144 * t202 * t135 + (-((t201 - t223 + (t136 * t170 * t226 + t223) * t162) * t159 + (t200 * t229 - t136 * t173 + (-t168 * t207 + (t136 - 0.2e1 * t222) * t173) * t162) * t160) * t210 + (-t144 * t223 + t173 * t212) * t135 + (-t143 + ((-t180 + t181) * t197 + t241 * t208) * t144) * t173 * qJD(1)) * t139) * t183, t128, t128, 0, 0, 0; 0.2e1 * (t152 * t193 + t156 * t233) * t239 + (0.2e1 * t156 * t209 - t199 * t152 * t218 + t191 * t234 + (-t199 * t157 * t219 + t156 * t137 + t138 * t193 - t191 * t232) * t153) * t148, t129, t129, 0, 0, t211 + 0.2e1 * (-t137 * t153 * t148 + (-t148 * t238 - t153 * t239) * t157) * t157;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end