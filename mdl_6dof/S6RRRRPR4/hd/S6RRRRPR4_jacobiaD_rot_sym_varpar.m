% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR4
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
%   Wie in S6RRRRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:36
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:52
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (3645->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t173 = sin(qJ(1));
	t235 = 0.2e1 * t173;
	t169 = t173 ^ 2;
	t171 = qJ(2) + qJ(3);
	t166 = sin(t171);
	t162 = t166 ^ 2;
	t167 = cos(t171);
	t164 = 0.1e1 / t167 ^ 2;
	t220 = t162 * t164;
	t157 = t169 * t220 + 0.1e1;
	t155 = 0.1e1 / t157;
	t163 = 0.1e1 / t167;
	t175 = cos(qJ(1));
	t206 = qJD(1) * t175;
	t196 = t166 * t206;
	t168 = qJD(2) + qJD(3);
	t214 = t168 * t173;
	t199 = t164 * t214;
	t129 = (-(-t167 * t214 - t196) * t163 + t162 * t199) * t155;
	t234 = t129 - t214;
	t174 = cos(qJ(4));
	t208 = t174 * t175;
	t172 = sin(qJ(4));
	t210 = t173 * t172;
	t151 = t167 * t208 + t210;
	t211 = t173 * t166;
	t154 = atan2(-t211, -t167);
	t153 = cos(t154);
	t152 = sin(t154);
	t200 = t152 * t211;
	t139 = -t153 * t167 - t200;
	t136 = 0.1e1 / t139;
	t145 = 0.1e1 / t151;
	t137 = 0.1e1 / t139 ^ 2;
	t146 = 0.1e1 / t151 ^ 2;
	t233 = t155 - 0.1e1;
	t222 = t153 * t166;
	t124 = (-t129 * t173 + t168) * t222 + (t234 * t167 - t196) * t152;
	t232 = t124 * t136 * t137;
	t184 = t167 * t210 + t208;
	t213 = t168 * t175;
	t197 = t166 * t213;
	t133 = t184 * qJD(1) - t151 * qJD(4) + t172 * t197;
	t209 = t173 * t174;
	t212 = t172 * t175;
	t150 = t167 * t212 - t209;
	t144 = t150 ^ 2;
	t143 = t144 * t146 + 0.1e1;
	t225 = t146 * t150;
	t190 = -qJD(1) * t167 + qJD(4);
	t191 = qJD(4) * t167 - qJD(1);
	t134 = -t191 * t212 + (t190 * t173 - t197) * t174;
	t230 = t134 * t145 * t146;
	t231 = (-t133 * t225 - t144 * t230) / t143 ^ 2;
	t161 = t166 * t162;
	t217 = t163 * t166;
	t183 = t168 * (t161 * t163 * t164 + t217);
	t218 = t162 * t173;
	t188 = t206 * t218;
	t229 = (t164 * t188 + t169 * t183) / t157 ^ 2;
	t228 = t137 * t166;
	t227 = t137 * t175;
	t226 = t145 * t172;
	t224 = t150 * t174;
	t223 = t152 * t173;
	t221 = t162 * t163;
	t170 = t175 ^ 2;
	t219 = t162 * t170;
	t216 = t166 * t175;
	t215 = t167 * t168;
	t207 = qJD(1) * t173;
	t132 = t137 * t219 + 0.1e1;
	t205 = 0.2e1 * (-t219 * t232 + (t166 * t170 * t215 - t188) * t137) / t132 ^ 2;
	t204 = 0.2e1 * t232;
	t203 = -0.2e1 * t231;
	t202 = t150 * t230;
	t201 = t137 * t216;
	t195 = 0.1e1 + t220;
	t194 = t166 * t205;
	t193 = -0.2e1 * t166 * t229;
	t192 = t229 * t235;
	t189 = t153 * t155 * t221;
	t187 = t195 * t175;
	t186 = t190 * t175;
	t185 = t146 * t224 - t226;
	t149 = -t167 * t209 + t212;
	t141 = 0.1e1 / t143;
	t140 = t195 * t173 * t155;
	t130 = 0.1e1 / t132;
	t128 = (t233 * t166 * t152 - t173 * t189) * t175;
	t126 = -t167 * t223 + t222 + (t152 * t167 - t153 * t211) * t140;
	t125 = -t195 * t192 + (qJD(1) * t187 + t183 * t235) * t155;
	t122 = t185 * t203 * t216 + (t185 * t167 * t213 + (-t185 * t207 + ((-qJD(4) * t145 - 0.2e1 * t202) * t174 + (-t133 * t174 + (-qJD(4) * t150 + t134) * t172) * t146) * t175) * t166) * t141;
	t121 = (t126 * t228 - t136 * t167) * t175 * t205 + ((-t136 * t207 + (-t126 * t168 - t124) * t227) * t167 + (-t136 * t213 - (-t125 * t153 * t173 - t234 * t152 + (t129 * t223 - t152 * t168 - t153 * t206) * t140) * t201 + (t137 * t207 + t175 * t204) * t126 - ((t125 - t206) * t152 + ((-t140 * t173 + 0.1e1) * t168 + (t140 - t173) * t129) * t153) * t167 * t227) * t166) * t130;
	t1 = [t163 * t175 * t193 + (t168 * t187 - t207 * t217) * t155, t125, t125, 0, 0, 0; (t136 * t194 + (-t136 * t215 + (qJD(1) * t128 + t124) * t228) * t130) * t173 + (t137 * t194 * t128 + (-((t193 - t215 + (t129 * t163 * t218 + t215) * t155) * t152 + (t192 * t221 - t129 * t166 + (-t161 * t199 + (t129 - 0.2e1 * t214) * t166) * t155) * t153) * t201 + (-t137 * t215 + t166 * t204) * t128 + (-t136 + ((-t169 + t170) * t189 + t233 * t200) * t137) * t166 * qJD(1)) * t130) * t175, t121, t121, 0, 0, 0; 0.2e1 * (t145 * t184 + t149 * t225) * t231 + (0.2e1 * t149 * t202 - t191 * t145 * t209 + (t168 * t211 + t186) * t226 + (t149 * t133 + t184 * t134 - t186 * t224 - (t166 * t168 * t174 + t191 * t172) * t150 * t173) * t146) * t141, t122, t122, t203 + 0.2e1 * (-t133 * t141 * t146 + (-t141 * t230 - t146 * t231) * t150) * t150, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:52
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (4131->97), mult. (3810->205), div. (753->12), fcn. (4455->9), ass. (0->96)
	t176 = sin(qJ(1));
	t237 = 0.2e1 * t176;
	t173 = t176 ^ 2;
	t175 = qJ(2) + qJ(3);
	t169 = sin(t175);
	t163 = t169 ^ 2;
	t170 = cos(t175);
	t165 = 0.1e1 / t170 ^ 2;
	t222 = t163 * t165;
	t159 = t173 * t222 + 0.1e1;
	t156 = 0.1e1 / t159;
	t164 = 0.1e1 / t170;
	t177 = cos(qJ(1));
	t208 = qJD(1) * t177;
	t198 = t169 * t208;
	t171 = qJD(2) + qJD(3);
	t216 = t171 * t176;
	t201 = t165 * t216;
	t130 = (-(-t170 * t216 - t198) * t164 + t163 * t201) * t156;
	t236 = t130 - t216;
	t172 = qJ(4) + pkin(11);
	t168 = cos(t172);
	t210 = t177 * t168;
	t167 = sin(t172);
	t214 = t176 * t167;
	t152 = t170 * t210 + t214;
	t212 = t176 * t169;
	t155 = atan2(-t212, -t170);
	t154 = cos(t155);
	t153 = sin(t155);
	t202 = t153 * t212;
	t140 = -t154 * t170 - t202;
	t137 = 0.1e1 / t140;
	t146 = 0.1e1 / t152;
	t138 = 0.1e1 / t140 ^ 2;
	t147 = 0.1e1 / t152 ^ 2;
	t235 = t156 - 0.1e1;
	t224 = t154 * t169;
	t125 = (-t130 * t176 + t171) * t224 + (t236 * t170 - t198) * t153;
	t234 = t125 * t137 * t138;
	t187 = t170 * t214 + t210;
	t215 = t171 * t177;
	t199 = t169 * t215;
	t131 = t187 * qJD(1) - t152 * qJD(4) + t167 * t199;
	t211 = t177 * t167;
	t213 = t176 * t168;
	t151 = t170 * t211 - t213;
	t145 = t151 ^ 2;
	t143 = t145 * t147 + 0.1e1;
	t227 = t147 * t151;
	t192 = -qJD(1) * t170 + qJD(4);
	t193 = qJD(4) * t170 - qJD(1);
	t132 = -t193 * t211 + (t192 * t176 - t199) * t168;
	t232 = t132 * t146 * t147;
	t233 = (-t131 * t227 - t145 * t232) / t143 ^ 2;
	t162 = t169 * t163;
	t219 = t164 * t169;
	t186 = t171 * (t162 * t164 * t165 + t219);
	t220 = t163 * t176;
	t190 = t208 * t220;
	t231 = (t165 * t190 + t173 * t186) / t159 ^ 2;
	t230 = t138 * t169;
	t229 = t138 * t177;
	t228 = t146 * t167;
	t226 = t151 * t168;
	t225 = t153 * t176;
	t223 = t163 * t164;
	t174 = t177 ^ 2;
	t221 = t163 * t174;
	t218 = t169 * t177;
	t217 = t170 * t171;
	t209 = qJD(1) * t176;
	t135 = t138 * t221 + 0.1e1;
	t207 = 0.2e1 * (-t221 * t234 + (t169 * t174 * t217 - t190) * t138) / t135 ^ 2;
	t206 = 0.2e1 * t234;
	t205 = -0.2e1 * t233;
	t204 = t138 * t218;
	t203 = t151 * t232;
	t197 = 0.1e1 + t222;
	t196 = t169 * t207;
	t195 = -0.2e1 * t169 * t231;
	t194 = t231 * t237;
	t191 = t154 * t156 * t223;
	t189 = t197 * t177;
	t188 = t147 * t226 - t228;
	t185 = t171 * t212 + t192 * t177;
	t150 = -t170 * t213 + t211;
	t144 = t197 * t176 * t156;
	t141 = 0.1e1 / t143;
	t133 = 0.1e1 / t135;
	t129 = (t235 * t169 * t153 - t176 * t191) * t177;
	t128 = -t170 * t225 + t224 + (t153 * t170 - t154 * t212) * t144;
	t126 = -t197 * t194 + (qJD(1) * t189 + t186 * t237) * t156;
	t123 = t188 * t205 * t218 + (t188 * t170 * t215 + (-t188 * t209 + ((-qJD(4) * t146 - 0.2e1 * t203) * t168 + (-t131 * t168 + (-qJD(4) * t151 + t132) * t167) * t147) * t177) * t169) * t141;
	t122 = (t128 * t230 - t137 * t170) * t177 * t207 + ((-t137 * t209 + (-t128 * t171 - t125) * t229) * t170 + (-t137 * t215 - (-t126 * t154 * t176 - t236 * t153 + (t130 * t225 - t153 * t171 - t154 * t208) * t144) * t204 + (t138 * t209 + t177 * t206) * t128 - ((t126 - t208) * t153 + ((-t144 * t176 + 0.1e1) * t171 + (t144 - t176) * t130) * t154) * t170 * t229) * t169) * t133;
	t1 = [t177 * t164 * t195 + (t171 * t189 - t209 * t219) * t156, t126, t126, 0, 0, 0; (t137 * t196 + (-t137 * t217 + (qJD(1) * t129 + t125) * t230) * t133) * t176 + (t138 * t196 * t129 + (-((t195 - t217 + (t130 * t164 * t220 + t217) * t156) * t153 + (t194 * t223 - t130 * t169 + (-t162 * t201 + (t130 - 0.2e1 * t216) * t169) * t156) * t154) * t204 + (-t138 * t217 + t169 * t206) * t129 + (-t137 + ((-t173 + t174) * t191 + t235 * t202) * t138) * t169 * qJD(1)) * t133) * t177, t122, t122, 0, 0, 0; 0.2e1 * (t146 * t187 + t150 * t227) * t233 + (0.2e1 * t150 * t203 - t193 * t146 * t213 + t185 * t228 + (-t151 * t193 * t214 + t150 * t131 + t132 * t187 - t185 * t226) * t147) * t141, t123, t123, t205 + 0.2e1 * (-t131 * t147 * t141 + (-t141 * t232 - t147 * t233) * t151) * t151, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:52
	% DurationCPUTime: 1.10s
	% Computational Cost: add. (5000->98), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->98)
	t196 = sin(qJ(1));
	t257 = 0.2e1 * t196;
	t193 = t196 ^ 2;
	t195 = qJ(2) + qJ(3);
	t189 = sin(t195);
	t184 = t189 ^ 2;
	t190 = cos(t195);
	t186 = 0.1e1 / t190 ^ 2;
	t242 = t184 * t186;
	t178 = t193 * t242 + 0.1e1;
	t176 = 0.1e1 / t178;
	t185 = 0.1e1 / t190;
	t197 = cos(qJ(1));
	t228 = qJD(1) * t197;
	t218 = t189 * t228;
	t192 = qJD(2) + qJD(3);
	t236 = t192 * t196;
	t221 = t186 * t236;
	t151 = (-(-t190 * t236 - t218) * t185 + t184 * t221) * t176;
	t256 = t151 - t236;
	t188 = qJ(4) + pkin(11) + qJ(6);
	t182 = cos(t188);
	t230 = t197 * t182;
	t181 = sin(t188);
	t234 = t196 * t181;
	t171 = t190 * t230 + t234;
	t232 = t196 * t189;
	t175 = atan2(-t232, -t190);
	t174 = cos(t175);
	t173 = sin(t175);
	t222 = t173 * t232;
	t162 = -t174 * t190 - t222;
	t159 = 0.1e1 / t162;
	t165 = 0.1e1 / t171;
	t160 = 0.1e1 / t162 ^ 2;
	t166 = 0.1e1 / t171 ^ 2;
	t255 = t176 - 0.1e1;
	t244 = t174 * t189;
	t144 = (-t151 * t196 + t192) * t244 + (t256 * t190 - t218) * t173;
	t254 = t144 * t159 * t160;
	t191 = qJD(4) + qJD(6);
	t207 = t190 * t234 + t230;
	t235 = t192 * t197;
	t219 = t189 * t235;
	t149 = t207 * qJD(1) - t171 * t191 + t181 * t219;
	t231 = t197 * t181;
	t233 = t196 * t182;
	t170 = t190 * t231 - t233;
	t164 = t170 ^ 2;
	t157 = t164 * t166 + 0.1e1;
	t247 = t166 * t170;
	t212 = -qJD(1) * t190 + t191;
	t213 = t190 * t191 - qJD(1);
	t150 = -t213 * t231 + (t212 * t196 - t219) * t182;
	t252 = t150 * t165 * t166;
	t253 = (-t149 * t247 - t164 * t252) / t157 ^ 2;
	t183 = t189 * t184;
	t239 = t185 * t189;
	t206 = t192 * (t183 * t185 * t186 + t239);
	t240 = t184 * t196;
	t210 = t228 * t240;
	t251 = (t186 * t210 + t193 * t206) / t178 ^ 2;
	t250 = t160 * t189;
	t249 = t160 * t197;
	t248 = t165 * t181;
	t246 = t170 * t182;
	t245 = t173 * t196;
	t243 = t184 * t185;
	t194 = t197 ^ 2;
	t241 = t184 * t194;
	t238 = t189 * t197;
	t237 = t190 * t192;
	t229 = qJD(1) * t196;
	t154 = t160 * t241 + 0.1e1;
	t227 = 0.2e1 * (-t241 * t254 + (t189 * t194 * t237 - t210) * t160) / t154 ^ 2;
	t226 = 0.2e1 * t254;
	t225 = -0.2e1 * t253;
	t224 = t160 * t238;
	t223 = t170 * t252;
	t217 = 0.1e1 + t242;
	t216 = t189 * t227;
	t215 = -0.2e1 * t189 * t251;
	t214 = t251 * t257;
	t211 = t174 * t176 * t243;
	t209 = t217 * t197;
	t208 = t166 * t246 - t248;
	t205 = t192 * t232 + t212 * t197;
	t169 = -t190 * t233 + t231;
	t163 = t217 * t196 * t176;
	t155 = 0.1e1 / t157;
	t152 = 0.1e1 / t154;
	t148 = (t255 * t189 * t173 - t196 * t211) * t197;
	t147 = -t190 * t245 + t244 + (t173 * t190 - t174 * t232) * t163;
	t145 = -t217 * t214 + (qJD(1) * t209 + t206 * t257) * t176;
	t142 = t225 + 0.2e1 * (-t149 * t166 * t155 + (-t155 * t252 - t166 * t253) * t170) * t170;
	t141 = t208 * t225 * t238 + (t208 * t190 * t235 + (-t208 * t229 + ((-t165 * t191 - 0.2e1 * t223) * t182 + (-t149 * t182 + (-t170 * t191 + t150) * t181) * t166) * t197) * t189) * t155;
	t140 = (t147 * t250 - t159 * t190) * t197 * t227 + ((-t159 * t229 + (-t147 * t192 - t144) * t249) * t190 + (-t159 * t235 - (-t145 * t174 * t196 - t256 * t173 + (t151 * t245 - t173 * t192 - t174 * t228) * t163) * t224 + (t160 * t229 + t197 * t226) * t147 - ((t145 - t228) * t173 + ((-t163 * t196 + 0.1e1) * t192 + (t163 - t196) * t151) * t174) * t190 * t249) * t189) * t152;
	t1 = [t197 * t185 * t215 + (t192 * t209 - t229 * t239) * t176, t145, t145, 0, 0, 0; (t159 * t216 + (-t159 * t237 + (qJD(1) * t148 + t144) * t250) * t152) * t196 + (t160 * t216 * t148 + (-((t215 - t237 + (t151 * t185 * t240 + t237) * t176) * t173 + (t214 * t243 - t151 * t189 + (-t183 * t221 + (t151 - 0.2e1 * t236) * t189) * t176) * t174) * t224 + (-t160 * t237 + t189 * t226) * t148 + (-t159 + ((-t193 + t194) * t211 + t255 * t222) * t160) * t189 * qJD(1)) * t152) * t197, t140, t140, 0, 0, 0; 0.2e1 * (t165 * t207 + t169 * t247) * t253 + (0.2e1 * t169 * t223 - t213 * t165 * t233 + t205 * t248 + (-t213 * t170 * t234 + t169 * t149 + t150 * t207 - t205 * t246) * t166) * t155, t141, t141, t142, 0, t142;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end