% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP1
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
%   Wie in S6RRRRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:29
	% DurationCPUTime: 1.10s
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
	t124 = (-t129 * t173 + t168) * t222 + (t167 * t234 - t196) * t152;
	t232 = t124 * t136 * t137;
	t184 = t167 * t210 + t208;
	t213 = t168 * t175;
	t197 = t166 * t213;
	t133 = t184 * qJD(1) - qJD(4) * t151 + t172 * t197;
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
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:29
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
	t172 = qJ(4) + pkin(10);
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
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:30
	% DurationCPUTime: 1.72s
	% Computational Cost: add. (9575->126), mult. (8382->275), div. (1515->15), fcn. (10508->9), ass. (0->118)
	t202 = qJ(4) + pkin(10);
	t197 = sin(t202);
	t204 = qJ(2) + qJ(3);
	t200 = cos(t204);
	t198 = cos(t202);
	t205 = cos(qJ(1));
	t254 = t205 * t198;
	t276 = sin(qJ(1));
	t179 = t197 * t276 + t200 * t254;
	t173 = 0.1e1 / t179 ^ 2;
	t199 = sin(t204);
	t193 = t199 ^ 2;
	t203 = t205 ^ 2;
	t261 = t193 * t203;
	t241 = t173 * t261;
	t169 = 0.1e1 + t241;
	t229 = qJD(1) * t276;
	t201 = qJD(2) + qJD(3);
	t257 = t201 * t205;
	t235 = t199 * t257;
	t215 = t200 * t229 + t235;
	t228 = t276 * qJD(4);
	t255 = t205 * t197;
	t158 = (-qJD(4) * t200 + qJD(1)) * t255 + (t228 - t215) * t198;
	t172 = 0.1e1 / t179;
	t271 = t158 * t172 * t173;
	t223 = t261 * t271;
	t236 = t199 * t201 * t203;
	t279 = (-t223 + (-t193 * t205 * t229 + t200 * t236) * t173) / t169 ^ 2;
	t259 = t199 * t205;
	t233 = t276 * t200;
	t175 = t197 * t233 + t254;
	t220 = t197 * t228;
	t250 = qJD(4) * t205;
	t231 = t198 * t250;
	t157 = qJD(1) * t175 + t197 * t235 - t200 * t231 - t220;
	t178 = -t198 * t276 + t200 * t255;
	t190 = 0.1e1 / t197;
	t191 = 0.1e1 / t197 ^ 2;
	t194 = 0.1e1 / t199;
	t195 = 0.1e1 / t199 ^ 2;
	t258 = t200 * t201;
	t237 = t195 * t258;
	t252 = qJD(4) * t198;
	t264 = t190 * t194;
	t278 = (t191 * t194 * t252 + t190 * t237) * t178 + t157 * t264;
	t260 = t199 * t197;
	t165 = atan2(-t175, t260);
	t162 = cos(t165);
	t161 = sin(t165);
	t270 = t161 * t175;
	t156 = t162 * t260 - t270;
	t153 = 0.1e1 / t156;
	t154 = 0.1e1 / t156 ^ 2;
	t277 = 0.2e1 * t178;
	t170 = t175 ^ 2;
	t263 = t191 * t195;
	t166 = t170 * t263 + 0.1e1;
	t163 = 0.1e1 / t166;
	t251 = qJD(4) * t199;
	t216 = t197 * t258 + t198 * t251;
	t239 = t175 * t263;
	t221 = t198 * t229;
	t234 = t276 * t199;
	t222 = t201 * t234;
	t253 = qJD(1) * t205;
	t159 = t198 * t228 * t200 - t221 + (t200 * t253 - t222 - t250) * t197;
	t242 = t159 * t264;
	t145 = (t216 * t239 - t242) * t163;
	t213 = -t145 * t175 + t216;
	t141 = (-t145 * t260 - t159) * t161 + t213 * t162;
	t155 = t153 * t154;
	t275 = t141 * t155;
	t192 = t190 * t191;
	t196 = t194 / t193;
	t232 = t195 * t252;
	t274 = (t159 * t239 + (-t191 * t196 * t258 - t192 * t232) * t170) / t166 ^ 2;
	t273 = t154 * t178;
	t272 = t157 * t154;
	t269 = t161 * t178;
	t268 = t161 * t199;
	t267 = t162 * t175;
	t266 = t162 * t178;
	t265 = t162 * t200;
	t262 = t191 * t198;
	t256 = t205 * t153;
	t171 = t178 ^ 2;
	t151 = t154 * t171 + 0.1e1;
	t249 = 0.2e1 * (-t171 * t275 - t178 * t272) / t151 ^ 2;
	t248 = -0.2e1 * t274;
	t247 = 0.2e1 * t279;
	t246 = t155 * t277;
	t245 = t194 * t274;
	t244 = t154 * t269;
	t240 = t175 * t264;
	t238 = t190 * t195 * t200;
	t218 = t175 * t238 + t276;
	t152 = t218 * t163;
	t230 = t276 - t152;
	t227 = t153 * t249;
	t226 = t154 * t249;
	t225 = t259 * t277;
	t224 = t190 * t245;
	t177 = t198 * t233 - t255;
	t219 = t175 * t262 - t177 * t190;
	t217 = t173 * t177 * t205 - t172 * t276;
	t167 = 0.1e1 / t169;
	t160 = qJD(1) * t179 - t198 * t222 - t200 * t220 - t231;
	t149 = 0.1e1 / t151;
	t148 = t219 * t194 * t163;
	t144 = (-t161 + (t162 * t240 + t161) * t163) * t178;
	t143 = -t152 * t267 + (t230 * t268 + t265) * t197;
	t142 = t162 * t198 * t199 - t161 * t177 + (-t161 * t260 - t267) * t148;
	t140 = t218 * t248 + (t159 * t238 + t253 + (-t191 * t200 * t232 + (-0.2e1 * t196 * t200 ^ 2 - t194) * t201 * t190) * t175) * t163;
	t138 = (t172 * t200 * t205 + t198 * t241) * t247 + (0.2e1 * t198 * t223 + t215 * t172 + ((t158 * t205 - 0.2e1 * t198 * t236) * t200 + (qJD(4) * t197 * t203 + 0.2e1 * t205 * t221) * t193) * t173) * t167;
	t137 = -0.2e1 * t219 * t245 + (-t219 * t237 + (t159 * t262 - t160 * t190 + (t177 * t262 + (-0.2e1 * t192 * t198 ^ 2 - t190) * t175) * qJD(4)) * t194) * t163;
	t136 = t143 * t178 * t226 + (-(-t140 * t267 + (t145 * t270 - t159 * t162) * t152) * t273 + (t141 * t246 + t272) * t143 + (-t199 * t256 - (-t152 * t268 + t161 * t234 + t265) * t273) * t252) * t149 + (t227 * t259 + ((-t201 * t256 - (t201 * t230 - t145) * t244) * t200 + (t153 * t229 + (t205 * t141 - (-t140 + t253) * t269 - (t145 * t230 - t201) * t266) * t154) * t199) * t149) * t197;
	t1 = [t163 * t278 + t224 * t277, t140, t140, t137, 0, 0; t175 * t227 + (-t159 * t153 + (t141 * t175 + t144 * t157) * t154) * t149 + (t144 * t226 + (0.2e1 * t144 * t275 + (t157 * t163 - t157 - (-t145 * t163 * t240 + t248) * t178) * t154 * t161 + (-(-0.2e1 * t175 * t224 - t145) * t273 + (-(t145 + t242) * t178 + t278 * t175) * t154 * t163) * t162) * t149) * t178, t136, t136, (t142 * t273 - t153 * t179) * t249 + (t142 * t272 + t158 * t153 + (t142 * t246 - t154 * t179) * t141 - (-t197 * t251 + t198 * t258 - t137 * t175 - t148 * t159 + (-t148 * t260 - t177) * t145) * t154 * t266 - (-t160 + (-t137 * t197 - t145 * t198) * t199 - t213 * t148) * t244) * t149, 0, 0; t217 * t199 * t247 + (-t217 * t258 + ((qJD(1) * t172 + 0.2e1 * t177 * t271) * t205 + (-t158 * t276 - t160 * t205 + t177 * t229) * t173) * t199) * t167, t138, t138, t173 * t225 * t279 + (t225 * t271 + (t157 * t259 + (t199 * t229 - t200 * t257) * t178) * t173) * t167, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end