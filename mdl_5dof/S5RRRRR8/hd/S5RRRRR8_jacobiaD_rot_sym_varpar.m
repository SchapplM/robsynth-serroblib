% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR8
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
%   Wie in S5RRRRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:56
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRR8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:55
	% DurationCPUTime: 1.68s
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
	t1 = [t163 * t175 * t193 + (t168 * t187 - t207 * t217) * t155, t125, t125, 0, 0; (t136 * t194 + (-t136 * t215 + (qJD(1) * t128 + t124) * t228) * t130) * t173 + (t137 * t194 * t128 + (-((t193 - t215 + (t129 * t163 * t218 + t215) * t155) * t152 + (t192 * t221 - t129 * t166 + (-t161 * t199 + (t129 - 0.2e1 * t214) * t166) * t155) * t153) * t201 + (-t137 * t215 + t166 * t204) * t128 + (-t136 + ((-t169 + t170) * t189 + t233 * t200) * t137) * t166 * qJD(1)) * t130) * t175, t121, t121, 0, 0; 0.2e1 * (t145 * t184 + t149 * t225) * t231 + (0.2e1 * t149 * t202 - t191 * t145 * t209 + (t168 * t211 + t186) * t226 + (t149 * t133 + t184 * t134 - t186 * t224 - (t166 * t168 * t174 + t191 * t172) * t150 * t173) * t146) * t141, t122, t122, t203 + 0.2e1 * (-t133 * t141 * t146 + (-t141 * t230 - t146 * t231) * t150) * t150, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:55
	% DurationCPUTime: 1.73s
	% Computational Cost: add. (4423->98), mult. (4025->206), div. (771->12), fcn. (4686->9), ass. (0->97)
	t196 = sin(qJ(1));
	t256 = 0.2e1 * t196;
	t192 = t196 ^ 2;
	t195 = qJ(2) + qJ(3);
	t187 = sin(t195);
	t182 = t187 ^ 2;
	t189 = cos(t195);
	t184 = 0.1e1 / t189 ^ 2;
	t242 = t182 * t184;
	t178 = t192 * t242 + 0.1e1;
	t176 = 0.1e1 / t178;
	t183 = 0.1e1 / t189;
	t197 = cos(qJ(1));
	t228 = qJD(1) * t197;
	t218 = t187 * t228;
	t191 = qJD(2) + qJD(3);
	t234 = t191 * t196;
	t221 = t184 * t234;
	t149 = (-(-t189 * t234 - t218) * t183 + t182 * t221) * t176;
	t255 = t149 - t234;
	t194 = qJ(4) + qJ(5);
	t188 = cos(t194);
	t186 = sin(t194);
	t232 = t196 * t186;
	t235 = t189 * t197;
	t171 = t188 * t235 + t232;
	t231 = t196 * t187;
	t175 = atan2(-t231, -t189);
	t173 = cos(t175);
	t172 = sin(t175);
	t222 = t172 * t231;
	t159 = -t173 * t189 - t222;
	t156 = 0.1e1 / t159;
	t165 = 0.1e1 / t171;
	t157 = 0.1e1 / t159 ^ 2;
	t166 = 0.1e1 / t171 ^ 2;
	t254 = t176 - 0.1e1;
	t244 = t173 * t187;
	t144 = (-t149 * t196 + t191) * t244 + (t189 * t255 - t218) * t172;
	t253 = t144 * t156 * t157;
	t190 = qJD(4) + qJD(5);
	t207 = t188 * t197 + t189 * t232;
	t233 = t191 * t197;
	t219 = t187 * t233;
	t150 = t207 * qJD(1) - t171 * t190 + t186 * t219;
	t230 = t196 * t188;
	t170 = t186 * t235 - t230;
	t164 = t170 ^ 2;
	t163 = t164 * t166 + 0.1e1;
	t247 = t166 * t170;
	t212 = -qJD(1) * t189 + t190;
	t213 = t189 * t190 - qJD(1);
	t238 = t186 * t197;
	t151 = -t213 * t238 + (t212 * t196 - t219) * t188;
	t251 = t151 * t165 * t166;
	t252 = (-t150 * t247 - t164 * t251) / t163 ^ 2;
	t181 = t187 * t182;
	t239 = t183 * t187;
	t206 = t191 * (t181 * t183 * t184 + t239);
	t240 = t182 * t196;
	t210 = t228 * t240;
	t250 = (t184 * t210 + t192 * t206) / t178 ^ 2;
	t249 = t157 * t187;
	t248 = t165 * t186;
	t246 = t170 * t188;
	t245 = t172 * t196;
	t243 = t182 * t183;
	t193 = t197 ^ 2;
	t241 = t182 * t193;
	t237 = t187 * t197;
	t236 = t189 * t191;
	t229 = qJD(1) * t196;
	t154 = t157 * t241 + 0.1e1;
	t227 = 0.2e1 * (-t241 * t253 + (t187 * t193 * t236 - t210) * t157) / t154 ^ 2;
	t226 = 0.2e1 * t253;
	t225 = -0.2e1 * t252;
	t224 = t157 * t237;
	t223 = t170 * t251;
	t217 = 0.1e1 + t242;
	t216 = t187 * t227;
	t215 = -0.2e1 * t187 * t250;
	t214 = t250 * t256;
	t211 = t173 * t176 * t243;
	t209 = t217 * t197;
	t208 = t166 * t246 - t248;
	t205 = t191 * t231 + t212 * t197;
	t169 = -t189 * t230 + t238;
	t162 = t217 * t196 * t176;
	t160 = 0.1e1 / t163;
	t152 = 0.1e1 / t154;
	t148 = (t254 * t187 * t172 - t196 * t211) * t197;
	t147 = -t189 * t245 + t244 + (t172 * t189 - t173 * t231) * t162;
	t145 = -t217 * t214 + (qJD(1) * t209 + t206 * t256) * t176;
	t142 = t225 + 0.2e1 * (-t150 * t160 * t166 + (-t160 * t251 - t166 * t252) * t170) * t170;
	t141 = t208 * t225 * t237 + (t208 * t189 * t233 + (-t208 * t229 + ((-t165 * t190 - 0.2e1 * t223) * t188 + (-t150 * t188 + (-t170 * t190 + t151) * t186) * t166) * t197) * t187) * t160;
	t140 = (t147 * t249 - t156 * t189) * t197 * t227 + ((-t156 * t229 + (-t147 * t191 - t144) * t197 * t157) * t189 + (-t156 * t233 - (-t145 * t173 * t196 - t255 * t172 + (t149 * t245 - t172 * t191 - t173 * t228) * t162) * t224 + (t157 * t229 + t197 * t226) * t147 - ((t145 - t228) * t172 + ((-t162 * t196 + 0.1e1) * t191 + (t162 - t196) * t149) * t173) * t157 * t235) * t187) * t152;
	t1 = [t183 * t197 * t215 + (t191 * t209 - t229 * t239) * t176, t145, t145, 0, 0; (t156 * t216 + (-t156 * t236 + (qJD(1) * t148 + t144) * t249) * t152) * t196 + (t157 * t216 * t148 + (-((t215 - t236 + (t149 * t183 * t240 + t236) * t176) * t172 + (t214 * t243 - t149 * t187 + (-t181 * t221 + (t149 - 0.2e1 * t234) * t187) * t176) * t173) * t224 + (-t157 * t236 + t187 * t226) * t148 + (-t156 + ((-t192 + t193) * t211 + t254 * t222) * t157) * t187 * qJD(1)) * t152) * t197, t140, t140, 0, 0; 0.2e1 * (t165 * t207 + t169 * t247) * t252 + (0.2e1 * t169 * t223 - t213 * t165 * t230 + t205 * t248 + (-t213 * t170 * t232 + t169 * t150 + t151 * t207 - t205 * t246) * t166) * t160, t141, t141, t142, t142;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end