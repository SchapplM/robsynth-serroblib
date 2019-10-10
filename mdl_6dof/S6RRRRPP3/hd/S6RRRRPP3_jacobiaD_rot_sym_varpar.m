% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP3
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
%   Wie in S6RRRRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:01
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
	t124 = (-t129 * t173 + t168) * t222 + (t234 * t167 - t196) * t152;
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
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:02
	% DurationCPUTime: 1.50s
	% Computational Cost: add. (5428->122), mult. (8127->262), div. (1567->14), fcn. (10302->9), ass. (0->115)
	t200 = qJ(2) + qJ(3);
	t191 = sin(t200);
	t204 = cos(qJ(1));
	t277 = t191 * t204;
	t187 = 0.1e1 / t191;
	t188 = 0.1e1 / t191 ^ 2;
	t189 = t187 * t188;
	t192 = cos(t200);
	t193 = qJD(2) + qJD(3);
	t276 = t193 * (0.2e1 * t189 * t192 ^ 2 + t187);
	t203 = cos(qJ(4));
	t250 = t204 * t203;
	t201 = sin(qJ(4));
	t202 = sin(qJ(1));
	t254 = t202 * t201;
	t172 = t192 * t254 + t250;
	t244 = qJD(4) * t204;
	t224 = t203 * t244;
	t246 = qJD(4) * t201;
	t225 = t202 * t246;
	t228 = t193 * t277;
	t156 = t172 * qJD(1) - t192 * t224 + t201 * t228 - t225;
	t251 = t204 * t201;
	t253 = t202 * t203;
	t175 = t192 * t251 - t253;
	t195 = 0.1e1 / t201 ^ 2;
	t245 = qJD(4) * t203;
	t226 = t195 * t245;
	t194 = 0.1e1 / t201;
	t262 = t188 * t192;
	t232 = t194 * t262;
	t263 = t187 * t194;
	t275 = t175 * (t187 * t226 + t193 * t232) + t156 * t263;
	t259 = t191 * t201;
	t166 = atan2(-t172, t259);
	t161 = cos(t166);
	t160 = sin(t166);
	t269 = t160 * t172;
	t155 = t161 * t259 - t269;
	t152 = 0.1e1 / t155;
	t197 = 0.1e1 / t204;
	t153 = 0.1e1 / t155 ^ 2;
	t198 = 0.1e1 / t204 ^ 2;
	t274 = 0.2e1 * t175;
	t169 = t172 ^ 2;
	t261 = t188 * t195;
	t167 = t169 * t261 + 0.1e1;
	t162 = 0.1e1 / t167;
	t257 = t192 * t193;
	t215 = t191 * t245 + t201 * t257;
	t234 = t172 * t261;
	t258 = t191 * t202;
	t229 = t193 * t258;
	t247 = qJD(1) * t204;
	t248 = qJD(1) * t202;
	t158 = t202 * t245 * t192 - t203 * t248 + (t247 * t192 - t229 - t244) * t201;
	t236 = t158 * t263;
	t144 = (t215 * t234 - t236) * t162;
	t213 = -t144 * t172 + t215;
	t140 = (-t144 * t259 - t158) * t160 + t213 * t161;
	t154 = t152 * t153;
	t273 = t140 * t154;
	t196 = t194 * t195;
	t230 = t189 * t257;
	t272 = (t158 * t234 + (-t188 * t196 * t245 - t195 * t230) * t169) / t167 ^ 2;
	t271 = t153 * t175;
	t270 = t156 * t153;
	t268 = t160 * t175;
	t267 = t160 * t191;
	t266 = t161 * t172;
	t265 = t161 * t175;
	t264 = t161 * t192;
	t260 = t188 * t198;
	t256 = t195 * t203;
	t255 = t198 * t202;
	t252 = t204 * t152;
	t218 = t172 * t232 + t202;
	t151 = t218 * t162;
	t249 = -t151 + t202;
	t170 = t175 ^ 2;
	t150 = t170 * t153 + 0.1e1;
	t243 = 0.2e1 * (-t170 * t273 - t175 * t270) / t150 ^ 2;
	t242 = -0.2e1 * t272;
	t157 = (-qJD(4) * t192 + qJD(1)) * t251 + (-t228 + (-qJD(1) * t192 + qJD(4)) * t202) * t203;
	t176 = t192 * t250 + t254;
	t171 = t176 ^ 2;
	t168 = t171 * t260 + 0.1e1;
	t199 = t197 * t198;
	t241 = 0.2e1 * (t176 * t157 * t260 + (t188 * t199 * t248 - t198 * t230) * t171) / t168 ^ 2;
	t240 = t154 * t274;
	t239 = t187 * t272;
	t238 = t153 * t268;
	t235 = t172 * t263;
	t233 = t188 * t257;
	t231 = t197 * t262;
	t227 = t198 * t248;
	t223 = t152 * t243;
	t222 = t153 * t243;
	t221 = t187 * t241;
	t219 = t194 * t239;
	t174 = t192 * t253 - t251;
	t217 = t172 * t256 - t174 * t194;
	t216 = t174 * t197 - t176 * t255;
	t164 = 0.1e1 / t168;
	t159 = qJD(1) * t176 - t192 * t225 - t203 * t229 - t224;
	t148 = 0.1e1 / t150;
	t147 = t217 * t187 * t162;
	t143 = (-t160 + (t161 * t235 + t160) * t162) * t175;
	t142 = -t151 * t266 + (t249 * t267 + t264) * t201;
	t141 = t161 * t191 * t203 - t160 * t174 + (-t160 * t259 - t266) * t147;
	t139 = t218 * t242 + (t158 * t232 + t247 + (-t194 * t276 - t226 * t262) * t172) * t162;
	t138 = (t176 * t231 + t203) * t241 + (-t157 * t231 + t246 + (t197 * t276 - t227 * t262) * t176) * t164;
	t136 = -0.2e1 * t217 * t239 + (-t217 * t233 + (t158 * t256 - t159 * t194 + (t174 * t256 + (-0.2e1 * t196 * t203 ^ 2 - t194) * t172) * qJD(4)) * t187) * t162;
	t135 = t142 * t175 * t222 + (-(-t139 * t266 + (t144 * t269 - t158 * t161) * t151) * t271 + (t140 * t240 + t270) * t142 + (-t191 * t252 - (-t151 * t267 + t160 * t258 + t264) * t271) * t245) * t148 + (t223 * t277 + ((-t193 * t252 - (t193 * t249 - t144) * t238) * t192 + (t152 * t248 + (t204 * t140 - (-t139 + t247) * t268 - (t144 * t249 - t193) * t265) * t153) * t191) * t148) * t201;
	t1 = [t275 * t162 + t219 * t274, t139, t139, t136, 0, 0; t172 * t223 + (-t158 * t152 + (t140 * t172 + t143 * t156) * t153) * t148 + (t143 * t222 + (0.2e1 * t143 * t273 + (t156 * t162 - t156 - (-t144 * t162 * t235 + t242) * t175) * t153 * t160 + (-(-0.2e1 * t172 * t219 - t144) * t271 + (-(t144 + t236) * t175 + t275 * t172) * t153 * t162) * t161) * t148) * t175, t135, t135, (t141 * t271 - t152 * t176) * t243 + (t141 * t270 + t157 * t152 + (t141 * t240 - t176 * t153) * t140 - (-t191 * t246 + t203 * t257 - t136 * t172 - t147 * t158 + (-t147 * t259 - t174) * t144) * t153 * t265 - (-t159 + (-t136 * t201 - t144 * t203) * t191 - t213 * t147) * t238) * t148, 0, 0; t216 * t221 + (t216 * t233 + (t157 * t255 - t159 * t197 + (-t174 * t255 + (0.2e1 * t199 * t202 ^ 2 + t197) * t176) * qJD(1)) * t187) * t164, t138, t138, t175 * t197 * t221 + (t156 * t187 * t197 + (-t187 * t227 + t193 * t231) * t175) * t164, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:02
	% DurationCPUTime: 1.53s
	% Computational Cost: add. (5428->121), mult. (8127->261), div. (1567->14), fcn. (10302->9), ass. (0->116)
	t197 = qJ(2) + qJ(3);
	t188 = sin(t197);
	t201 = cos(qJ(1));
	t275 = t188 * t201;
	t184 = 0.1e1 / t188;
	t185 = 0.1e1 / t188 ^ 2;
	t186 = t184 * t185;
	t189 = cos(t197);
	t190 = qJD(2) + qJD(3);
	t274 = t190 * (0.2e1 * t186 * t189 ^ 2 + t184);
	t198 = sin(qJ(4));
	t199 = sin(qJ(1));
	t200 = cos(qJ(4));
	t226 = t190 * t275;
	t242 = qJD(4) * t201;
	t243 = qJD(4) * t200;
	t245 = qJD(1) * t201;
	t246 = qJD(1) * t199;
	t158 = t200 * t226 - t199 * t243 - t198 * t245 + (t198 * t242 + t200 * t246) * t189;
	t250 = t200 * t201;
	t252 = t199 * t198;
	t177 = t189 * t250 + t252;
	t192 = 0.1e1 / t200 ^ 2;
	t244 = qJD(4) * t198;
	t223 = t192 * t244;
	t191 = 0.1e1 / t200;
	t260 = t185 * t189;
	t230 = t191 * t260;
	t261 = t184 * t191;
	t273 = -(t184 * t223 - t190 * t230) * t177 + t158 * t261;
	t248 = t201 * t198;
	t251 = t199 * t200;
	t174 = t189 * t251 - t248;
	t256 = t188 * t200;
	t167 = atan2(-t174, t256);
	t162 = cos(t167);
	t161 = sin(t167);
	t267 = t161 * t174;
	t156 = t162 * t256 - t267;
	t153 = 0.1e1 / t156;
	t194 = 0.1e1 / t201;
	t154 = 0.1e1 / t156 ^ 2;
	t195 = 0.1e1 / t201 ^ 2;
	t272 = 0.2e1 * t177;
	t170 = t174 ^ 2;
	t259 = t185 * t192;
	t168 = t170 * t259 + 0.1e1;
	t163 = 0.1e1 / t168;
	t255 = t189 * t190;
	t212 = -t188 * t244 + t200 * t255;
	t232 = t174 * t259;
	t225 = t189 * t252;
	t257 = t188 * t199;
	t227 = t190 * t257;
	t160 = t177 * qJD(1) - qJD(4) * t225 + (-t227 - t242) * t200;
	t234 = t160 * t261;
	t145 = (t212 * t232 - t234) * t163;
	t210 = t145 * t174 - t212;
	t141 = (-t145 * t256 - t160) * t161 - t210 * t162;
	t155 = t153 * t154;
	t271 = t141 * t155;
	t193 = t191 * t192;
	t228 = t186 * t255;
	t270 = (t160 * t232 + (t185 * t193 * t244 - t192 * t228) * t170) / t168 ^ 2;
	t269 = t154 * t177;
	t268 = t158 * t154;
	t266 = t161 * t177;
	t265 = t161 * t188;
	t264 = t162 * t174;
	t263 = t162 * t177;
	t262 = t162 * t189;
	t258 = t185 * t195;
	t254 = t192 * t198;
	t253 = t195 * t199;
	t249 = t201 * t153;
	t215 = t174 * t230 + t199;
	t152 = t215 * t163;
	t247 = -t152 + t199;
	t172 = t177 ^ 2;
	t151 = t154 * t172 + 0.1e1;
	t241 = 0.2e1 * (-t172 * t271 - t177 * t268) / t151 ^ 2;
	t240 = -0.2e1 * t270;
	t217 = qJD(1) * t189 - qJD(4);
	t218 = qJD(4) * t189 - qJD(1);
	t157 = -t218 * t250 + (t199 * t217 + t226) * t198;
	t176 = -t189 * t248 + t251;
	t171 = t176 ^ 2;
	t169 = t171 * t258 + 0.1e1;
	t196 = t194 * t195;
	t239 = 0.2e1 * (t157 * t176 * t258 + (t185 * t196 * t246 - t195 * t228) * t171) / t169 ^ 2;
	t238 = t155 * t272;
	t237 = t184 * t270;
	t236 = t154 * t266;
	t233 = t174 * t261;
	t231 = t185 * t255;
	t229 = t194 * t260;
	t224 = t195 * t246;
	t222 = t153 * t241;
	t221 = t154 * t241;
	t220 = t184 * t239;
	t216 = t191 * t237;
	t173 = t225 + t250;
	t214 = -t173 * t191 + t174 * t254;
	t213 = -t173 * t194 - t176 * t253;
	t165 = 0.1e1 / t169;
	t159 = t218 * t251 + (t201 * t217 - t227) * t198;
	t149 = 0.1e1 / t151;
	t148 = t214 * t184 * t163;
	t144 = (-t161 + (t162 * t233 + t161) * t163) * t177;
	t143 = -t152 * t264 + (t247 * t265 + t262) * t200;
	t142 = -t162 * t188 * t198 + t161 * t173 - (-t161 * t256 - t264) * t148;
	t140 = t215 * t240 + (t160 * t230 + t245 + (-t191 * t274 + t223 * t260) * t174) * t163;
	t139 = (t176 * t229 - t198) * t239 + (-t157 * t229 + t243 + (t194 * t274 - t224 * t260) * t176) * t165;
	t137 = 0.2e1 * t214 * t237 + (t214 * t231 + (-t160 * t254 + t159 * t191 + (t173 * t254 + (-0.2e1 * t193 * t198 ^ 2 - t191) * t174) * qJD(4)) * t184) * t163;
	t136 = t143 * t177 * t221 + (-(-t140 * t264 + (t145 * t267 - t160 * t162) * t152) * t269 + (t141 * t238 + t268) * t143 + (t188 * t249 - (t152 * t265 - t161 * t257 - t262) * t269) * t244) * t149 + (t222 * t275 + ((-t190 * t249 - (t190 * t247 - t145) * t236) * t189 + (t153 * t246 + (t201 * t141 - (-t140 + t245) * t266 - (t145 * t247 - t190) * t263) * t154) * t188) * t149) * t200;
	t1 = [t163 * t273 + t216 * t272, t140, t140, t137, 0, 0; t174 * t222 + (-t160 * t153 + (t141 * t174 + t144 * t158) * t154) * t149 + (t144 * t221 + (0.2e1 * t144 * t271 + (t158 * t163 - t158 - (-t145 * t163 * t233 + t240) * t177) * t154 * t161 + (-(-0.2e1 * t174 * t216 - t145) * t269 + (-(t145 + t234) * t177 + t273 * t174) * t154 * t163) * t162) * t149) * t177, t136, t136, (t142 * t269 - t153 * t176) * t241 + (t142 * t268 + t157 * t153 + (t142 * t238 - t154 * t176) * t141 - (-t188 * t243 - t198 * t255 - t137 * t174 + t148 * t160 + (t148 * t256 + t173) * t145) * t154 * t263 - (t159 + (-t137 * t200 + t145 * t198) * t188 - t210 * t148) * t236) * t149, 0, 0; t213 * t220 + (t213 * t231 + (t157 * t253 + t159 * t194 + (t173 * t253 + (0.2e1 * t196 * t199 ^ 2 + t194) * t176) * qJD(1)) * t184) * t165, t139, t139, t177 * t194 * t220 + (t158 * t184 * t194 + (-t184 * t224 + t190 * t229) * t177) * t165, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end