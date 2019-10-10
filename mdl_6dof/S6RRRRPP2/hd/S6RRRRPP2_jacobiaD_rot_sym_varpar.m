% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP2
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
%   Wie in S6RRRRPP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:16
	% EndTime: 2019-10-10 12:22:17
	% DurationCPUTime: 1.43s
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
	% StartTime: 2019-10-10 12:22:16
	% EndTime: 2019-10-10 12:22:17
	% DurationCPUTime: 1.60s
	% Computational Cost: add. (5529->125), mult. (8382->273), div. (1515->15), fcn. (10508->9), ass. (0->119)
	t201 = qJ(2) + qJ(3);
	t195 = cos(t201);
	t202 = sin(qJ(4));
	t277 = sin(qJ(1));
	t233 = t277 * t202;
	t203 = cos(qJ(4));
	t204 = cos(qJ(1));
	t254 = t204 * t203;
	t179 = t195 * t254 + t233;
	t173 = 0.1e1 / t179 ^ 2;
	t194 = sin(t201);
	t190 = t194 ^ 2;
	t200 = t204 ^ 2;
	t265 = t190 * t200;
	t235 = t173 * t265;
	t168 = 0.1e1 + t235;
	t228 = qJD(1) * t277;
	t196 = qJD(2) + qJD(3);
	t258 = t196 * t204;
	t237 = t194 * t258;
	t214 = t195 * t228 + t237;
	t227 = t277 * qJD(4);
	t255 = t204 * t202;
	t158 = (-qJD(4) * t195 + qJD(1)) * t255 + (t227 - t214) * t203;
	t172 = 0.1e1 / t179;
	t272 = t158 * t172 * t173;
	t222 = t265 * t272;
	t238 = t194 * t196 * t200;
	t280 = (-t222 + (-t190 * t204 * t228 + t195 * t238) * t173) / t168 ^ 2;
	t261 = t194 * t204;
	t175 = t195 * t233 + t254;
	t219 = t202 * t227;
	t250 = qJD(4) * t204;
	t230 = t203 * t250;
	t157 = t175 * qJD(1) - t195 * t230 + t202 * t237 - t219;
	t232 = t277 * t203;
	t178 = t195 * t255 - t232;
	t191 = 0.1e1 / t194;
	t192 = 0.1e1 / t194 ^ 2;
	t198 = 0.1e1 / t202 ^ 2;
	t251 = qJD(4) * t203;
	t231 = t198 * t251;
	t197 = 0.1e1 / t202;
	t259 = t196 * t197;
	t236 = t195 * t259;
	t264 = t191 * t197;
	t279 = (t191 * t231 + t192 * t236) * t178 + t157 * t264;
	t262 = t194 * t202;
	t167 = atan2(-t175, t262);
	t162 = cos(t167);
	t161 = sin(t167);
	t271 = t161 * t175;
	t156 = t162 * t262 - t271;
	t153 = 0.1e1 / t156;
	t154 = 0.1e1 / t156 ^ 2;
	t278 = 0.2e1 * t178;
	t170 = t175 ^ 2;
	t263 = t192 * t198;
	t169 = t170 * t263 + 0.1e1;
	t165 = 0.1e1 / t169;
	t260 = t195 * t196;
	t215 = t194 * t251 + t202 * t260;
	t240 = t175 * t263;
	t220 = t203 * t228;
	t234 = t277 * t194;
	t221 = t196 * t234;
	t253 = qJD(1) * t204;
	t159 = t203 * t227 * t195 - t220 + (t253 * t195 - t221 - t250) * t202;
	t242 = t159 * t264;
	t145 = (t215 * t240 - t242) * t165;
	t212 = -t145 * t175 + t215;
	t141 = (-t145 * t262 - t159) * t161 + t212 * t162;
	t155 = t153 * t154;
	t276 = t141 * t155;
	t193 = t191 / t190;
	t199 = t197 * t198;
	t275 = (t159 * t240 + (-t192 * t199 * t251 - t193 * t198 * t260) * t170) / t169 ^ 2;
	t274 = t154 * t178;
	t273 = t157 * t154;
	t270 = t161 * t178;
	t269 = t161 * t194;
	t268 = t162 * t175;
	t267 = t162 * t178;
	t266 = t162 * t195;
	t257 = t198 * t203;
	t256 = t204 * t153;
	t252 = qJD(4) * t202;
	t171 = t178 ^ 2;
	t151 = t171 * t154 + 0.1e1;
	t249 = 0.2e1 * (-t171 * t276 - t178 * t273) / t151 ^ 2;
	t248 = 0.2e1 * t280;
	t247 = -0.2e1 * t275;
	t246 = t155 * t278;
	t245 = t191 * t275;
	t244 = t154 * t270;
	t241 = t175 * t264;
	t239 = t192 * t195 * t197;
	t217 = t175 * t239 + t277;
	t152 = t217 * t165;
	t229 = t277 - t152;
	t226 = t153 * t249;
	t225 = t154 * t249;
	t224 = t261 * t278;
	t223 = t197 * t245;
	t177 = t195 * t232 - t255;
	t218 = t175 * t257 - t177 * t197;
	t216 = t173 * t177 * t204 - t277 * t172;
	t163 = 0.1e1 / t168;
	t160 = t179 * qJD(1) - t195 * t219 - t203 * t221 - t230;
	t149 = 0.1e1 / t151;
	t148 = t218 * t191 * t165;
	t144 = (-t161 + (t162 * t241 + t161) * t165) * t178;
	t143 = -t152 * t268 + (t229 * t269 + t266) * t202;
	t142 = t162 * t194 * t203 - t161 * t177 + (-t161 * t262 - t268) * t148;
	t140 = t217 * t247 + (t159 * t239 + t253 + (-t191 * t259 + (-t192 * t231 - 0.2e1 * t193 * t236) * t195) * t175) * t165;
	t138 = -0.2e1 * t218 * t245 + (-t218 * t192 * t260 + (t159 * t257 - t160 * t197 + (t177 * t257 + (-0.2e1 * t199 * t203 ^ 2 - t197) * t175) * qJD(4)) * t191) * t165;
	t137 = (t172 * t195 * t204 + t203 * t235) * t248 + (0.2e1 * t203 * t222 + t214 * t172 + ((t158 * t204 - 0.2e1 * t203 * t238) * t195 + (t200 * t252 + 0.2e1 * t204 * t220) * t190) * t173) * t163;
	t136 = t143 * t178 * t225 + (-(-t140 * t268 + (t145 * t271 - t159 * t162) * t152) * t274 + (t141 * t246 + t273) * t143 + (-t194 * t256 - (-t152 * t269 + t161 * t234 + t266) * t274) * t251) * t149 + (t226 * t261 + ((-t196 * t256 - (t229 * t196 - t145) * t244) * t195 + (t153 * t228 + (t204 * t141 - (-t140 + t253) * t270 - (t229 * t145 - t196) * t267) * t154) * t194) * t149) * t202;
	t1 = [t279 * t165 + t223 * t278, t140, t140, t138, 0, 0; t175 * t226 + (-t159 * t153 + (t141 * t175 + t144 * t157) * t154) * t149 + (t144 * t225 + (0.2e1 * t144 * t276 + (t157 * t165 - t157 - (-t145 * t165 * t241 + t247) * t178) * t154 * t161 + (-(-0.2e1 * t175 * t223 - t145) * t274 + (-(t145 + t242) * t178 + t279 * t175) * t154 * t165) * t162) * t149) * t178, t136, t136, (t142 * t274 - t153 * t179) * t249 + (t142 * t273 + t158 * t153 + (t142 * t246 - t179 * t154) * t141 - (-t194 * t252 + t203 * t260 - t138 * t175 - t148 * t159 + (-t148 * t262 - t177) * t145) * t154 * t267 - (-t160 + (-t138 * t202 - t145 * t203) * t194 - t212 * t148) * t244) * t149, 0, 0; t216 * t194 * t248 + (-t216 * t260 + ((qJD(1) * t172 + 0.2e1 * t177 * t272) * t204 + (-t277 * t158 - t160 * t204 + t177 * t228) * t173) * t194) * t163, t137, t137, t173 * t224 * t280 + (t224 * t272 + (t157 * t261 + (t194 * t228 - t195 * t258) * t178) * t173) * t163, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:16
	% EndTime: 2019-10-10 12:22:17
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (3091->99), mult. (3836->212), div. (757->12), fcn. (4483->9), ass. (0->95)
	t173 = sin(qJ(1));
	t169 = t173 ^ 2;
	t171 = qJ(2) + qJ(3);
	t166 = sin(t171);
	t162 = t166 ^ 2;
	t167 = cos(t171);
	t164 = 0.1e1 / t167 ^ 2;
	t216 = t162 * t164;
	t158 = t169 * t216 + 0.1e1;
	t161 = t166 * t162;
	t163 = 0.1e1 / t167;
	t168 = qJD(2) + qJD(3);
	t213 = t163 * t166;
	t183 = t168 * (t161 * t163 * t164 + t213);
	t175 = cos(qJ(1));
	t202 = qJD(1) * t175;
	t214 = t162 * t173;
	t186 = t202 * t214;
	t225 = (t164 * t186 + t169 * t183) / t158 ^ 2;
	t232 = -0.2e1 * t225;
	t172 = sin(qJ(4));
	t174 = cos(qJ(4));
	t209 = t168 * t175;
	t194 = t166 * t209;
	t203 = qJD(1) * t173;
	t205 = t173 * t174;
	t208 = t172 * t175;
	t135 = t174 * t194 - qJD(4) * t205 - t172 * t202 + (qJD(4) * t208 + t174 * t203) * t167;
	t151 = -t167 * t208 + t205;
	t145 = t151 ^ 2;
	t204 = t174 * t175;
	t206 = t173 * t172;
	t152 = t167 * t204 + t206;
	t147 = 0.1e1 / t152 ^ 2;
	t222 = t145 * t147;
	t144 = 0.1e1 + t222;
	t188 = qJD(1) * t167 - qJD(4);
	t189 = qJD(4) * t167 - qJD(1);
	t134 = -t189 * t204 + (t188 * t173 + t194) * t172;
	t219 = t147 * t151;
	t199 = t134 * t219;
	t146 = 0.1e1 / t152;
	t148 = t146 * t147;
	t221 = t145 * t148;
	t231 = (t135 * t221 + t199) / t144 ^ 2;
	t192 = 0.1e1 + t216;
	t230 = t173 * t192;
	t155 = 0.1e1 / t158;
	t193 = t166 * t202;
	t210 = t168 * t173;
	t195 = t164 * t210;
	t130 = ((t167 * t210 + t193) * t163 + t162 * t195) * t155;
	t229 = t130 - t210;
	t207 = t173 * t166;
	t157 = atan2(t207, t167);
	t154 = cos(t157);
	t153 = sin(t157);
	t196 = t153 * t207;
	t140 = t154 * t167 + t196;
	t137 = 0.1e1 / t140;
	t138 = 0.1e1 / t140 ^ 2;
	t228 = t155 - 0.1e1;
	t170 = t175 ^ 2;
	t215 = t162 * t170;
	t133 = t138 * t215 + 0.1e1;
	t211 = t167 * t168;
	t217 = t154 * t166;
	t125 = (t130 * t173 - t168) * t217 + (-t229 * t167 + t193) * t153;
	t226 = t125 * t137 * t138;
	t227 = (-t215 * t226 + (t166 * t170 * t211 - t186) * t138) / t133 ^ 2;
	t224 = t138 * t166;
	t223 = t138 * t175;
	t220 = t146 * t172;
	t218 = t153 * t173;
	t212 = t166 * t175;
	t201 = -0.2e1 * t226;
	t200 = 0.2e1 * t231;
	t198 = t135 * t148 * t151;
	t197 = t138 * t212;
	t191 = -0.2e1 * t166 * t227;
	t190 = t163 * t232;
	t187 = t154 * t155 * t162 * t163;
	t185 = t192 * t175;
	t184 = t174 * t219 + t220;
	t150 = -t167 * t205 + t208;
	t149 = t167 * t206 + t204;
	t142 = 0.1e1 / t144;
	t141 = t155 * t230;
	t131 = 0.1e1 / t133;
	t129 = (-t228 * t166 * t153 + t173 * t187) * t175;
	t127 = t167 * t218 - t217 + (-t153 * t167 + t154 * t207) * t141;
	t126 = t230 * t232 + (qJD(1) * t185 + 0.2e1 * t173 * t183) * t155;
	t123 = -0.2e1 * t184 * t212 * t231 + (t184 * t167 * t209 + (-t184 * t203 + ((qJD(4) * t146 + 0.2e1 * t198) * t174 + (t134 * t174 + (-qJD(4) * t151 + t135) * t172) * t147) * t175) * t166) * t142;
	t122 = 0.2e1 * (-t127 * t224 + t137 * t167) * t175 * t227 + ((t137 * t203 + (t127 * t168 + t125) * t223) * t167 + (t137 * t209 + (t126 * t154 * t173 + t229 * t153 + (-t130 * t218 + t153 * t168 + t154 * t202) * t141) * t197 + (-t138 * t203 + t175 * t201) * t127 + ((-t126 + t202) * t153 + ((t141 * t173 - 0.1e1) * t168 + (-t141 + t173) * t130) * t154) * t167 * t223) * t166) * t131;
	t1 = [t190 * t212 + (t168 * t185 - t203 * t213) * t155, t126, t126, 0, 0, 0; (t137 * t191 + (t137 * t211 + (-qJD(1) * t129 - t125) * t224) * t131) * t173 + (t138 * t191 * t129 + (((0.2e1 * t166 * t225 + t211 + (-t130 * t163 * t214 - t211) * t155) * t153 + (t190 * t214 + t130 * t166 + (t161 * t195 + (-t130 + 0.2e1 * t210) * t166) * t155) * t154) * t197 + (t138 * t211 + t166 * t201) * t129 + (t137 + ((-t169 + t170) * t187 + t228 * t196) * t138) * t166 * qJD(1)) * t131) * t175, t122, t122, 0, 0, 0; (-t146 * t149 + t150 * t219) * t200 + (-0.2e1 * t150 * t198 + t189 * t146 * t205 + (-t168 * t207 + t188 * t175) * t220 + (-t150 * t134 + t149 * t135 + (t188 * t204 - (t166 * t168 * t174 + t189 * t172) * t173) * t151) * t147) * t142, t123, t123, (t146 * t152 + t222) * t200 + (-0.2e1 * t199 + (-t147 * t152 + t146 - 0.2e1 * t221) * t135) * t142, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end