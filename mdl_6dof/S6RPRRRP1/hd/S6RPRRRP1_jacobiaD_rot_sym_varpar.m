% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP1
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
%   Wie in S6RPRRRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:45
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:25
	% DurationCPUTime: 1.17s
	% Computational Cost: add. (5102->98), mult. (3810->203), div. (753->12), fcn. (4455->9), ass. (0->97)
	t173 = qJ(3) + qJ(4);
	t169 = sin(t173);
	t163 = t169 ^ 2;
	t170 = cos(t173);
	t165 = 0.1e1 / t170 ^ 2;
	t220 = t163 * t165;
	t172 = qJ(1) + pkin(10);
	t167 = sin(t172);
	t239 = 0.2e1 * t167;
	t238 = t169 * t220;
	t168 = cos(t172);
	t174 = sin(qJ(5));
	t175 = cos(qJ(5));
	t211 = t170 * t175;
	t150 = t167 * t174 + t168 * t211;
	t192 = qJD(5) * t170 - qJD(1);
	t171 = qJD(3) + qJD(4);
	t214 = t169 * t171;
	t237 = t192 * t174 + t175 * t214;
	t218 = t167 * t169;
	t153 = atan2(-t218, -t170);
	t152 = cos(t153);
	t151 = sin(t153);
	t201 = t151 * t218;
	t138 = -t152 * t170 - t201;
	t135 = 0.1e1 / t138;
	t144 = 0.1e1 / t150;
	t164 = 0.1e1 / t170;
	t136 = 0.1e1 / t138 ^ 2;
	t145 = 0.1e1 / t150 ^ 2;
	t236 = -0.2e1 * t169;
	t160 = t167 ^ 2;
	t156 = t160 * t220 + 0.1e1;
	t154 = 0.1e1 / t156;
	t235 = t154 - 0.1e1;
	t207 = qJD(1) * t169;
	t197 = t168 * t207;
	t213 = t170 * t171;
	t217 = t167 * t171;
	t128 = (-(-t167 * t213 - t197) * t164 + t217 * t220) * t154;
	t223 = t152 * t169;
	t123 = (-t128 * t167 + t171) * t223 + (-t197 + (t128 - t217) * t170) * t151;
	t234 = t123 * t135 * t136;
	t212 = t170 * t174;
	t185 = t167 * t212 + t168 * t175;
	t199 = t174 * t214;
	t132 = t185 * qJD(1) - t150 * qJD(5) + t168 * t199;
	t149 = -t167 * t175 + t168 * t212;
	t143 = t149 ^ 2;
	t142 = t143 * t145 + 0.1e1;
	t225 = t145 * t149;
	t191 = -qJD(1) * t170 + qJD(5);
	t187 = t175 * t191;
	t133 = t167 * t187 - t237 * t168;
	t230 = t133 * t144 * t145;
	t233 = (-t132 * t225 - t143 * t230) / t142 ^ 2;
	t232 = t128 * t151;
	t231 = t128 * t169;
	t184 = (t169 + t238) * t164 * t171;
	t208 = qJD(1) * t168;
	t189 = t163 * t167 * t208;
	t229 = (t160 * t184 + t165 * t189) / t156 ^ 2;
	t228 = t136 * t168;
	t227 = t136 * t169;
	t196 = 0.1e1 + t220;
	t139 = t196 * t167 * t154;
	t226 = t139 * t167;
	t224 = t151 * t170;
	t161 = t168 ^ 2;
	t222 = t161 * t163;
	t221 = t163 * t164;
	t219 = t164 * t167;
	t215 = t168 * t174;
	t210 = t171 * t135;
	t209 = qJD(1) * t167;
	t131 = t136 * t222 + 0.1e1;
	t206 = 0.2e1 * (-t222 * t234 + (t161 * t169 * t213 - t189) * t136) / t131 ^ 2;
	t205 = 0.2e1 * t234;
	t204 = 0.2e1 * t233;
	t203 = t168 * t227;
	t202 = t149 * t230;
	t195 = t169 * t206;
	t194 = t229 * t239;
	t193 = t229 * t236;
	t190 = t152 * t154 * t221;
	t188 = t196 * t168;
	t186 = -t144 * t174 + t175 * t225;
	t183 = t186 * t169;
	t148 = -t167 * t211 + t215;
	t140 = 0.1e1 / t142;
	t129 = 0.1e1 / t131;
	t127 = (t235 * t169 * t151 - t167 * t190) * t168;
	t126 = -t167 * t224 + t223 + (-t152 * t218 + t224) * t139;
	t124 = -t196 * t194 + (qJD(1) * t188 + t184 * t239) * t154;
	t121 = -t168 * t183 * t204 + (-t183 * t209 + (t186 * t213 + ((-qJD(5) * t144 - 0.2e1 * t202) * t175 + (-t132 * t175 + (-qJD(5) * t149 + t133) * t174) * t145) * t169) * t168) * t140;
	t120 = (t126 * t227 - t135 * t170) * t168 * t206 + ((-t135 * t209 + (-t126 * t171 - t123) * t228) * t170 + (-t168 * t210 - (-t124 * t152 * t167 + t151 * t217 + t226 * t232 - t232 + (-t151 * t171 - t152 * t208) * t139) * t203 + (t136 * t209 + t168 * t205) * t126 - ((t124 - t208) * t151 + ((0.1e1 - t226) * t171 + (t139 - t167) * t128) * t152) * t170 * t228) * t169) * t129;
	t1 = [t168 * t164 * t193 + (t171 * t188 - t207 * t219) * t154, 0, t124, t124, 0, 0; (t135 * t195 + (-t170 * t210 + (qJD(1) * t127 + t123) * t227) * t129) * t167 + (t136 * t195 * t127 + (-((t193 - t213 + (t128 * t163 * t219 + t213) * t154) * t151 + (t194 * t221 - t231 + (t231 + (t236 - t238) * t217) * t154) * t152) * t203 + (-t136 * t213 + t169 * t205) * t127 + (-t135 + ((-t160 + t161) * t190 + t235 * t201) * t136) * t207) * t129) * t168, 0, t120, t120, 0, 0; (t144 * t185 + t148 * t225) * t204 + (0.2e1 * t148 * t202 + (t148 * t132 + t185 * t133 + (-t237 * t167 - t168 * t187) * t149) * t145 + (t191 * t215 + (-t192 * t175 + t199) * t167) * t144) * t140, 0, t121, t121, -0.2e1 * t233 + 0.2e1 * (-t132 * t145 * t140 + (-t140 * t230 - t145 * t233) * t149) * t149, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:25
	% DurationCPUTime: 1.67s
	% Computational Cost: add. (8522->127), mult. (8382->276), div. (1515->15), fcn. (10508->9), ass. (0->119)
	t198 = qJ(3) + qJ(4);
	t193 = cos(t198);
	t199 = sin(qJ(5));
	t247 = qJ(1) + pkin(10);
	t229 = sin(t247);
	t222 = t229 * t199;
	t191 = cos(t247);
	t200 = cos(qJ(5));
	t257 = t191 * t200;
	t175 = t193 * t257 + t222;
	t169 = 0.1e1 / t175 ^ 2;
	t186 = t191 ^ 2;
	t192 = sin(t198);
	t187 = t192 ^ 2;
	t262 = t186 * t187;
	t238 = t169 * t262;
	t159 = 0.1e1 + t238;
	t219 = qJD(1) * t229;
	t215 = t193 * t219;
	t218 = t229 * qJD(5);
	t194 = qJD(3) + qJD(4);
	t252 = t194 * t200;
	t154 = (-t215 + t218) * t200 + (-t192 * t252 + (-qJD(5) * t193 + qJD(1)) * t199) * t191;
	t168 = 0.1e1 / t175;
	t269 = t154 * t168 * t169;
	t224 = t262 * t269;
	t255 = t193 * t194;
	t276 = (-t224 + (t186 * t192 * t255 - t187 * t191 * t219) * t169) / t159 ^ 2;
	t259 = t191 * t192;
	t171 = t193 * t222 + t257;
	t214 = t199 * t218;
	t248 = qJD(5) * t200;
	t230 = t191 * t248;
	t258 = t191 * t194;
	t234 = t192 * t258;
	t153 = t171 * qJD(1) - t193 * t230 + t199 * t234 - t214;
	t221 = t229 * t200;
	t254 = t193 * t199;
	t174 = t191 * t254 - t221;
	t188 = 0.1e1 / t192;
	t189 = 0.1e1 / t192 ^ 2;
	t196 = 0.1e1 / t199 ^ 2;
	t231 = t196 * t248;
	t195 = 0.1e1 / t199;
	t253 = t194 * t195;
	t233 = t193 * t253;
	t261 = t188 * t195;
	t275 = (t188 * t231 + t189 * t233) * t174 + t153 * t261;
	t256 = t192 * t199;
	t164 = atan2(-t171, t256);
	t161 = cos(t164);
	t160 = sin(t164);
	t268 = t160 * t171;
	t152 = t161 * t256 - t268;
	t149 = 0.1e1 / t152;
	t150 = 0.1e1 / t152 ^ 2;
	t274 = 0.2e1 * t174;
	t166 = t171 ^ 2;
	t260 = t189 * t196;
	t165 = t166 * t260 + 0.1e1;
	t162 = 0.1e1 / t165;
	t212 = t192 * t248 + t194 * t254;
	t236 = t171 * t260;
	t216 = t200 * t219;
	t223 = t229 * t192;
	t217 = t194 * t223;
	t249 = qJD(5) * t199;
	t250 = qJD(1) * t191;
	t155 = -t199 * t217 - t191 * t249 - t216 + (t199 * t250 + t200 * t218) * t193;
	t239 = t155 * t261;
	t141 = (t212 * t236 - t239) * t162;
	t208 = -t141 * t171 + t212;
	t137 = (-t141 * t256 - t155) * t160 + t208 * t161;
	t151 = t149 * t150;
	t273 = t137 * t151;
	t190 = t188 / t187;
	t197 = t195 * t196;
	t272 = (t155 * t236 + (-t189 * t197 * t248 - t190 * t196 * t255) * t166) / t165 ^ 2;
	t271 = t150 * t174;
	t270 = t153 * t150;
	t267 = t160 * t174;
	t266 = t160 * t192;
	t265 = t161 * t171;
	t264 = t161 * t174;
	t263 = t161 * t193;
	t251 = t196 * t200;
	t167 = t174 ^ 2;
	t147 = t150 * t167 + 0.1e1;
	t246 = 0.2e1 * (-t167 * t273 - t174 * t270) / t147 ^ 2;
	t245 = 0.2e1 * t276;
	t244 = -0.2e1 * t272;
	t243 = t151 * t274;
	t242 = t188 * t272;
	t241 = t150 * t267;
	t237 = t171 * t261;
	t235 = t189 * t193 * t195;
	t232 = t193 * t252;
	t228 = t149 * t246;
	t227 = t150 * t246;
	t226 = t259 * t274;
	t225 = t195 * t242;
	t211 = t171 * t235 + t229;
	t148 = t211 * t162;
	t220 = t229 - t148;
	t173 = -t191 * t199 + t193 * t221;
	t213 = t171 * t251 - t173 * t195;
	t210 = t169 * t173 * t191 - t229 * t168;
	t157 = 0.1e1 / t159;
	t156 = t175 * qJD(1) - t193 * t214 - t200 * t217 - t230;
	t145 = 0.1e1 / t147;
	t144 = t213 * t188 * t162;
	t140 = (-t160 + (t161 * t237 + t160) * t162) * t174;
	t139 = -t148 * t265 + (t220 * t266 + t263) * t199;
	t138 = t161 * t192 * t200 - t160 * t173 + (-t160 * t256 - t265) * t144;
	t136 = t211 * t244 + (t155 * t235 + t250 + (-t188 * t253 + (-t189 * t231 - 0.2e1 * t190 * t233) * t193) * t171) * t162;
	t134 = -0.2e1 * t213 * t242 + (-t213 * t189 * t255 + (t155 * t251 - t156 * t195 + (t173 * t251 + (-0.2e1 * t197 * t200 ^ 2 - t195) * t171) * qJD(5)) * t188) * t162;
	t133 = (t168 * t191 * t193 + t200 * t238) * t245 + (0.2e1 * t200 * t224 + (t215 + t234) * t168 + ((t154 * t193 + 0.2e1 * t187 * t216) * t191 + (t187 * t249 - 0.2e1 * t192 * t232) * t186) * t169) * t157;
	t132 = t139 * t174 * t227 + (-(-t136 * t265 + (t141 * t268 - t155 * t161) * t148) * t271 + (t137 * t243 + t270) * t139 + (-t149 * t259 - (-t148 * t266 + t160 * t223 + t263) * t271) * t248) * t145 + (t228 * t259 + ((-t149 * t258 - (t220 * t194 - t141) * t241) * t193 + (t149 * t219 + (t191 * t137 - (-t136 + t250) * t267 - (t220 * t141 - t194) * t264) * t150) * t192) * t145) * t199;
	t1 = [t162 * t275 + t225 * t274, 0, t136, t136, t134, 0; t171 * t228 + (-t155 * t149 + (t137 * t171 + t140 * t153) * t150) * t145 + (t140 * t227 + (0.2e1 * t140 * t273 + (t153 * t162 - t153 - (-t141 * t162 * t237 + t244) * t174) * t150 * t160 + (-(-0.2e1 * t171 * t225 - t141) * t271 + (-(t141 + t239) * t174 + t275 * t171) * t150 * t162) * t161) * t145) * t174, 0, t132, t132, (t138 * t271 - t149 * t175) * t246 + (t138 * t270 + t154 * t149 + (t138 * t243 - t175 * t150) * t137 - (-t192 * t249 + t232 - t134 * t171 - t144 * t155 + (-t144 * t256 - t173) * t141) * t150 * t264 - (-t156 + (-t134 * t199 - t141 * t200) * t192 - t208 * t144) * t241) * t145, 0; t210 * t192 * t245 + (-t210 * t255 + ((qJD(1) * t168 + 0.2e1 * t173 * t269) * t191 + (-t229 * t154 - t156 * t191 + t173 * t219) * t169) * t192) * t157, 0, t133, t133, t169 * t226 * t276 + (t226 * t269 + (t153 * t259 + (-t191 * t255 + t192 * t219) * t174) * t169) * t157, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end