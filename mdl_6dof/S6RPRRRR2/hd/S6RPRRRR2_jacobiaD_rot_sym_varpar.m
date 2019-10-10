% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR2
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
%   Wie in S6RPRRRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:38
	% DurationCPUTime: 1.17s
	% Computational Cost: add. (5102->98), mult. (3810->203), div. (753->12), fcn. (4455->9), ass. (0->97)
	t173 = qJ(3) + qJ(4);
	t169 = sin(t173);
	t163 = t169 ^ 2;
	t170 = cos(t173);
	t165 = 0.1e1 / t170 ^ 2;
	t220 = t163 * t165;
	t172 = qJ(1) + pkin(11);
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
	t187 = t191 * t175;
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
	% StartTime: 2019-10-10 09:00:37
	% EndTime: 2019-10-10 09:00:38
	% DurationCPUTime: 1.20s
	% Computational Cost: add. (5971->100), mult. (4025->203), div. (771->12), fcn. (4686->9), ass. (0->100)
	t196 = qJ(3) + qJ(4);
	t189 = sin(t196);
	t182 = t189 ^ 2;
	t191 = cos(t196);
	t184 = 0.1e1 / t191 ^ 2;
	t241 = t182 * t184;
	t194 = qJ(1) + pkin(11);
	t186 = sin(t194);
	t260 = 0.2e1 * t186;
	t259 = t189 * t241;
	t187 = cos(t194);
	t195 = qJ(5) + qJ(6);
	t188 = sin(t195);
	t190 = cos(t195);
	t233 = t190 * t191;
	t169 = t186 * t188 + t187 * t233;
	t192 = qJD(5) + qJD(6);
	t213 = t191 * t192 - qJD(1);
	t193 = qJD(3) + qJD(4);
	t234 = t189 * t193;
	t258 = t213 * t188 + t190 * t234;
	t238 = t186 * t189;
	t173 = atan2(-t238, -t191);
	t171 = cos(t173);
	t170 = sin(t173);
	t222 = t170 * t238;
	t157 = -t171 * t191 - t222;
	t154 = 0.1e1 / t157;
	t163 = 0.1e1 / t169;
	t183 = 0.1e1 / t191;
	t155 = 0.1e1 / t157 ^ 2;
	t164 = 0.1e1 / t169 ^ 2;
	t257 = -0.2e1 * t189;
	t179 = t186 ^ 2;
	t176 = t179 * t241 + 0.1e1;
	t174 = 0.1e1 / t176;
	t256 = t174 - 0.1e1;
	t228 = qJD(1) * t189;
	t218 = t187 * t228;
	t232 = t191 * t193;
	t237 = t186 * t193;
	t147 = (-(-t186 * t232 - t218) * t183 + t237 * t241) * t174;
	t244 = t171 * t189;
	t142 = (-t147 * t186 + t193) * t244 + (-t218 + (t147 - t237) * t191) * t170;
	t255 = t142 * t154 * t155;
	t235 = t188 * t191;
	t206 = t186 * t235 + t187 * t190;
	t220 = t188 * t234;
	t148 = t206 * qJD(1) - t169 * t192 + t187 * t220;
	t168 = -t186 * t190 + t187 * t235;
	t162 = t168 ^ 2;
	t161 = t162 * t164 + 0.1e1;
	t246 = t164 * t168;
	t212 = -qJD(1) * t191 + t192;
	t208 = t212 * t190;
	t149 = t186 * t208 - t258 * t187;
	t251 = t149 * t163 * t164;
	t254 = (-t148 * t246 - t162 * t251) / t161 ^ 2;
	t253 = t147 * t170;
	t252 = t147 * t189;
	t205 = (t189 + t259) * t183 * t193;
	t229 = qJD(1) * t187;
	t210 = t182 * t186 * t229;
	t250 = (t179 * t205 + t184 * t210) / t176 ^ 2;
	t249 = t155 * t187;
	t248 = t155 * t189;
	t217 = 0.1e1 + t241;
	t160 = t217 * t186 * t174;
	t247 = t160 * t186;
	t245 = t170 * t191;
	t180 = t187 ^ 2;
	t243 = t180 * t182;
	t242 = t182 * t183;
	t240 = t183 * t186;
	t236 = t187 * t188;
	t231 = t193 * t154;
	t230 = qJD(1) * t186;
	t152 = t155 * t243 + 0.1e1;
	t227 = 0.2e1 * (-t243 * t255 + (t180 * t189 * t232 - t210) * t155) / t152 ^ 2;
	t226 = 0.2e1 * t255;
	t225 = 0.2e1 * t254;
	t224 = t187 * t248;
	t223 = t168 * t251;
	t216 = t189 * t227;
	t215 = t250 * t260;
	t214 = t250 * t257;
	t211 = t171 * t174 * t242;
	t209 = t217 * t187;
	t207 = -t163 * t188 + t190 * t246;
	t204 = t207 * t189;
	t167 = -t186 * t233 + t236;
	t158 = 0.1e1 / t161;
	t150 = 0.1e1 / t152;
	t146 = (t256 * t189 * t170 - t186 * t211) * t187;
	t145 = -t186 * t245 + t244 + (-t171 * t238 + t245) * t160;
	t143 = -t217 * t215 + (qJD(1) * t209 + t205 * t260) * t174;
	t140 = -0.2e1 * t254 + 0.2e1 * (-t148 * t164 * t158 + (-t158 * t251 - t164 * t254) * t168) * t168;
	t139 = -t187 * t204 * t225 + (-t204 * t230 + (t207 * t232 + ((-t163 * t192 - 0.2e1 * t223) * t190 + (-t148 * t190 + (-t168 * t192 + t149) * t188) * t164) * t189) * t187) * t158;
	t138 = (t145 * t248 - t154 * t191) * t187 * t227 + ((-t154 * t230 + (-t145 * t193 - t142) * t249) * t191 + (-t187 * t231 - (-t143 * t171 * t186 + t170 * t237 + t247 * t253 - t253 + (-t170 * t193 - t171 * t229) * t160) * t224 + (t155 * t230 + t187 * t226) * t145 - ((t143 - t229) * t170 + ((0.1e1 - t247) * t193 + (t160 - t186) * t147) * t171) * t191 * t249) * t189) * t150;
	t1 = [t187 * t183 * t214 + (t193 * t209 - t228 * t240) * t174, 0, t143, t143, 0, 0; (t154 * t216 + (-t191 * t231 + (qJD(1) * t146 + t142) * t248) * t150) * t186 + (t155 * t216 * t146 + (-((t214 - t232 + (t147 * t182 * t240 + t232) * t174) * t170 + (t215 * t242 - t252 + (t252 + (t257 - t259) * t237) * t174) * t171) * t224 + (-t155 * t232 + t189 * t226) * t146 + (-t154 + ((-t179 + t180) * t211 + t256 * t222) * t155) * t228) * t150) * t187, 0, t138, t138, 0, 0; (t163 * t206 + t167 * t246) * t225 + (0.2e1 * t167 * t223 + (t167 * t148 + t206 * t149 + (-t258 * t186 - t187 * t208) * t168) * t164 + (t212 * t236 + (-t213 * t190 + t220) * t186) * t163) * t158, 0, t139, t139, t140, t140;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end