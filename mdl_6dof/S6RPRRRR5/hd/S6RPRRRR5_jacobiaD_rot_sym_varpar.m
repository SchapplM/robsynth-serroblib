% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR5
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
%   Wie in S6RPRRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:08
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (5587->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t174 = sin(qJ(1));
	t236 = 0.2e1 * t174;
	t171 = t174 ^ 2;
	t169 = pkin(11) + qJ(3) + qJ(4);
	t166 = sin(t169);
	t162 = t166 ^ 2;
	t167 = cos(t169);
	t164 = 0.1e1 / t167 ^ 2;
	t221 = t162 * t164;
	t158 = t171 * t221 + 0.1e1;
	t156 = 0.1e1 / t158;
	t163 = 0.1e1 / t167;
	t176 = cos(qJ(1));
	t207 = qJD(1) * t176;
	t197 = t166 * t207;
	t170 = qJD(3) + qJD(4);
	t215 = t170 * t174;
	t200 = t164 * t215;
	t130 = (-(-t167 * t215 - t197) * t163 + t162 * t200) * t156;
	t235 = t130 - t215;
	t175 = cos(qJ(5));
	t209 = t175 * t176;
	t173 = sin(qJ(5));
	t211 = t174 * t173;
	t155 = t167 * t209 + t211;
	t212 = t174 * t166;
	t151 = atan2(-t212, -t167);
	t146 = cos(t151);
	t145 = sin(t151);
	t202 = t145 * t212;
	t140 = -t146 * t167 - t202;
	t137 = 0.1e1 / t140;
	t148 = 0.1e1 / t155;
	t138 = 0.1e1 / t140 ^ 2;
	t149 = 0.1e1 / t155 ^ 2;
	t234 = t156 - 0.1e1;
	t226 = t146 * t166;
	t125 = (-t130 * t174 + t170) * t226 + (t235 * t167 - t197) * t145;
	t233 = t125 * t137 * t138;
	t185 = t167 * t211 + t209;
	t214 = t170 * t176;
	t198 = t166 * t214;
	t135 = t185 * qJD(1) - t155 * qJD(5) + t173 * t198;
	t210 = t174 * t175;
	t213 = t173 * t176;
	t154 = t167 * t213 - t210;
	t147 = t154 ^ 2;
	t144 = t147 * t149 + 0.1e1;
	t224 = t149 * t154;
	t191 = -qJD(1) * t167 + qJD(5);
	t192 = qJD(5) * t167 - qJD(1);
	t136 = -t192 * t213 + (t191 * t174 - t198) * t175;
	t230 = t136 * t148 * t149;
	t232 = (-t135 * t224 - t147 * t230) / t144 ^ 2;
	t161 = t166 * t162;
	t218 = t163 * t166;
	t184 = t170 * (t161 * t163 * t164 + t218);
	t219 = t162 * t174;
	t189 = t207 * t219;
	t231 = (t164 * t189 + t171 * t184) / t158 ^ 2;
	t229 = t138 * t166;
	t228 = t138 * t176;
	t227 = t145 * t174;
	t225 = t148 * t173;
	t223 = t154 * t175;
	t222 = t162 * t163;
	t172 = t176 ^ 2;
	t220 = t162 * t172;
	t217 = t166 * t176;
	t216 = t167 * t170;
	t208 = qJD(1) * t174;
	t133 = t138 * t220 + 0.1e1;
	t206 = 0.2e1 * (-t220 * t233 + (t166 * t172 * t216 - t189) * t138) / t133 ^ 2;
	t205 = 0.2e1 * t233;
	t204 = -0.2e1 * t232;
	t203 = t138 * t217;
	t201 = t154 * t230;
	t196 = 0.1e1 + t221;
	t195 = t166 * t206;
	t194 = -0.2e1 * t166 * t231;
	t193 = t231 * t236;
	t190 = t146 * t156 * t222;
	t188 = t196 * t176;
	t187 = t191 * t176;
	t186 = t149 * t223 - t225;
	t153 = -t167 * t210 + t213;
	t142 = 0.1e1 / t144;
	t141 = t196 * t174 * t156;
	t131 = 0.1e1 / t133;
	t129 = (t234 * t166 * t145 - t174 * t190) * t176;
	t127 = -t167 * t227 + t226 + (t145 * t167 - t146 * t212) * t141;
	t126 = -t196 * t193 + (qJD(1) * t188 + t184 * t236) * t156;
	t123 = t186 * t204 * t217 + (t186 * t167 * t214 + (-t186 * t208 + ((-qJD(5) * t148 - 0.2e1 * t201) * t175 + (-t135 * t175 + (-qJD(5) * t154 + t136) * t173) * t149) * t176) * t166) * t142;
	t122 = (t127 * t229 - t137 * t167) * t176 * t206 + ((-t137 * t208 + (-t127 * t170 - t125) * t228) * t167 + (-t137 * t214 - (-t126 * t146 * t174 - t235 * t145 + (t130 * t227 - t145 * t170 - t146 * t207) * t141) * t203 + (t138 * t208 + t176 * t205) * t127 - ((t126 - t207) * t145 + ((-t141 * t174 + 0.1e1) * t170 + (t141 - t174) * t130) * t146) * t167 * t228) * t166) * t131;
	t1 = [t163 * t176 * t194 + (t170 * t188 - t208 * t218) * t156, 0, t126, t126, 0, 0; (t137 * t195 + (-t137 * t216 + (qJD(1) * t129 + t125) * t229) * t131) * t174 + (t138 * t195 * t129 + (-((t194 - t216 + (t130 * t163 * t219 + t216) * t156) * t145 + (t193 * t222 - t130 * t166 + (-t161 * t200 + (t130 - 0.2e1 * t215) * t166) * t156) * t146) * t203 + (-t138 * t216 + t166 * t205) * t129 + (-t137 + ((-t171 + t172) * t190 + t234 * t202) * t138) * t166 * qJD(1)) * t131) * t176, 0, t122, t122, 0, 0; 0.2e1 * (t148 * t185 + t153 * t224) * t232 + (0.2e1 * t153 * t201 - t192 * t148 * t210 + (t170 * t212 + t187) * t225 + (t153 * t135 + t185 * t136 - t187 * t223 - (t166 * t170 * t175 + t192 * t173) * t154 * t174) * t149) * t142, 0, t123, t123, t204 + 0.2e1 * (-t135 * t142 * t149 + (-t142 * t230 - t149 * t232) * t154) * t154, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:08
	% DurationCPUTime: 1.10s
	% Computational Cost: add. (6414->98), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->98)
	t197 = sin(qJ(1));
	t258 = 0.2e1 * t197;
	t194 = t197 ^ 2;
	t189 = pkin(11) + qJ(3) + qJ(4);
	t187 = sin(t189);
	t183 = t187 ^ 2;
	t188 = cos(t189);
	t185 = 0.1e1 / t188 ^ 2;
	t243 = t183 * t185;
	t179 = t194 * t243 + 0.1e1;
	t176 = 0.1e1 / t179;
	t184 = 0.1e1 / t188;
	t198 = cos(qJ(1));
	t229 = qJD(1) * t198;
	t219 = t187 * t229;
	t193 = qJD(3) + qJD(4);
	t237 = t193 * t197;
	t222 = t185 * t237;
	t150 = (-(-t188 * t237 - t219) * t184 + t183 * t222) * t176;
	t257 = t150 - t237;
	t196 = qJ(5) + qJ(6);
	t191 = cos(t196);
	t231 = t198 * t191;
	t190 = sin(t196);
	t234 = t197 * t190;
	t172 = t188 * t231 + t234;
	t235 = t197 * t187;
	t175 = atan2(-t235, -t188);
	t174 = cos(t175);
	t173 = sin(t175);
	t223 = t173 * t235;
	t160 = -t174 * t188 - t223;
	t157 = 0.1e1 / t160;
	t166 = 0.1e1 / t172;
	t158 = 0.1e1 / t160 ^ 2;
	t167 = 0.1e1 / t172 ^ 2;
	t256 = t176 - 0.1e1;
	t245 = t174 * t187;
	t145 = (-t150 * t197 + t193) * t245 + (t257 * t188 - t219) * t173;
	t255 = t145 * t157 * t158;
	t192 = qJD(5) + qJD(6);
	t208 = t188 * t234 + t231;
	t236 = t193 * t198;
	t220 = t187 * t236;
	t151 = t208 * qJD(1) - t172 * t192 + t190 * t220;
	t232 = t198 * t190;
	t233 = t197 * t191;
	t171 = t188 * t232 - t233;
	t165 = t171 ^ 2;
	t164 = t165 * t167 + 0.1e1;
	t248 = t167 * t171;
	t213 = -qJD(1) * t188 + t192;
	t214 = t188 * t192 - qJD(1);
	t152 = -t214 * t232 + (t213 * t197 - t220) * t191;
	t253 = t152 * t166 * t167;
	t254 = (-t151 * t248 - t165 * t253) / t164 ^ 2;
	t182 = t187 * t183;
	t240 = t184 * t187;
	t207 = t193 * (t182 * t184 * t185 + t240);
	t241 = t183 * t197;
	t211 = t229 * t241;
	t252 = (t185 * t211 + t194 * t207) / t179 ^ 2;
	t251 = t158 * t187;
	t250 = t158 * t198;
	t249 = t166 * t190;
	t247 = t171 * t191;
	t246 = t173 * t197;
	t244 = t183 * t184;
	t195 = t198 ^ 2;
	t242 = t183 * t195;
	t239 = t187 * t198;
	t238 = t188 * t193;
	t230 = qJD(1) * t197;
	t155 = t158 * t242 + 0.1e1;
	t228 = 0.2e1 * (-t242 * t255 + (t187 * t195 * t238 - t211) * t158) / t155 ^ 2;
	t227 = 0.2e1 * t255;
	t226 = -0.2e1 * t254;
	t225 = t158 * t239;
	t224 = t171 * t253;
	t218 = 0.1e1 + t243;
	t217 = t187 * t228;
	t216 = -0.2e1 * t187 * t252;
	t215 = t252 * t258;
	t212 = t174 * t176 * t244;
	t210 = t218 * t198;
	t209 = t167 * t247 - t249;
	t206 = t193 * t235 + t213 * t198;
	t170 = -t188 * t233 + t232;
	t162 = 0.1e1 / t164;
	t161 = t218 * t197 * t176;
	t153 = 0.1e1 / t155;
	t149 = (t256 * t187 * t173 - t197 * t212) * t198;
	t148 = -t188 * t246 + t245 + (t173 * t188 - t174 * t235) * t161;
	t146 = -t218 * t215 + (qJD(1) * t210 + t207 * t258) * t176;
	t143 = t226 + 0.2e1 * (-t151 * t167 * t162 + (-t162 * t253 - t167 * t254) * t171) * t171;
	t142 = t209 * t226 * t239 + (t209 * t188 * t236 + (-t209 * t230 + ((-t166 * t192 - 0.2e1 * t224) * t191 + (-t151 * t191 + (-t171 * t192 + t152) * t190) * t167) * t198) * t187) * t162;
	t141 = (t148 * t251 - t157 * t188) * t198 * t228 + ((-t157 * t230 + (-t148 * t193 - t145) * t250) * t188 + (-t157 * t236 - (-t146 * t174 * t197 - t257 * t173 + (t150 * t246 - t173 * t193 - t174 * t229) * t161) * t225 + (t158 * t230 + t198 * t227) * t148 - ((t146 - t229) * t173 + ((-t161 * t197 + 0.1e1) * t193 + (t161 - t197) * t150) * t174) * t188 * t250) * t187) * t153;
	t1 = [t198 * t184 * t216 + (t193 * t210 - t230 * t240) * t176, 0, t146, t146, 0, 0; (t157 * t217 + (-t157 * t238 + (qJD(1) * t149 + t145) * t251) * t153) * t197 + (t158 * t217 * t149 + (-((t216 - t238 + (t150 * t184 * t241 + t238) * t176) * t173 + (t215 * t244 - t150 * t187 + (-t182 * t222 + (t150 - 0.2e1 * t237) * t187) * t176) * t174) * t225 + (-t158 * t238 + t187 * t227) * t149 + (-t157 + ((-t194 + t195) * t212 + t256 * t223) * t158) * t187 * qJD(1)) * t153) * t198, 0, t141, t141, 0, 0; 0.2e1 * (t166 * t208 + t170 * t248) * t254 + (0.2e1 * t170 * t224 - t214 * t166 * t233 + t206 * t249 + (-t214 * t171 * t234 + t170 * t151 + t152 * t208 - t206 * t247) * t167) * t162, 0, t142, t142, t143, t143;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end