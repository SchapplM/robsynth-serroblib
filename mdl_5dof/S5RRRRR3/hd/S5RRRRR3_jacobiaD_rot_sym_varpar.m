% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR3
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
%   Wie in S5RRRRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiaD_rot_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:40
	% EndTime: 2019-10-09 21:04:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:41
	% EndTime: 2019-10-09 21:04:42
	% DurationCPUTime: 1.14s
	% Computational Cost: add. (3645->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t172 = sin(qJ(1));
	t234 = 0.2e1 * t172;
	t168 = t172 ^ 2;
	t170 = qJ(2) + qJ(3);
	t165 = sin(t170);
	t161 = t165 ^ 2;
	t166 = cos(t170);
	t163 = 0.1e1 / t166 ^ 2;
	t219 = t161 * t163;
	t156 = t168 * t219 + 0.1e1;
	t154 = 0.1e1 / t156;
	t162 = 0.1e1 / t166;
	t174 = cos(qJ(1));
	t205 = qJD(1) * t174;
	t195 = t165 * t205;
	t167 = qJD(2) + qJD(3);
	t213 = t167 * t172;
	t198 = t163 * t213;
	t128 = (-(-t166 * t213 - t195) * t162 + t161 * t198) * t154;
	t233 = t128 - t213;
	t173 = cos(qJ(4));
	t207 = t174 * t173;
	t171 = sin(qJ(4));
	t210 = t172 * t171;
	t150 = t166 * t207 + t210;
	t211 = t172 * t165;
	t153 = atan2(-t211, -t166);
	t152 = cos(t153);
	t151 = sin(t153);
	t199 = t151 * t211;
	t138 = -t152 * t166 - t199;
	t135 = 0.1e1 / t138;
	t144 = 0.1e1 / t150;
	t136 = 0.1e1 / t138 ^ 2;
	t145 = 0.1e1 / t150 ^ 2;
	t232 = t154 - 0.1e1;
	t221 = t152 * t165;
	t123 = (-t128 * t172 + t167) * t221 + (t233 * t166 - t195) * t151;
	t231 = t123 * t135 * t136;
	t183 = t166 * t210 + t207;
	t212 = t167 * t174;
	t196 = t165 * t212;
	t132 = t183 * qJD(1) - t150 * qJD(4) + t171 * t196;
	t208 = t174 * t171;
	t209 = t172 * t173;
	t149 = t166 * t208 - t209;
	t143 = t149 ^ 2;
	t142 = t143 * t145 + 0.1e1;
	t224 = t145 * t149;
	t189 = -qJD(1) * t166 + qJD(4);
	t190 = qJD(4) * t166 - qJD(1);
	t133 = -t190 * t208 + (t189 * t172 - t196) * t173;
	t229 = t133 * t144 * t145;
	t230 = (-t132 * t224 - t143 * t229) / t142 ^ 2;
	t160 = t165 * t161;
	t216 = t162 * t165;
	t182 = t167 * (t160 * t162 * t163 + t216);
	t217 = t161 * t172;
	t187 = t205 * t217;
	t228 = (t163 * t187 + t168 * t182) / t156 ^ 2;
	t227 = t136 * t165;
	t226 = t136 * t174;
	t225 = t144 * t171;
	t223 = t149 * t173;
	t222 = t151 * t172;
	t220 = t161 * t162;
	t169 = t174 ^ 2;
	t218 = t161 * t169;
	t215 = t165 * t174;
	t214 = t166 * t167;
	t206 = qJD(1) * t172;
	t131 = t136 * t218 + 0.1e1;
	t204 = 0.2e1 * (-t218 * t231 + (t165 * t169 * t214 - t187) * t136) / t131 ^ 2;
	t203 = 0.2e1 * t231;
	t202 = -0.2e1 * t230;
	t201 = t136 * t215;
	t200 = t149 * t229;
	t194 = 0.1e1 + t219;
	t193 = t165 * t204;
	t192 = -0.2e1 * t165 * t228;
	t191 = t228 * t234;
	t188 = t152 * t154 * t220;
	t186 = t194 * t174;
	t185 = t189 * t174;
	t184 = t145 * t223 - t225;
	t148 = -t166 * t209 + t208;
	t140 = 0.1e1 / t142;
	t139 = t194 * t172 * t154;
	t129 = 0.1e1 / t131;
	t127 = (t232 * t165 * t151 - t172 * t188) * t174;
	t125 = -t166 * t222 + t221 + (t151 * t166 - t152 * t211) * t139;
	t124 = -t194 * t191 + (qJD(1) * t186 + t182 * t234) * t154;
	t121 = t184 * t202 * t215 + (t184 * t166 * t212 + (-t184 * t206 + ((-qJD(4) * t144 - 0.2e1 * t200) * t173 + (-t132 * t173 + (-qJD(4) * t149 + t133) * t171) * t145) * t174) * t165) * t140;
	t120 = (t125 * t227 - t135 * t166) * t174 * t204 + ((-t135 * t206 + (-t125 * t167 - t123) * t226) * t166 + (-t135 * t212 - (-t124 * t152 * t172 - t233 * t151 + (t128 * t222 - t151 * t167 - t152 * t205) * t139) * t201 + (t136 * t206 + t174 * t203) * t125 - ((t124 - t205) * t151 + ((-t139 * t172 + 0.1e1) * t167 + (t139 - t172) * t128) * t152) * t166 * t226) * t165) * t129;
	t1 = [t174 * t162 * t192 + (t167 * t186 - t206 * t216) * t154, t124, t124, 0, 0; (t135 * t193 + (-t135 * t214 + (qJD(1) * t127 + t123) * t227) * t129) * t172 + (t136 * t193 * t127 + (-((t192 - t214 + (t128 * t162 * t217 + t214) * t154) * t151 + (t191 * t220 - t128 * t165 + (-t160 * t198 + (t128 - 0.2e1 * t213) * t165) * t154) * t152) * t201 + (-t136 * t214 + t165 * t203) * t127 + (-t135 + ((-t168 + t169) * t188 + t232 * t199) * t136) * t165 * qJD(1)) * t129) * t174, t120, t120, 0, 0; 0.2e1 * (t144 * t183 + t148 * t224) * t230 + (0.2e1 * t148 * t200 - t190 * t144 * t209 + (t167 * t211 + t185) * t225 + (t148 * t132 + t183 * t133 - t185 * t223 - (t165 * t167 * t173 + t190 * t171) * t149 * t172) * t145) * t140, t121, t121, t202 + 0.2e1 * (-t132 * t145 * t140 + (-t140 * t229 - t145 * t230) * t149) * t149, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:04:41
	% EndTime: 2019-10-09 21:04:42
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (4423->98), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->98)
	t196 = sin(qJ(1));
	t257 = 0.2e1 * t196;
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
	t236 = t191 * t196;
	t221 = t184 * t236;
	t149 = (-(-t189 * t236 - t218) * t183 + t182 * t221) * t176;
	t256 = t149 - t236;
	t194 = qJ(4) + qJ(5);
	t188 = cos(t194);
	t230 = t197 * t188;
	t186 = sin(t194);
	t234 = t196 * t186;
	t171 = t189 * t230 + t234;
	t233 = t196 * t187;
	t175 = atan2(-t233, -t189);
	t173 = cos(t175);
	t172 = sin(t175);
	t222 = t172 * t233;
	t159 = -t173 * t189 - t222;
	t156 = 0.1e1 / t159;
	t165 = 0.1e1 / t171;
	t157 = 0.1e1 / t159 ^ 2;
	t166 = 0.1e1 / t171 ^ 2;
	t255 = t176 - 0.1e1;
	t244 = t173 * t187;
	t144 = (-t149 * t196 + t191) * t244 + (t256 * t189 - t218) * t172;
	t254 = t144 * t156 * t157;
	t190 = qJD(4) + qJD(5);
	t207 = t189 * t234 + t230;
	t235 = t191 * t197;
	t219 = t187 * t235;
	t150 = t207 * qJD(1) - t171 * t190 + t186 * t219;
	t231 = t197 * t186;
	t232 = t196 * t188;
	t170 = t189 * t231 - t232;
	t164 = t170 ^ 2;
	t163 = t164 * t166 + 0.1e1;
	t247 = t166 * t170;
	t212 = -qJD(1) * t189 + t190;
	t213 = t189 * t190 - qJD(1);
	t151 = -t213 * t231 + (t212 * t196 - t219) * t188;
	t252 = t151 * t165 * t166;
	t253 = (-t150 * t247 - t164 * t252) / t163 ^ 2;
	t181 = t187 * t182;
	t239 = t183 * t187;
	t206 = t191 * (t181 * t183 * t184 + t239);
	t240 = t182 * t196;
	t210 = t228 * t240;
	t251 = (t184 * t210 + t192 * t206) / t178 ^ 2;
	t250 = t157 * t187;
	t249 = t157 * t197;
	t248 = t165 * t186;
	t246 = t170 * t188;
	t245 = t172 * t196;
	t243 = t182 * t183;
	t193 = t197 ^ 2;
	t241 = t182 * t193;
	t238 = t187 * t197;
	t237 = t189 * t191;
	t229 = qJD(1) * t196;
	t154 = t157 * t241 + 0.1e1;
	t227 = 0.2e1 * (-t241 * t254 + (t187 * t193 * t237 - t210) * t157) / t154 ^ 2;
	t226 = 0.2e1 * t254;
	t225 = -0.2e1 * t253;
	t224 = t157 * t238;
	t223 = t170 * t252;
	t217 = 0.1e1 + t242;
	t216 = t187 * t227;
	t215 = -0.2e1 * t187 * t251;
	t214 = t251 * t257;
	t211 = t173 * t176 * t243;
	t209 = t217 * t197;
	t208 = t166 * t246 - t248;
	t205 = t191 * t233 + t212 * t197;
	t169 = -t189 * t232 + t231;
	t162 = t217 * t196 * t176;
	t160 = 0.1e1 / t163;
	t152 = 0.1e1 / t154;
	t148 = (t255 * t187 * t172 - t196 * t211) * t197;
	t147 = -t189 * t245 + t244 + (t172 * t189 - t173 * t233) * t162;
	t145 = -t217 * t214 + (qJD(1) * t209 + t206 * t257) * t176;
	t142 = t225 + 0.2e1 * (-t150 * t166 * t160 + (-t160 * t252 - t166 * t253) * t170) * t170;
	t141 = t208 * t225 * t238 + (t208 * t189 * t235 + (-t208 * t229 + ((-t165 * t190 - 0.2e1 * t223) * t188 + (-t150 * t188 + (-t170 * t190 + t151) * t186) * t166) * t197) * t187) * t160;
	t140 = (t147 * t250 - t156 * t189) * t197 * t227 + ((-t156 * t229 + (-t147 * t191 - t144) * t249) * t189 + (-t156 * t235 - (-t145 * t173 * t196 - t256 * t172 + (t149 * t245 - t172 * t191 - t173 * t228) * t162) * t224 + (t157 * t229 + t197 * t226) * t147 - ((t145 - t228) * t172 + ((-t162 * t196 + 0.1e1) * t191 + (t162 - t196) * t149) * t173) * t189 * t249) * t187) * t152;
	t1 = [t197 * t183 * t215 + (t191 * t209 - t229 * t239) * t176, t145, t145, 0, 0; (t156 * t216 + (-t156 * t237 + (qJD(1) * t148 + t144) * t250) * t152) * t196 + (t157 * t216 * t148 + (-((t215 - t237 + (t149 * t183 * t240 + t237) * t176) * t172 + (t214 * t243 - t149 * t187 + (-t181 * t221 + (t149 - 0.2e1 * t236) * t187) * t176) * t173) * t224 + (-t157 * t237 + t187 * t226) * t148 + (-t156 + ((-t192 + t193) * t211 + t255 * t222) * t157) * t187 * qJD(1)) * t152) * t197, t140, t140, 0, 0; 0.2e1 * (t165 * t207 + t169 * t247) * t253 + (0.2e1 * t169 * t223 - t213 * t165 * t232 + t205 * t248 + (-t213 * t170 * t234 + t169 * t150 + t151 * t207 - t205 * t246) * t166) * t160, t141, t141, t142, t142;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end