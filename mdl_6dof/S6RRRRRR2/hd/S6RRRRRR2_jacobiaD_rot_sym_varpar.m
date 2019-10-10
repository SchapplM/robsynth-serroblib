% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR2
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
%   Wie in S6RRRRRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:03
	% EndTime: 2019-10-10 13:18:04
	% DurationCPUTime: 1.16s
	% Computational Cost: add. (7770->97), mult. (5101->208), div. (1026->12), fcn. (5942->9), ass. (0->95)
	t179 = sin(qJ(1));
	t241 = 0.2e1 * t179;
	t176 = t179 ^ 2;
	t175 = qJ(2) + qJ(3) + qJ(4);
	t172 = sin(t175);
	t167 = t172 ^ 2;
	t173 = cos(t175);
	t169 = 0.1e1 / t173 ^ 2;
	t226 = t167 * t169;
	t163 = t176 * t226 + 0.1e1;
	t161 = 0.1e1 / t163;
	t168 = 0.1e1 / t173;
	t181 = cos(qJ(1));
	t212 = qJD(1) * t181;
	t202 = t172 * t212;
	t174 = qJD(2) + qJD(3) + qJD(4);
	t220 = t174 * t179;
	t205 = t169 * t220;
	t135 = (-(-t173 * t220 - t202) * t168 + t167 * t205) * t161;
	t240 = t135 - t220;
	t180 = cos(qJ(5));
	t214 = t181 * t180;
	t178 = sin(qJ(5));
	t217 = t179 * t178;
	t159 = t173 * t214 + t217;
	t218 = t179 * t172;
	t160 = atan2(-t218, -t173);
	t151 = cos(t160);
	t150 = sin(t160);
	t207 = t150 * t218;
	t145 = -t151 * t173 - t207;
	t142 = 0.1e1 / t145;
	t153 = 0.1e1 / t159;
	t143 = 0.1e1 / t145 ^ 2;
	t154 = 0.1e1 / t159 ^ 2;
	t239 = t161 - 0.1e1;
	t231 = t151 * t172;
	t130 = (-t135 * t179 + t174) * t231 + (t240 * t173 - t202) * t150;
	t238 = t130 * t142 * t143;
	t190 = t173 * t217 + t214;
	t219 = t174 * t181;
	t203 = t172 * t219;
	t140 = t190 * qJD(1) - t159 * qJD(5) + t178 * t203;
	t215 = t181 * t178;
	t216 = t179 * t180;
	t158 = t173 * t215 - t216;
	t152 = t158 ^ 2;
	t149 = t152 * t154 + 0.1e1;
	t229 = t154 * t158;
	t196 = -qJD(1) * t173 + qJD(5);
	t197 = qJD(5) * t173 - qJD(1);
	t141 = -t197 * t215 + (t196 * t179 - t203) * t180;
	t235 = t141 * t153 * t154;
	t237 = (-t140 * t229 - t152 * t235) / t149 ^ 2;
	t166 = t172 * t167;
	t223 = t168 * t172;
	t189 = t174 * (t166 * t168 * t169 + t223);
	t224 = t167 * t179;
	t194 = t212 * t224;
	t236 = (t169 * t194 + t176 * t189) / t163 ^ 2;
	t234 = t143 * t172;
	t233 = t143 * t181;
	t232 = t150 * t179;
	t230 = t153 * t178;
	t228 = t158 * t180;
	t227 = t167 * t168;
	t177 = t181 ^ 2;
	t225 = t167 * t177;
	t222 = t172 * t181;
	t221 = t173 * t174;
	t213 = qJD(1) * t179;
	t138 = t143 * t225 + 0.1e1;
	t211 = 0.2e1 * (-t225 * t238 + (t172 * t177 * t221 - t194) * t143) / t138 ^ 2;
	t210 = 0.2e1 * t238;
	t209 = -0.2e1 * t237;
	t208 = t143 * t222;
	t206 = t158 * t235;
	t201 = 0.1e1 + t226;
	t200 = t172 * t211;
	t199 = -0.2e1 * t172 * t236;
	t198 = t236 * t241;
	t195 = t151 * t161 * t227;
	t193 = t201 * t181;
	t192 = t196 * t181;
	t191 = t154 * t228 - t230;
	t157 = -t173 * t216 + t215;
	t147 = 0.1e1 / t149;
	t146 = t201 * t179 * t161;
	t136 = 0.1e1 / t138;
	t134 = (t239 * t172 * t150 - t179 * t195) * t181;
	t132 = -t173 * t232 + t231 + (t150 * t173 - t151 * t218) * t146;
	t131 = -t201 * t198 + (qJD(1) * t193 + t189 * t241) * t161;
	t128 = t191 * t209 * t222 + (t191 * t173 * t219 + (-t191 * t213 + ((-qJD(5) * t153 - 0.2e1 * t206) * t180 + (-t140 * t180 + (-qJD(5) * t158 + t141) * t178) * t154) * t181) * t172) * t147;
	t127 = (t132 * t234 - t142 * t173) * t181 * t211 + ((-t142 * t213 + (-t132 * t174 - t130) * t233) * t173 + (-t142 * t219 - (-t131 * t151 * t179 - t240 * t150 + (t135 * t232 - t150 * t174 - t151 * t212) * t146) * t208 + (t143 * t213 + t181 * t210) * t132 - ((t131 - t212) * t150 + ((-t146 * t179 + 0.1e1) * t174 + (t146 - t179) * t135) * t151) * t173 * t233) * t172) * t136;
	t1 = [t181 * t168 * t199 + (t174 * t193 - t213 * t223) * t161, t131, t131, t131, 0, 0; (t142 * t200 + (-t142 * t221 + (qJD(1) * t134 + t130) * t234) * t136) * t179 + (t143 * t200 * t134 + (-((t199 - t221 + (t135 * t168 * t224 + t221) * t161) * t150 + (t198 * t227 - t135 * t172 + (-t166 * t205 + (t135 - 0.2e1 * t220) * t172) * t161) * t151) * t208 + (-t143 * t221 + t172 * t210) * t134 + (-t142 + ((-t176 + t177) * t195 + t239 * t207) * t143) * t172 * qJD(1)) * t136) * t181, t127, t127, t127, 0, 0; 0.2e1 * (t153 * t190 + t157 * t229) * t237 + (0.2e1 * t157 * t206 - t197 * t153 * t216 + (t174 * t218 + t192) * t230 + (t157 * t140 + t190 * t141 - t192 * t228 - (t172 * t174 * t180 + t197 * t178) * t158 * t179) * t154) * t147, t128, t128, t128, t209 + 0.2e1 * (-t140 * t154 * t147 + (-t147 * t235 - t154 * t237) * t158) * t158, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:03
	% EndTime: 2019-10-10 13:18:04
	% DurationCPUTime: 1.22s
	% Computational Cost: add. (8758->98), mult. (5316->205), div. (1044->12), fcn. (6173->9), ass. (0->98)
	t202 = sin(qJ(1));
	t263 = 0.2e1 * t202;
	t199 = t202 ^ 2;
	t197 = qJ(2) + qJ(3) + qJ(4);
	t192 = sin(t197);
	t188 = t192 ^ 2;
	t193 = cos(t197);
	t190 = 0.1e1 / t193 ^ 2;
	t248 = t188 * t190;
	t184 = t199 * t248 + 0.1e1;
	t182 = 0.1e1 / t184;
	t189 = 0.1e1 / t193;
	t203 = cos(qJ(1));
	t234 = qJD(1) * t203;
	t224 = t192 * t234;
	t194 = qJD(2) + qJD(3) + qJD(4);
	t242 = t194 * t202;
	t227 = t190 * t242;
	t155 = (-(-t193 * t242 - t224) * t189 + t188 * t227) * t182;
	t262 = t155 - t242;
	t201 = qJ(5) + qJ(6);
	t196 = cos(t201);
	t236 = t203 * t196;
	t195 = sin(t201);
	t239 = t202 * t195;
	t177 = t193 * t236 + t239;
	t240 = t202 * t192;
	t180 = atan2(-t240, -t193);
	t179 = cos(t180);
	t178 = sin(t180);
	t228 = t178 * t240;
	t165 = -t179 * t193 - t228;
	t162 = 0.1e1 / t165;
	t171 = 0.1e1 / t177;
	t163 = 0.1e1 / t165 ^ 2;
	t172 = 0.1e1 / t177 ^ 2;
	t261 = t182 - 0.1e1;
	t250 = t179 * t192;
	t150 = (-t155 * t202 + t194) * t250 + (t262 * t193 - t224) * t178;
	t260 = t150 * t162 * t163;
	t198 = qJD(5) + qJD(6);
	t213 = t193 * t239 + t236;
	t241 = t194 * t203;
	t225 = t192 * t241;
	t156 = t213 * qJD(1) - t177 * t198 + t195 * t225;
	t237 = t203 * t195;
	t238 = t202 * t196;
	t176 = t193 * t237 - t238;
	t170 = t176 ^ 2;
	t169 = t170 * t172 + 0.1e1;
	t253 = t172 * t176;
	t218 = -qJD(1) * t193 + t198;
	t219 = t193 * t198 - qJD(1);
	t157 = -t219 * t237 + (t218 * t202 - t225) * t196;
	t258 = t157 * t171 * t172;
	t259 = (-t156 * t253 - t170 * t258) / t169 ^ 2;
	t187 = t192 * t188;
	t245 = t189 * t192;
	t212 = t194 * (t187 * t189 * t190 + t245);
	t246 = t188 * t202;
	t216 = t234 * t246;
	t257 = (t190 * t216 + t199 * t212) / t184 ^ 2;
	t256 = t163 * t192;
	t255 = t163 * t203;
	t254 = t171 * t195;
	t252 = t176 * t196;
	t251 = t178 * t202;
	t249 = t188 * t189;
	t200 = t203 ^ 2;
	t247 = t188 * t200;
	t244 = t192 * t203;
	t243 = t193 * t194;
	t235 = qJD(1) * t202;
	t160 = t163 * t247 + 0.1e1;
	t233 = 0.2e1 * (-t247 * t260 + (t192 * t200 * t243 - t216) * t163) / t160 ^ 2;
	t232 = 0.2e1 * t260;
	t231 = -0.2e1 * t259;
	t230 = t163 * t244;
	t229 = t176 * t258;
	t223 = 0.1e1 + t248;
	t222 = t192 * t233;
	t221 = -0.2e1 * t192 * t257;
	t220 = t257 * t263;
	t217 = t179 * t182 * t249;
	t215 = t223 * t203;
	t214 = t172 * t252 - t254;
	t211 = t194 * t240 + t218 * t203;
	t175 = -t193 * t238 + t237;
	t167 = 0.1e1 / t169;
	t166 = t223 * t202 * t182;
	t158 = 0.1e1 / t160;
	t154 = (t261 * t192 * t178 - t202 * t217) * t203;
	t153 = -t193 * t251 + t250 + (t178 * t193 - t179 * t240) * t166;
	t151 = -t223 * t220 + (qJD(1) * t215 + t212 * t263) * t182;
	t148 = t231 + 0.2e1 * (-t156 * t172 * t167 + (-t167 * t258 - t172 * t259) * t176) * t176;
	t147 = t214 * t231 * t244 + (t214 * t193 * t241 + (-t214 * t235 + ((-t171 * t198 - 0.2e1 * t229) * t196 + (-t156 * t196 + (-t176 * t198 + t157) * t195) * t172) * t203) * t192) * t167;
	t146 = (t153 * t256 - t162 * t193) * t203 * t233 + ((-t162 * t235 + (-t153 * t194 - t150) * t255) * t193 + (-t162 * t241 - (-t151 * t179 * t202 - t262 * t178 + (t155 * t251 - t178 * t194 - t179 * t234) * t166) * t230 + (t163 * t235 + t203 * t232) * t153 - ((t151 - t234) * t178 + ((-t166 * t202 + 0.1e1) * t194 + (t166 - t202) * t155) * t179) * t193 * t255) * t192) * t158;
	t1 = [t203 * t189 * t221 + (t194 * t215 - t235 * t245) * t182, t151, t151, t151, 0, 0; (t162 * t222 + (-t162 * t243 + (qJD(1) * t154 + t150) * t256) * t158) * t202 + (t163 * t222 * t154 + (-((t221 - t243 + (t155 * t189 * t246 + t243) * t182) * t178 + (t220 * t249 - t155 * t192 + (-t187 * t227 + (t155 - 0.2e1 * t242) * t192) * t182) * t179) * t230 + (-t163 * t243 + t192 * t232) * t154 + (-t162 + ((-t199 + t200) * t217 + t261 * t228) * t163) * t192 * qJD(1)) * t158) * t203, t146, t146, t146, 0, 0; 0.2e1 * (t171 * t213 + t175 * t253) * t259 + (0.2e1 * t175 * t229 - t219 * t171 * t238 + t211 * t254 + (-t219 * t176 * t239 + t175 * t156 + t157 * t213 - t211 * t252) * t172) * t167, t147, t147, t147, t148, t148;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end