% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP2
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
%   Wie in S6RRPRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:52
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (5587->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t176 = sin(qJ(1));
	t238 = 0.2e1 * t176;
	t173 = t176 ^ 2;
	t171 = qJ(2) + pkin(10) + qJ(4);
	t168 = sin(t171);
	t164 = t168 ^ 2;
	t169 = cos(t171);
	t166 = 0.1e1 / t169 ^ 2;
	t223 = t164 * t166;
	t160 = t173 * t223 + 0.1e1;
	t158 = 0.1e1 / t160;
	t165 = 0.1e1 / t169;
	t178 = cos(qJ(1));
	t209 = qJD(1) * t178;
	t199 = t168 * t209;
	t172 = qJD(2) + qJD(4);
	t217 = t172 * t176;
	t202 = t166 * t217;
	t132 = (-(-t169 * t217 - t199) * t165 + t164 * t202) * t158;
	t237 = t132 - t217;
	t177 = cos(qJ(5));
	t211 = t177 * t178;
	t175 = sin(qJ(5));
	t213 = t176 * t175;
	t157 = t169 * t211 + t213;
	t214 = t176 * t168;
	t153 = atan2(-t214, -t169);
	t148 = cos(t153);
	t147 = sin(t153);
	t204 = t147 * t214;
	t142 = -t148 * t169 - t204;
	t139 = 0.1e1 / t142;
	t150 = 0.1e1 / t157;
	t140 = 0.1e1 / t142 ^ 2;
	t151 = 0.1e1 / t157 ^ 2;
	t236 = t158 - 0.1e1;
	t228 = t148 * t168;
	t127 = (-t132 * t176 + t172) * t228 + (t237 * t169 - t199) * t147;
	t235 = t127 * t139 * t140;
	t187 = t169 * t213 + t211;
	t216 = t172 * t178;
	t200 = t168 * t216;
	t137 = t187 * qJD(1) - t157 * qJD(5) + t175 * t200;
	t212 = t176 * t177;
	t215 = t175 * t178;
	t156 = t169 * t215 - t212;
	t149 = t156 ^ 2;
	t146 = t149 * t151 + 0.1e1;
	t226 = t151 * t156;
	t193 = -qJD(1) * t169 + qJD(5);
	t194 = qJD(5) * t169 - qJD(1);
	t138 = -t194 * t215 + (t193 * t176 - t200) * t177;
	t232 = t138 * t150 * t151;
	t234 = (-t137 * t226 - t149 * t232) / t146 ^ 2;
	t163 = t168 * t164;
	t220 = t165 * t168;
	t186 = t172 * (t163 * t165 * t166 + t220);
	t221 = t164 * t176;
	t191 = t209 * t221;
	t233 = (t166 * t191 + t173 * t186) / t160 ^ 2;
	t231 = t140 * t168;
	t230 = t140 * t178;
	t229 = t147 * t176;
	t227 = t150 * t175;
	t225 = t156 * t177;
	t224 = t164 * t165;
	t174 = t178 ^ 2;
	t222 = t164 * t174;
	t219 = t168 * t178;
	t218 = t169 * t172;
	t210 = qJD(1) * t176;
	t135 = t140 * t222 + 0.1e1;
	t208 = 0.2e1 * (-t222 * t235 + (t168 * t174 * t218 - t191) * t140) / t135 ^ 2;
	t207 = 0.2e1 * t235;
	t206 = -0.2e1 * t234;
	t205 = t140 * t219;
	t203 = t156 * t232;
	t198 = 0.1e1 + t223;
	t197 = t168 * t208;
	t196 = -0.2e1 * t168 * t233;
	t195 = t233 * t238;
	t192 = t148 * t158 * t224;
	t190 = t198 * t178;
	t189 = t193 * t178;
	t188 = t151 * t225 - t227;
	t155 = -t169 * t212 + t215;
	t144 = 0.1e1 / t146;
	t143 = t198 * t176 * t158;
	t133 = 0.1e1 / t135;
	t131 = (t236 * t168 * t147 - t176 * t192) * t178;
	t129 = -t169 * t229 + t228 + (t147 * t169 - t148 * t214) * t143;
	t128 = -t198 * t195 + (qJD(1) * t190 + t186 * t238) * t158;
	t125 = t188 * t206 * t219 + (t188 * t169 * t216 + (-t188 * t210 + ((-qJD(5) * t150 - 0.2e1 * t203) * t177 + (-t137 * t177 + (-qJD(5) * t156 + t138) * t175) * t151) * t178) * t168) * t144;
	t124 = (t129 * t231 - t139 * t169) * t178 * t208 + ((-t139 * t210 + (-t129 * t172 - t127) * t230) * t169 + (-t139 * t216 - (-t128 * t148 * t176 - t237 * t147 + (t132 * t229 - t147 * t172 - t148 * t209) * t143) * t205 + (t140 * t210 + t178 * t207) * t129 - ((t128 - t209) * t147 + ((-t143 * t176 + 0.1e1) * t172 + (t143 - t176) * t132) * t148) * t169 * t230) * t168) * t133;
	t1 = [t165 * t178 * t196 + (t172 * t190 - t210 * t220) * t158, t128, 0, t128, 0, 0; (t139 * t197 + (-t139 * t218 + (qJD(1) * t131 + t127) * t231) * t133) * t176 + (t140 * t197 * t131 + (-((t196 - t218 + (t132 * t165 * t221 + t218) * t158) * t147 + (t195 * t224 - t132 * t168 + (-t163 * t202 + (t132 - 0.2e1 * t217) * t168) * t158) * t148) * t205 + (-t140 * t218 + t168 * t207) * t131 + (-t139 + ((-t173 + t174) * t192 + t236 * t204) * t140) * t168 * qJD(1)) * t133) * t178, t124, 0, t124, 0, 0; 0.2e1 * (t150 * t187 + t155 * t226) * t234 + (0.2e1 * t155 * t203 - t194 * t150 * t212 + (t172 * t214 + t189) * t227 + (t155 * t137 + t187 * t138 - t189 * t225 - (t168 * t172 * t177 + t194 * t175) * t156 * t176) * t151) * t144, t125, 0, t125, t206 + 0.2e1 * (-t137 * t144 * t151 + (-t144 * t232 - t151 * t234) * t156) * t156, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:53
	% DurationCPUTime: 1.62s
	% Computational Cost: add. (8326->125), mult. (8382->273), div. (1515->15), fcn. (10508->9), ass. (0->119)
	t199 = qJ(2) + pkin(10) + qJ(4);
	t196 = cos(t199);
	t205 = sin(qJ(5));
	t280 = sin(qJ(1));
	t236 = t280 * t205;
	t206 = cos(qJ(5));
	t207 = cos(qJ(1));
	t257 = t207 * t206;
	t182 = t196 * t257 + t236;
	t176 = 0.1e1 / t182 ^ 2;
	t195 = sin(t199);
	t191 = t195 ^ 2;
	t204 = t207 ^ 2;
	t268 = t191 * t204;
	t244 = t176 * t268;
	t171 = 0.1e1 + t244;
	t231 = qJD(1) * t280;
	t200 = qJD(2) + qJD(4);
	t261 = t200 * t207;
	t239 = t195 * t261;
	t217 = t196 * t231 + t239;
	t230 = t280 * qJD(5);
	t258 = t207 * t205;
	t161 = (-qJD(5) * t196 + qJD(1)) * t258 + (t230 - t217) * t206;
	t175 = 0.1e1 / t182;
	t275 = t161 * t175 * t176;
	t225 = t268 * t275;
	t240 = t195 * t200 * t204;
	t283 = (-t225 + (-t191 * t207 * t231 + t196 * t240) * t176) / t171 ^ 2;
	t264 = t195 * t207;
	t178 = t196 * t236 + t257;
	t222 = t205 * t230;
	t253 = qJD(5) * t207;
	t233 = t206 * t253;
	t160 = t178 * qJD(1) - t196 * t233 + t205 * t239 - t222;
	t235 = t280 * t206;
	t181 = t196 * t258 - t235;
	t192 = 0.1e1 / t195;
	t193 = 0.1e1 / t195 ^ 2;
	t202 = 0.1e1 / t205 ^ 2;
	t254 = qJD(5) * t206;
	t234 = t202 * t254;
	t201 = 0.1e1 / t205;
	t262 = t200 * t201;
	t238 = t196 * t262;
	t267 = t192 * t201;
	t282 = (t192 * t234 + t193 * t238) * t181 + t160 * t267;
	t265 = t195 * t205;
	t170 = atan2(-t178, t265);
	t165 = cos(t170);
	t164 = sin(t170);
	t274 = t164 * t178;
	t159 = t165 * t265 - t274;
	t156 = 0.1e1 / t159;
	t157 = 0.1e1 / t159 ^ 2;
	t281 = 0.2e1 * t181;
	t173 = t178 ^ 2;
	t266 = t193 * t202;
	t172 = t173 * t266 + 0.1e1;
	t168 = 0.1e1 / t172;
	t263 = t196 * t200;
	t218 = t195 * t254 + t205 * t263;
	t242 = t178 * t266;
	t223 = t206 * t231;
	t237 = t280 * t195;
	t224 = t200 * t237;
	t256 = qJD(1) * t207;
	t162 = t206 * t230 * t196 - t223 + (t256 * t196 - t224 - t253) * t205;
	t245 = t162 * t267;
	t148 = (t218 * t242 - t245) * t168;
	t215 = -t148 * t178 + t218;
	t144 = (-t148 * t265 - t162) * t164 + t215 * t165;
	t158 = t156 * t157;
	t279 = t144 * t158;
	t194 = t192 / t191;
	t203 = t201 * t202;
	t278 = (t162 * t242 + (-t193 * t203 * t254 - t194 * t202 * t263) * t173) / t172 ^ 2;
	t277 = t157 * t181;
	t276 = t160 * t157;
	t273 = t164 * t181;
	t272 = t164 * t195;
	t271 = t165 * t178;
	t270 = t165 * t181;
	t269 = t165 * t196;
	t260 = t202 * t206;
	t259 = t207 * t156;
	t255 = qJD(5) * t205;
	t174 = t181 ^ 2;
	t154 = t157 * t174 + 0.1e1;
	t252 = 0.2e1 * (-t174 * t279 - t181 * t276) / t154 ^ 2;
	t251 = 0.2e1 * t283;
	t250 = -0.2e1 * t278;
	t249 = t158 * t281;
	t248 = t192 * t278;
	t247 = t157 * t273;
	t243 = t178 * t267;
	t241 = t193 * t196 * t201;
	t220 = t178 * t241 + t280;
	t155 = t220 * t168;
	t232 = t280 - t155;
	t229 = t156 * t252;
	t228 = t157 * t252;
	t227 = t264 * t281;
	t226 = t201 * t248;
	t180 = t196 * t235 - t258;
	t221 = t178 * t260 - t180 * t201;
	t219 = t176 * t180 * t207 - t280 * t175;
	t166 = 0.1e1 / t171;
	t163 = t182 * qJD(1) - t196 * t222 - t206 * t224 - t233;
	t152 = 0.1e1 / t154;
	t151 = t221 * t192 * t168;
	t147 = (-t164 + (t165 * t243 + t164) * t168) * t181;
	t146 = -t155 * t271 + (t232 * t272 + t269) * t205;
	t145 = t165 * t195 * t206 - t164 * t180 + (-t164 * t265 - t271) * t151;
	t143 = t220 * t250 + (t162 * t241 + t256 + (-t192 * t262 + (-t193 * t234 - 0.2e1 * t194 * t238) * t196) * t178) * t168;
	t141 = -0.2e1 * t221 * t248 + (-t221 * t193 * t263 + (t162 * t260 - t163 * t201 + (t180 * t260 + (-0.2e1 * t203 * t206 ^ 2 - t201) * t178) * qJD(5)) * t192) * t168;
	t140 = (t175 * t196 * t207 + t206 * t244) * t251 + (0.2e1 * t206 * t225 + t217 * t175 + ((t161 * t207 - 0.2e1 * t206 * t240) * t196 + (t204 * t255 + 0.2e1 * t207 * t223) * t191) * t176) * t166;
	t139 = t146 * t181 * t228 + (-(-t143 * t271 + (t148 * t274 - t162 * t165) * t155) * t277 + (t144 * t249 + t276) * t146 + (-t195 * t259 - (-t155 * t272 + t164 * t237 + t269) * t277) * t254) * t152 + (t229 * t264 + ((-t200 * t259 - (t232 * t200 - t148) * t247) * t196 + (t156 * t231 + (t207 * t144 - (-t143 + t256) * t273 - (t232 * t148 - t200) * t270) * t157) * t195) * t152) * t205;
	t1 = [t168 * t282 + t226 * t281, t143, 0, t143, t141, 0; t178 * t229 + (-t162 * t156 + (t144 * t178 + t147 * t160) * t157) * t152 + (t147 * t228 + (0.2e1 * t147 * t279 + (t160 * t168 - t160 - (-t148 * t168 * t243 + t250) * t181) * t157 * t164 + (-(-0.2e1 * t178 * t226 - t148) * t277 + (-(t148 + t245) * t181 + t282 * t178) * t157 * t168) * t165) * t152) * t181, t139, 0, t139, (t145 * t277 - t156 * t182) * t252 + (t145 * t276 + t161 * t156 + (t145 * t249 - t157 * t182) * t144 - (-t195 * t255 + t206 * t263 - t141 * t178 - t151 * t162 + (-t151 * t265 - t180) * t148) * t157 * t270 - (-t163 + (-t141 * t205 - t148 * t206) * t195 - t215 * t151) * t247) * t152, 0; t219 * t195 * t251 + (-t219 * t263 + ((qJD(1) * t175 + 0.2e1 * t180 * t275) * t207 + (-t280 * t161 - t163 * t207 + t180 * t231) * t176) * t195) * t166, t140, 0, t140, t176 * t227 * t283 + (t227 * t275 + (t160 * t264 + (t195 * t231 - t196 * t261) * t181) * t176) * t166, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end