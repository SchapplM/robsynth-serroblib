% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:43
% EndTime: 2019-02-26 22:34:44
% DurationCPUTime: 1.12s
% Computational Cost: add. (6985->123), mult. (8378->263), div. (1558->15), fcn. (10537->9), ass. (0->117)
t191 = qJ(3) + qJ(4);
t184 = cos(t191);
t269 = 0.2e1 * t184;
t183 = sin(t191);
t263 = sin(qJ(1));
t223 = t263 * t183;
t193 = cos(qJ(2));
t194 = cos(qJ(1));
t242 = t194 * t184;
t225 = t193 * t242;
t169 = t223 + t225;
t163 = 0.1e1 / t169 ^ 2;
t192 = sin(qJ(2));
t186 = t192 ^ 2;
t190 = t194 ^ 2;
t247 = t186 * t190;
t226 = t163 * t247;
t159 = 0.1e1 + t226;
t217 = qJD(1) * t263;
t240 = qJD(2) * t193;
t202 = t186 * t194 * t217 - t190 * t192 * t240;
t185 = qJD(3) + qJD(4);
t239 = qJD(2) * t194;
t219 = t192 * t239;
t205 = t193 * t217 + t219;
t221 = t263 * t185;
t243 = t194 * t183;
t148 = (-t185 * t193 + qJD(1)) * t243 + (t221 - t205) * t184;
t162 = 0.1e1 / t169;
t258 = t148 * t162 * t163;
t212 = t247 * t258;
t268 = (-t202 * t163 - t212) / t159 ^ 2;
t165 = t193 * t223 + t242;
t267 = t165 * t185;
t245 = t192 * t194;
t147 = t165 * qJD(1) - t185 * t225 + (t219 - t221) * t183;
t222 = t263 * t184;
t168 = t193 * t243 - t222;
t180 = 0.1e1 / t183;
t181 = 0.1e1 / t183 ^ 2;
t187 = 0.1e1 / t192;
t188 = 0.1e1 / t192 ^ 2;
t220 = t188 * t240;
t249 = t184 * t185;
t251 = t180 * t187;
t266 = (t181 * t187 * t249 + t180 * t220) * t168 + t147 * t251;
t246 = t192 * t183;
t155 = atan2(-t165, t246);
t152 = cos(t155);
t151 = sin(t155);
t257 = t151 * t165;
t146 = t152 * t246 - t257;
t143 = 0.1e1 / t146;
t144 = 0.1e1 / t146 ^ 2;
t265 = -0.2e1 * t165;
t264 = 0.2e1 * t168;
t160 = t165 ^ 2;
t250 = t181 * t188;
t156 = t160 * t250 + 0.1e1;
t153 = 0.1e1 / t156;
t248 = t184 * t192;
t206 = t183 * t240 + t185 * t248;
t229 = t165 * t250;
t224 = t192 * t263;
t210 = qJD(2) * t224;
t241 = qJD(1) * t194;
t149 = -t183 * t210 - t185 * t243 - t184 * t217 + (t183 * t241 + t184 * t221) * t193;
t231 = t149 * t251;
t135 = (t206 * t229 - t231) * t153;
t203 = -t135 * t165 + t206;
t130 = (-t135 * t246 - t149) * t151 + t203 * t152;
t145 = t143 * t144;
t262 = t130 * t145;
t182 = t180 * t181;
t189 = t187 / t186;
t227 = t188 * t249;
t261 = (t149 * t229 + (-t181 * t189 * t240 - t182 * t227) * t160) / t156 ^ 2;
t260 = t144 * t168;
t259 = t147 * t144;
t256 = t151 * t168;
t255 = t151 * t192;
t254 = t152 * t165;
t253 = t152 * t168;
t252 = t152 * t193;
t244 = t193 * t194;
t161 = t168 ^ 2;
t141 = t161 * t144 + 0.1e1;
t238 = 0.2e1 * (-t161 * t262 - t168 * t259) / t141 ^ 2;
t237 = -0.2e1 * t261;
t236 = 0.2e1 * t268;
t235 = t145 * t264;
t234 = t187 * t261;
t233 = t144 * t256;
t230 = t165 * t251;
t228 = t180 * t188 * t193;
t208 = t165 * t228 + t263;
t142 = t208 * t153;
t218 = t263 - t142;
t216 = t143 * t238;
t215 = t144 * t238;
t214 = t245 * t264;
t213 = t180 * t234;
t167 = t193 * t222 - t243;
t209 = t165 * t181 * t184 - t167 * t180;
t207 = t163 * t167 * t194 - t263 * t162;
t157 = 0.1e1 / t159;
t150 = t169 * qJD(1) - t184 * t210 - t267;
t139 = 0.1e1 / t141;
t138 = t209 * t187 * t153;
t134 = (-t151 + (t152 * t230 + t151) * t153) * t168;
t133 = -t142 * t254 + (t218 * t255 + t252) * t183;
t132 = t152 * t248 - t151 * t167 + (-t151 * t246 - t254) * t138;
t131 = t163 * t214 * t268 + (t214 * t258 + (t147 * t245 + (t192 * t217 - t193 * t239) * t168) * t163) * t157;
t129 = t208 * t237 + (t149 * t228 + t241 + (-t181 * t193 * t227 + (-0.2e1 * t189 * t193 ^ 2 - t187) * t180 * qJD(2)) * t165) * t153;
t127 = -0.2e1 * t209 * t234 + (-t209 * t220 + ((-t150 - t267) * t180 + (t182 * t249 * t265 + (t167 * t185 + t149) * t181) * t184) * t187) * t153;
t126 = (t132 * t260 - t143 * t169) * t238 + (t132 * t259 + t148 * t143 + (t132 * t235 - t169 * t144) * t130 - (t184 * t240 - t185 * t246 - t127 * t165 - t138 * t149 + (-t138 * t246 - t167) * t135) * t144 * t253 - (-t150 + (-t127 * t183 - t135 * t184) * t192 - t203 * t138) * t233) * t139;
t1 = [t266 * t153 + t213 * t264, t129, t127, t127, 0, 0; t165 * t216 + (-t149 * t143 + (t130 * t165 + t134 * t147) * t144) * t139 + (t134 * t215 + (0.2e1 * t134 * t262 + (t147 * t153 - t147 - (-t135 * t153 * t230 + t237) * t168) * t144 * t151 + (-(t213 * t265 - t135) * t260 + (-(t135 + t231) * t168 + t266 * t165) * t144 * t153) * t152) * t139) * t168, t133 * t168 * t215 + (-(-t129 * t254 + (t135 * t257 - t149 * t152) * t142) * t260 + (-t143 * t245 - (-t142 * t255 + t151 * t224 + t252) * t260) * t249 + (t130 * t235 + t259) * t133) * t139 + (t216 * t245 + ((-t143 * t239 - (t218 * qJD(2) - t135) * t233) * t193 + (t143 * t217 + (t194 * t130 - (-t129 + t241) * t256 - (t218 * t135 - qJD(2)) * t253) * t144) * t192) * t139) * t183, t126, t126, 0, 0; t207 * t192 * t236 + (-t207 * t240 + ((qJD(1) * t162 + 0.2e1 * t167 * t258) * t194 + (-t263 * t148 - t150 * t194 + t167 * t217) * t163) * t192) * t157 (t162 * t244 + t184 * t226) * t236 + (t212 * t269 + t205 * t162 + (t183 * t185 * t247 + t148 * t244 + t202 * t269) * t163) * t157, t131, t131, 0, 0;];
JaD_rot  = t1;
