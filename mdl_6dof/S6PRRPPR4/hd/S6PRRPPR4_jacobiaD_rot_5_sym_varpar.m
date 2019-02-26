% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:00
% EndTime: 2019-02-26 20:00:01
% DurationCPUTime: 0.97s
% Computational Cost: add. (4171->104), mult. (12384->230), div. (514->12), fcn. (16056->13), ass. (0->104)
t208 = sin(qJ(2));
t210 = cos(qJ(2));
t259 = sin(pkin(10));
t261 = cos(pkin(6));
t230 = t261 * t259;
t260 = cos(pkin(10));
t200 = -t208 * t230 + t260 * t210;
t204 = sin(pkin(11));
t206 = cos(pkin(11));
t207 = sin(qJ(3));
t209 = cos(qJ(3));
t231 = t261 * t260;
t218 = -t208 * t231 - t259 * t210;
t205 = sin(pkin(6));
t237 = t205 * t260;
t220 = t207 * t237 + t209 * t218;
t229 = t259 * t208 - t210 * t231;
t173 = -t220 * t204 - t229 * t206;
t168 = t173 ^ 2;
t244 = t208 * t209;
t221 = t205 * t244 + t261 * t207;
t245 = t206 * t210;
t186 = t221 * t204 + t205 * t245;
t183 = 0.1e1 / t186 ^ 2;
t159 = t168 * t183 + 0.1e1;
t187 = t207 * t218 - t209 * t237;
t195 = t229 * qJD(2);
t243 = qJD(2) * t206;
t161 = (t187 * qJD(3) - t195 * t209) * t204 + t218 * t243;
t235 = t261 * t209;
t241 = qJD(3) * t207;
t242 = qJD(2) * t210;
t177 = qJD(3) * t204 * t235 + ((-t208 * t241 + t209 * t242) * t204 - t208 * t243) * t205;
t182 = 0.1e1 / t186;
t253 = t177 * t182 * t183;
t254 = t173 * t183;
t262 = -0.2e1 * (t161 * t254 - t168 * t253) / t159 ^ 2;
t160 = atan2(-t173, t186);
t152 = sin(t160);
t153 = cos(t160);
t227 = t152 * t173 - t153 * t186;
t148 = 0.1e1 / t227;
t217 = -t260 * t208 - t210 * t230;
t236 = t205 * t259;
t219 = -t200 * t209 - t207 * t236;
t176 = -t204 * t217 - t206 * t219;
t170 = 0.1e1 / t176;
t149 = 0.1e1 / t227 ^ 2;
t171 = 0.1e1 / t176 ^ 2;
t156 = 0.1e1 / t159;
t140 = (-t161 * t182 + t177 * t254) * t156;
t228 = -t152 * t186 - t153 * t173;
t137 = t228 * t140 - t152 * t161 + t153 * t177;
t258 = t137 * t148 * t149;
t175 = -t204 * t219 + t206 * t217;
t257 = t149 * t175;
t188 = -t200 * t207 + t209 * t236;
t196 = t217 * qJD(2);
t167 = t188 * qJD(3) + t196 * t209;
t197 = t200 * qJD(2);
t162 = t167 * t204 - t197 * t206;
t256 = t162 * t149;
t163 = t167 * t206 + t197 * t204;
t172 = t170 * t171;
t255 = t163 * t172;
t249 = t217 * t209;
t180 = t200 * t204 + t206 * t249;
t252 = t180 * t188;
t185 = t188 ^ 2;
t251 = t185 * t206;
t250 = t217 * t207;
t248 = t204 * t209;
t247 = t204 * t210;
t246 = t205 * t207;
t169 = t175 ^ 2;
t145 = t149 * t169 + 0.1e1;
t240 = 0.2e1 * (t169 * t258 + t175 * t256) / t145 ^ 2;
t158 = t171 * t185 + 0.1e1;
t166 = t219 * qJD(3) - t196 * t207;
t238 = t166 * t171 * t188;
t239 = 0.2e1 * (-t185 * t255 + t238) / t158 ^ 2;
t233 = -0.2e1 * t173 * t253;
t232 = -0.2e1 * t175 * t258;
t222 = t204 * t229;
t178 = t206 * t218 - t209 * t222;
t192 = (-t206 * t208 + t209 * t247) * t205;
t225 = -t178 * t182 + t192 * t254;
t201 = -t208 * t246 + t235;
t224 = -t182 * t187 + t201 * t254;
t223 = -t197 * t209 - t217 * t241;
t190 = -t221 * qJD(3) - t242 * t246;
t181 = (-t241 * t247 + (-t204 * t244 - t245) * qJD(2)) * t205;
t179 = -t200 * t206 + t217 * t248;
t165 = t220 * qJD(3) + t195 * t207;
t164 = t218 * qJD(2) * t248 + t195 * t206 + t222 * t241;
t154 = 0.1e1 / t158;
t143 = 0.1e1 / t145;
t142 = t224 * t204 * t156;
t141 = t225 * t156;
t139 = (-t152 * t187 + t153 * t201) * t204 + t228 * t142;
t138 = t228 * t141 - t152 * t178 + t153 * t192;
t136 = (t224 * t262 + (t201 * t233 - t165 * t182 + (t161 * t201 + t173 * t190 + t177 * t187) * t183) * t156) * t204;
t135 = t225 * t262 + (t192 * t233 - t164 * t182 + (t161 * t192 + t173 * t181 + t177 * t178) * t183) * t156;
t1 = [0, t135, t136, 0, 0, 0; 0 (t138 * t257 + t148 * t179) * t240 + (-(-t196 * t206 + t223 * t204) * t148 + t138 * t232 + (-t179 * t137 - t138 * t162 + (-(-t135 * t173 - t141 * t161 + t181 + (-t141 * t186 - t178) * t140) * t153 - (-t135 * t186 - t141 * t177 - t164 + (t141 * t173 - t192) * t140) * t152) * t175) * t149) * t143 (t148 * t188 * t204 + t139 * t257) * t240 + (-(t228 * t136 + (t227 * t140 - t152 * t177 - t153 * t161) * t142) * t257 + (t232 - t256) * t139 + (-t166 * t148 + (-t188 * t137 - (-t152 * t165 + t153 * t190 + (-t152 * t201 - t153 * t187) * t140) * t175) * t149) * t204) * t143, 0, 0, 0; 0 (t170 * t250 + t171 * t252) * t239 + (0.2e1 * t252 * t255 + (-qJD(3) * t249 + t197 * t207) * t170 + (t163 * t250 - (t196 * t204 + t223 * t206) * t188 - t180 * t166) * t171) * t154 (-t170 * t219 + t171 * t251) * t239 + (-0.2e1 * t206 * t238 - t167 * t170 + (-t171 * t219 + 0.2e1 * t172 * t251) * t163) * t154, 0, 0, 0;];
JaD_rot  = t1;
