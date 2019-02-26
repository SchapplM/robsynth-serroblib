% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR10_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:01
% EndTime: 2019-02-26 21:43:02
% DurationCPUTime: 1.15s
% Computational Cost: add. (7177->123), mult. (11003->252), div. (691->12), fcn. (14028->11), ass. (0->106)
t207 = cos(pkin(6));
t209 = sin(qJ(1));
t210 = cos(qJ(2));
t266 = cos(qJ(1));
t231 = t266 * qJD(2);
t232 = t266 * qJD(1);
t206 = sin(pkin(6));
t237 = t206 * t266;
t208 = sin(qJ(2));
t250 = t209 * t208;
t238 = t207 * t250;
t248 = qJD(2) * t208;
t271 = -qJD(1) * t238 - t209 * t248 + (t207 * t231 + t232) * t210 - qJD(4) * t237;
t205 = pkin(11) + qJ(4);
t203 = sin(t205);
t204 = cos(t205);
t236 = t266 * t208;
t249 = t209 * t210;
t218 = -t207 * t236 - t249;
t174 = -t203 * t218 + t204 * t237;
t253 = t206 * t208;
t184 = t203 * t253 - t207 * t204;
t161 = atan2(-t174, t184);
t156 = sin(t161);
t157 = cos(t161);
t151 = -t156 * t174 + t157 * t184;
t149 = 0.1e1 / t151 ^ 2;
t235 = t266 * t210;
t192 = t235 - t238;
t252 = t206 * t209;
t222 = -t192 * t203 + t204 * t252;
t172 = t222 ^ 2;
t147 = t172 * t149 + 0.1e1;
t219 = -t207 * t249 - t236;
t166 = qJD(1) * t218 + qJD(2) * t219;
t180 = t192 * t204 + t203 * t252;
t226 = t206 * t232;
t152 = qJD(4) * t180 + t166 * t203 - t204 * t226;
t261 = t152 * t149;
t171 = t174 ^ 2;
t182 = 0.1e1 / t184 ^ 2;
t160 = t171 * t182 + 0.1e1;
t158 = 0.1e1 / t160;
t234 = qJD(1) * t252;
t247 = qJD(4) * t204;
t154 = t271 * t203 - t204 * t234 - t218 * t247;
t185 = t207 * t203 + t204 * t253;
t251 = t206 * t210;
t233 = qJD(2) * t251;
t169 = qJD(4) * t185 + t203 * t233;
t181 = 0.1e1 / t184;
t256 = t174 * t182;
t224 = -t154 * t181 + t169 * t256;
t140 = t224 * t158;
t225 = -t156 * t184 - t157 * t174;
t136 = t140 * t225 - t156 * t154 + t157 * t169;
t148 = 0.1e1 / t151;
t150 = t148 * t149;
t264 = t136 * t150;
t246 = 0.2e1 * (-t172 * t264 - t222 * t261) / t147 ^ 2;
t270 = t169 * t182;
t228 = t207 * t235;
t189 = t228 - t250;
t220 = -t181 * t189 + t251 * t256;
t269 = t203 * t220;
t155 = t203 * (qJD(4) * t218 + t234) + t271 * t204;
t186 = 0.1e1 / t219;
t187 = 0.1e1 / t219 ^ 2;
t268 = -0.2e1 * t174;
t267 = -0.2e1 * t222;
t258 = t181 * t270;
t263 = (t154 * t256 - t171 * t258) / t160 ^ 2;
t262 = t149 * t222;
t260 = t156 * t222;
t259 = t157 * t222;
t257 = t174 * t181;
t255 = t180 * t187;
t254 = t219 * t203;
t245 = -0.2e1 * t263;
t153 = qJD(4) * t222 + t166 * t204 + t203 * t226;
t173 = t180 ^ 2;
t164 = t173 * t187 + 0.1e1;
t165 = -qJD(1) * t228 - t210 * t231 + (qJD(2) * t207 + qJD(1)) * t250;
t188 = t186 * t187;
t244 = 0.2e1 * (-t173 * t188 * t165 + t153 * t255) / t164 ^ 2;
t243 = t150 * t267;
t242 = 0.2e1 * t180 * t188;
t241 = t181 * t263;
t240 = t149 * t260;
t239 = t149 * t259;
t230 = t258 * t268;
t176 = -t203 * t237 - t204 * t218;
t223 = -t176 * t181 + t185 * t256;
t217 = -t156 + (t157 * t257 + t156) * t158;
t170 = -qJD(4) * t184 + t204 * t233;
t167 = qJD(1) * t219 + qJD(2) * t218;
t162 = 0.1e1 / t164;
t145 = 0.1e1 / t147;
t144 = t158 * t269;
t141 = t223 * t158;
t139 = t217 * t222;
t138 = (-t156 * t189 + t157 * t251) * t203 + t225 * t144;
t137 = t141 * t225 - t156 * t176 + t157 * t185;
t135 = t223 * t245 + (t185 * t230 - t155 * t181 + (t154 * t185 + t169 * t176 + t170 * t174) * t182) * t158;
t133 = t245 * t269 + (t220 * t247 + (t230 * t251 - t167 * t181 + (t169 * t189 + (t154 * t210 - t174 * t248) * t206) * t182) * t203) * t158;
t1 = [t241 * t267 + (-t152 * t181 - t222 * t270) * t158, t133, 0, t135, 0, 0; t174 * t148 * t246 + (-t154 * t148 + (t136 * t174 + t139 * t152) * t149) * t145 - (-t139 * t149 * t246 + (-0.2e1 * t139 * t264 + (-t140 * t158 * t257 + t245) * t240 + (t241 * t268 - t140 + (t140 - t224) * t158) * t239 - t217 * t261) * t145) * t222 (-t138 * t262 - t148 * t254) * t246 + (-t138 * t261 + (t165 * t203 + t219 * t247) * t148 + (t138 * t243 - t149 * t254) * t136 + (-t133 * t174 - t144 * t154 + (-t203 * t248 + t210 * t247) * t206 + (-t144 * t184 - t189 * t203) * t140) * t239 + (-t189 * t247 - t133 * t184 - t144 * t169 - t167 * t203 + (t144 * t174 - t203 * t251) * t140) * t240) * t145, 0 (-t137 * t262 - t148 * t180) * t246 + (t137 * t136 * t243 + t153 * t148 + (-t180 * t136 - t137 * t152 + (-t135 * t174 - t141 * t154 + t170 + (-t141 * t184 - t176) * t140) * t259 + (-t135 * t184 - t141 * t169 - t155 + (t141 * t174 - t185) * t140) * t260) * t149) * t145, 0, 0; (-t176 * t186 + t189 * t255) * t244 + (t155 * t186 + t189 * t165 * t242 + (-t189 * t153 - t165 * t176 - t167 * t180) * t187) * t162 (t186 * t204 * t219 + t192 * t255) * t244 + (qJD(4) * t186 * t254 + (-t153 * t192 - t166 * t180) * t187 + (t192 * t242 + (t187 * t219 - t186) * t204) * t165) * t162, 0, t222 * t186 * t244 + (t165 * t187 * t222 + t152 * t186) * t162, 0, 0;];
JaD_rot  = t1;
