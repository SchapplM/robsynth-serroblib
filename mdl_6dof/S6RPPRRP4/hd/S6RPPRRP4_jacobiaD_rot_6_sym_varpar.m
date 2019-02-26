% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:05
% EndTime: 2019-02-26 20:32:06
% DurationCPUTime: 1.27s
% Computational Cost: add. (4280->124), mult. (10807->281), div. (1114->15), fcn. (14382->11), ass. (0->112)
t265 = sin(pkin(9));
t266 = cos(pkin(9));
t267 = sin(qJ(1));
t268 = cos(qJ(1));
t180 = -t267 * t265 - t268 * t266;
t181 = t268 * t265 - t267 * t266;
t196 = cos(qJ(5));
t194 = sin(qJ(5));
t197 = cos(qJ(4));
t244 = t194 * t197;
t164 = -t180 * t196 + t181 * t244;
t195 = sin(qJ(4));
t243 = t195 * t194;
t155 = atan2(t164, -t243);
t153 = sin(t155);
t154 = cos(t155);
t159 = t164 ^ 2;
t188 = 0.1e1 / t194 ^ 2;
t192 = 0.1e1 / t195 ^ 2;
t246 = t188 * t192;
t158 = t159 * t246 + 0.1e1;
t156 = 0.1e1 / t158;
t187 = 0.1e1 / t194;
t225 = t187 * t192 * t197;
t212 = -t164 * t225 - t181;
t141 = t212 * t156;
t241 = -t141 - t181;
t272 = t241 * t153 * t195 - t154 * t197;
t248 = t180 * t195;
t242 = t196 * t197;
t170 = -t180 * t242 + t181 * t194;
t177 = t180 * qJD(1);
t178 = t181 * qJD(1);
t240 = qJD(4) * t195;
t223 = t194 * t240;
t148 = t170 * qJD(5) - t177 * t196 + t178 * t244 + t180 * t223;
t169 = -t180 * t244 - t181 * t196;
t191 = 0.1e1 / t195;
t239 = qJD(4) * t197;
t224 = t192 * t239;
t237 = qJD(5) * t196;
t247 = t187 * t191;
t270 = -(t188 * t191 * t237 + t187 * t224) * t169 + t148 * t247;
t257 = t153 * t164;
t145 = -t154 * t243 + t257;
t142 = 0.1e1 / t145;
t161 = 0.1e1 / t170;
t143 = 0.1e1 / t145 ^ 2;
t162 = 0.1e1 / t170 ^ 2;
t269 = -0.2e1 * t169;
t160 = t169 ^ 2;
t140 = t143 * t160 + 0.1e1;
t259 = t148 * t143;
t208 = -t194 * t239 - t195 * t237;
t226 = t164 * t246;
t165 = t180 * t194 + t181 * t242;
t250 = t178 * t196;
t146 = -t165 * qJD(5) - t177 * t244 + t181 * t223 - t250;
t230 = t146 * t247;
t135 = (-t208 * t226 + t230) * t156;
t205 = -t135 * t164 - t208;
t130 = (t135 * t243 - t146) * t153 - t205 * t154;
t144 = t142 * t143;
t263 = t130 * t144;
t264 = (-t160 * t263 + t169 * t259) / t140 ^ 2;
t179 = t180 ^ 2;
t190 = t195 ^ 2;
t249 = t179 * t190;
t228 = t162 * t249;
t152 = 0.1e1 + t228;
t209 = -t178 * t197 - t180 * t240;
t236 = qJD(5) * t197;
t149 = (t180 * t236 + t177) * t194 + (qJD(5) * t181 - t209) * t196;
t258 = t149 * t161 * t162;
t216 = t249 * t258;
t262 = (-t216 + (-t178 * t180 * t190 + t179 * t195 * t239) * t162) / t152 ^ 2;
t189 = t187 * t188;
t193 = t191 / t190;
t221 = t192 * t237;
t261 = (-t146 * t226 + (-t188 * t193 * t239 - t189 * t221) * t159) / t158 ^ 2;
t260 = t143 * t169;
t256 = t153 * t169;
t254 = t154 * t164;
t253 = t154 * t169;
t251 = t162 * t180;
t245 = t188 * t196;
t238 = qJD(5) * t194;
t235 = 0.2e1 * t264;
t234 = 0.2e1 * t144 * t169;
t233 = t195 * t262;
t232 = t191 * t261;
t231 = t143 * t256;
t227 = t164 * t247;
t222 = t196 * t239;
t220 = -0.2e1 * t142 * t264;
t219 = t143 * t235;
t218 = t187 * t232;
t217 = t248 * t258;
t213 = t181 * t236 + t178;
t211 = -t164 * t245 + t165 * t187;
t210 = -t178 * t195 + t180 * t239;
t207 = -qJD(5) * t180 - t177 * t197 + t181 * t240;
t147 = t213 * t194 + t207 * t196;
t150 = 0.1e1 / t152;
t138 = 0.1e1 / t140;
t137 = t211 * t191 * t156;
t133 = (-t153 + (t154 * t227 + t153) * t156) * t169;
t132 = -t141 * t254 + t272 * t194;
t131 = -t154 * t195 * t196 + t153 * t165 - (t153 * t243 + t254) * t137;
t129 = 0.2e1 * t212 * t261 + (-t146 * t225 + t177 - (t188 * t197 * t221 + (0.2e1 * t193 * t197 ^ 2 + t191) * t187 * qJD(4)) * t164) * t156;
t127 = 0.2e1 * t211 * t232 + (t211 * t224 + (-t146 * t245 + t147 * t187 + (t165 * t245 - (0.2e1 * t189 * t196 ^ 2 + t187) * t164) * qJD(5)) * t191) * t156;
t1 = [t156 * t270 + t218 * t269, 0, 0, t129, t127, 0; t164 * t220 + ((-t207 * t194 + t213 * t196) * t142 + (-t130 * t164 - t133 * t148) * t143) * t138 + (t133 * t219 + (0.2e1 * t133 * t263 + (-t148 * t156 + t148 - (-t135 * t156 * t227 - 0.2e1 * t261) * t169) * t143 * t153 + (-(-0.2e1 * t164 * t218 - t135) * t260 + (-(t135 - t230) * t169 - t270 * t164) * t143 * t156) * t154) * t138) * t169, 0, 0, t132 * t169 * t219 + (-(t129 * t254 - (-t135 * t257 - t146 * t154) * t141) * t260 + (t130 * t234 - t259) * t132 + (t142 * t248 - t272 * t260) * t237) * t138 + (t220 * t248 + ((t180 * qJD(4) * t142 - (t241 * qJD(4) + t135) * t231) * t197 + (-t178 * t142 + (-t180 * t130 - (t129 - t177) * t256 - (t241 * t135 + qJD(4)) * t253) * t143) * t195) * t138) * t194 (t131 * t260 - t142 * t170) * t235 + (-t131 * t259 + t149 * t142 + (t131 * t234 - t143 * t170) * t130 - (-t222 + t195 * t238 + t127 * t164 + t137 * t146 + (-t137 * t243 + t165) * t135) * t143 * t253 - (-t147 + (t127 * t194 + t135 * t196) * t195 - t205 * t137) * t231) * t138, 0; 0.2e1 * (t161 * t181 + t165 * t251) * t233 + (0.2e1 * t165 * t217 + (-t177 * t195 - t181 * t239) * t161 + ((t147 * t180 + t181 * t149) * t195 - t210 * t165) * t162) * t150, 0, 0, 0.2e1 * (-t161 * t180 * t197 + t196 * t228) * t262 + (0.2e1 * t196 * t216 + t209 * t161 + ((-t149 * t197 + 0.2e1 * t190 * t250) * t180 + (t190 * t238 - 0.2e1 * t195 * t222) * t179) * t162) * t150, t233 * t251 * t269 + (t217 * t269 + (t148 * t248 + t210 * t169) * t162) * t150, 0;];
JaD_rot  = t1;
