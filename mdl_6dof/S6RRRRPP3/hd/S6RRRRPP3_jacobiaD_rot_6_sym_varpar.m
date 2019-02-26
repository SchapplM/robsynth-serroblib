% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:21
% EndTime: 2019-02-26 22:26:23
% DurationCPUTime: 1.09s
% Computational Cost: add. (5428->121), mult. (8127->261), div. (1567->14), fcn. (10302->9), ass. (0->116)
t198 = qJ(2) + qJ(3);
t189 = sin(t198);
t202 = cos(qJ(1));
t276 = t189 * t202;
t185 = 0.1e1 / t189;
t186 = 0.1e1 / t189 ^ 2;
t187 = t185 * t186;
t190 = cos(t198);
t191 = qJD(2) + qJD(3);
t275 = t191 * (0.2e1 * t187 * t190 ^ 2 + t185);
t199 = sin(qJ(4));
t200 = sin(qJ(1));
t201 = cos(qJ(4));
t227 = t191 * t276;
t243 = qJD(4) * t202;
t244 = qJD(4) * t201;
t246 = qJD(1) * t202;
t247 = qJD(1) * t200;
t159 = t201 * t227 - t200 * t244 - t199 * t246 + (t199 * t243 + t201 * t247) * t190;
t251 = t201 * t202;
t253 = t200 * t199;
t178 = t190 * t251 + t253;
t193 = 0.1e1 / t201 ^ 2;
t245 = qJD(4) * t199;
t224 = t193 * t245;
t192 = 0.1e1 / t201;
t261 = t186 * t190;
t231 = t192 * t261;
t262 = t185 * t192;
t274 = -(t185 * t224 - t191 * t231) * t178 + t159 * t262;
t249 = t202 * t199;
t252 = t200 * t201;
t175 = t190 * t252 - t249;
t257 = t189 * t201;
t168 = atan2(-t175, t257);
t163 = cos(t168);
t162 = sin(t168);
t268 = t162 * t175;
t157 = t163 * t257 - t268;
t154 = 0.1e1 / t157;
t195 = 0.1e1 / t202;
t155 = 0.1e1 / t157 ^ 2;
t196 = 0.1e1 / t202 ^ 2;
t273 = 0.2e1 * t178;
t171 = t175 ^ 2;
t260 = t186 * t193;
t169 = t171 * t260 + 0.1e1;
t164 = 0.1e1 / t169;
t256 = t190 * t191;
t213 = -t189 * t245 + t201 * t256;
t233 = t175 * t260;
t226 = t190 * t253;
t258 = t189 * t200;
t228 = t191 * t258;
t161 = t178 * qJD(1) - qJD(4) * t226 + (-t228 - t243) * t201;
t235 = t161 * t262;
t146 = (t213 * t233 - t235) * t164;
t211 = t146 * t175 - t213;
t142 = (-t146 * t257 - t161) * t162 - t211 * t163;
t156 = t154 * t155;
t272 = t142 * t156;
t194 = t192 * t193;
t229 = t187 * t256;
t271 = (t161 * t233 + (t186 * t194 * t245 - t193 * t229) * t171) / t169 ^ 2;
t270 = t155 * t178;
t269 = t159 * t155;
t267 = t162 * t178;
t266 = t162 * t189;
t265 = t163 * t175;
t264 = t163 * t178;
t263 = t163 * t190;
t259 = t186 * t196;
t255 = t193 * t199;
t254 = t196 * t200;
t250 = t202 * t154;
t216 = t175 * t231 + t200;
t153 = t216 * t164;
t248 = -t153 + t200;
t173 = t178 ^ 2;
t152 = t155 * t173 + 0.1e1;
t242 = 0.2e1 * (-t173 * t272 - t178 * t269) / t152 ^ 2;
t241 = -0.2e1 * t271;
t218 = qJD(1) * t190 - qJD(4);
t219 = qJD(4) * t190 - qJD(1);
t158 = -t219 * t251 + (t218 * t200 + t227) * t199;
t177 = -t190 * t249 + t252;
t172 = t177 ^ 2;
t170 = t172 * t259 + 0.1e1;
t197 = t195 * t196;
t240 = 0.2e1 * (t158 * t177 * t259 + (t186 * t197 * t247 - t196 * t229) * t172) / t170 ^ 2;
t239 = t156 * t273;
t238 = t185 * t271;
t237 = t155 * t267;
t234 = t175 * t262;
t232 = t186 * t256;
t230 = t195 * t261;
t225 = t196 * t247;
t223 = t154 * t242;
t222 = t155 * t242;
t221 = t185 * t240;
t217 = t192 * t238;
t174 = t226 + t251;
t215 = -t174 * t192 + t175 * t255;
t214 = -t174 * t195 - t177 * t254;
t166 = 0.1e1 / t170;
t160 = t219 * t252 + (t218 * t202 - t228) * t199;
t150 = 0.1e1 / t152;
t149 = t215 * t185 * t164;
t145 = (-t162 + (t163 * t234 + t162) * t164) * t178;
t144 = -t153 * t265 + (t248 * t266 + t263) * t201;
t143 = -t163 * t189 * t199 + t162 * t174 - (-t162 * t257 - t265) * t149;
t141 = t216 * t241 + (t161 * t231 + t246 + (-t192 * t275 + t224 * t261) * t175) * t164;
t140 = (t177 * t230 - t199) * t240 + (-t158 * t230 + t244 + (t195 * t275 - t225 * t261) * t177) * t166;
t138 = 0.2e1 * t215 * t238 + (t215 * t232 + (-t161 * t255 + t160 * t192 + (t174 * t255 + (-0.2e1 * t194 * t199 ^ 2 - t192) * t175) * qJD(4)) * t185) * t164;
t137 = t144 * t178 * t222 + (-(-t141 * t265 + (t146 * t268 - t161 * t163) * t153) * t270 + (t142 * t239 + t269) * t144 + (t189 * t250 - (t153 * t266 - t162 * t258 - t263) * t270) * t245) * t150 + (t223 * t276 + ((-t191 * t250 - (t248 * t191 - t146) * t237) * t190 + (t154 * t247 + (t202 * t142 - (-t141 + t246) * t267 - (t248 * t146 - t191) * t264) * t155) * t189) * t150) * t201;
t1 = [t164 * t274 + t217 * t273, t141, t141, t138, 0, 0; t175 * t223 + (-t161 * t154 + (t142 * t175 + t145 * t159) * t155) * t150 + (t145 * t222 + (0.2e1 * t145 * t272 + (t159 * t164 - t159 - (-t146 * t164 * t234 + t241) * t178) * t155 * t162 + (-(-0.2e1 * t175 * t217 - t146) * t270 + (-(t146 + t235) * t178 + t274 * t175) * t155 * t164) * t163) * t150) * t178, t137, t137 (t143 * t270 - t154 * t177) * t242 + (t143 * t269 + t158 * t154 + (t143 * t239 - t155 * t177) * t142 - (-t189 * t244 - t199 * t256 - t138 * t175 + t149 * t161 + (t149 * t257 + t174) * t146) * t155 * t264 - (t160 + (-t138 * t201 + t146 * t199) * t189 - t211 * t149) * t237) * t150, 0, 0; t214 * t221 + (t214 * t232 + (t158 * t254 + t160 * t195 + (t174 * t254 + (0.2e1 * t197 * t200 ^ 2 + t195) * t177) * qJD(1)) * t185) * t166, t140, t140, t178 * t195 * t221 + (t159 * t185 * t195 + (-t185 * t225 + t191 * t230) * t178) * t166, 0, 0;];
JaD_rot  = t1;
