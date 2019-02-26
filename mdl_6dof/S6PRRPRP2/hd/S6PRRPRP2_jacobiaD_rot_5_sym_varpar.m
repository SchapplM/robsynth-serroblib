% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRP2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:48
% EndTime: 2019-02-26 20:01:49
% DurationCPUTime: 0.99s
% Computational Cost: add. (5804->110), mult. (9085->227), div. (559->12), fcn. (11668->13), ass. (0->106)
t224 = sin(pkin(10));
t226 = cos(pkin(10));
t229 = sin(qJ(2));
t227 = cos(pkin(6));
t231 = cos(qJ(2));
t257 = t227 * t231;
t213 = -t224 * t229 + t226 * t257;
t209 = t213 * qJD(2);
t258 = t227 * t229;
t214 = t224 * t231 + t226 * t258;
t223 = qJ(3) + pkin(11);
t221 = sin(t223);
t225 = sin(pkin(6));
t261 = t225 * t226;
t247 = t221 * t261;
t222 = cos(t223);
t254 = qJD(3) * t222;
t176 = -qJD(3) * t247 + t209 * t221 + t214 * t254;
t198 = t214 * t221 + t222 * t261;
t196 = t198 ^ 2;
t260 = t225 * t229;
t206 = t221 * t260 - t227 * t222;
t204 = 0.1e1 / t206 ^ 2;
t190 = t196 * t204 + 0.1e1;
t188 = 0.1e1 / t190;
t207 = t227 * t221 + t222 * t260;
t255 = qJD(2) * t231;
t246 = t225 * t255;
t194 = t207 * qJD(3) + t221 * t246;
t203 = 0.1e1 / t206;
t266 = t198 * t204;
t160 = (-t176 * t203 + t194 * t266) * t188;
t191 = atan2(-t198, t206);
t184 = sin(t191);
t185 = cos(t191);
t243 = -t184 * t206 - t185 * t198;
t156 = t243 * t160 - t184 * t176 + t185 * t194;
t170 = -t184 * t198 + t185 * t206;
t167 = 0.1e1 / t170;
t168 = 0.1e1 / t170 ^ 2;
t280 = t156 * t167 * t168;
t248 = t224 * t258;
t216 = t226 * t231 - t248;
t262 = t224 * t225;
t240 = -t216 * t221 + t222 * t262;
t279 = -0.2e1 * t240 * t280;
t259 = t225 * t231;
t239 = -t203 * t213 + t259 * t266;
t278 = t221 * t239;
t267 = t194 * t203 * t204;
t277 = -0.2e1 * (t176 * t266 - t196 * t267) / t190 ^ 2;
t202 = t216 * t222 + t221 * t262;
t230 = cos(qJ(5));
t215 = t224 * t257 + t226 * t229;
t228 = sin(qJ(5));
t264 = t215 * t228;
t187 = t202 * t230 + t264;
t181 = 0.1e1 / t187;
t182 = 0.1e1 / t187 ^ 2;
t211 = t215 * qJD(2);
t179 = t240 * qJD(3) - t211 * t222;
t212 = -qJD(2) * t248 + t226 * t255;
t171 = t187 * qJD(5) + t179 * t228 - t212 * t230;
t263 = t215 * t230;
t186 = t202 * t228 - t263;
t180 = t186 ^ 2;
t175 = t180 * t182 + 0.1e1;
t271 = t182 * t186;
t253 = qJD(5) * t186;
t172 = t179 * t230 + t212 * t228 - t253;
t274 = t172 * t181 * t182;
t276 = (t171 * t271 - t180 * t274) / t175 ^ 2;
t275 = t168 * t240;
t178 = t202 * qJD(3) - t211 * t221;
t273 = t178 * t168;
t272 = t181 * t228;
t270 = t184 * t240;
t269 = t185 * t240;
t268 = t186 * t230;
t265 = t215 * t221;
t256 = qJD(2) * t229;
t197 = t240 ^ 2;
t166 = t197 * t168 + 0.1e1;
t252 = 0.2e1 * (-t197 * t280 - t240 * t273) / t166 ^ 2;
t251 = -0.2e1 * t276;
t249 = t186 * t274;
t245 = -0.2e1 * t198 * t267;
t244 = qJD(5) * t215 * t222 - t211;
t242 = t182 * t268 - t272;
t200 = t214 * t222 - t247;
t241 = -t200 * t203 + t207 * t266;
t238 = qJD(3) * t265 + qJD(5) * t216 - t212 * t222;
t210 = t214 * qJD(2);
t195 = -t206 * qJD(3) + t222 * t246;
t193 = t216 * t228 - t222 * t263;
t192 = -t216 * t230 - t222 * t264;
t177 = -t198 * qJD(3) + t209 * t222;
t173 = 0.1e1 / t175;
t163 = 0.1e1 / t166;
t162 = t188 * t278;
t161 = t241 * t188;
t158 = (-t184 * t213 + t185 * t259) * t221 + t243 * t162;
t157 = t243 * t161 - t184 * t200 + t185 * t207;
t155 = t241 * t277 + (t207 * t245 - t177 * t203 + (t176 * t207 + t194 * t200 + t195 * t198) * t204) * t188;
t153 = t277 * t278 + (t239 * t254 + (t245 * t259 + t203 * t210 + (t194 * t213 + (t176 * t231 - t198 * t256) * t225) * t204) * t221) * t188;
t1 = [0, t153, t155, 0, 0, 0; 0 (-t158 * t275 + t167 * t265) * t252 + ((-t212 * t221 - t215 * t254) * t167 + (-t273 + t279) * t158 + (t265 * t156 + (-t153 * t198 - t162 * t176 + (-t221 * t256 + t231 * t254) * t225 + (-t162 * t206 - t213 * t221) * t160) * t269 + (-t213 * t254 - t153 * t206 - t162 * t194 + t210 * t221 + (t162 * t198 - t221 * t259) * t160) * t270) * t168) * t163 (-t157 * t275 - t167 * t202) * t252 + (t157 * t279 + t179 * t167 + (-t202 * t156 - t157 * t178 + (-t155 * t198 - t161 * t176 + t195 + (-t161 * t206 - t200) * t160) * t269 + (-t155 * t206 - t161 * t194 - t177 + (t161 * t198 - t207) * t160) * t270) * t168) * t163, 0, 0, 0; 0, 0.2e1 * (-t181 * t192 + t193 * t271) * t276 + (0.2e1 * t193 * t249 - t244 * t181 * t230 + t238 * t272 + (-t244 * t186 * t228 - t193 * t171 - t192 * t172 - t238 * t268) * t182) * t173, -t242 * t240 * t251 + (t242 * t178 - ((-qJD(5) * t181 - 0.2e1 * t249) * t230 + (t171 * t230 + (t172 - t253) * t228) * t182) * t240) * t173, 0, t251 + 0.2e1 * (t171 * t182 * t173 + (-t173 * t274 - t182 * t276) * t186) * t186, 0;];
JaD_rot  = t1;
