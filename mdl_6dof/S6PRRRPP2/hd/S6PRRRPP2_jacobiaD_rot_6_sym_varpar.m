% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:09:33
% EndTime: 2019-02-26 20:09:34
% DurationCPUTime: 0.86s
% Computational Cost: add. (3027->111), mult. (9156->228), div. (563->12), fcn. (11754->13), ass. (0->106)
t223 = sin(pkin(10));
t225 = cos(pkin(10));
t229 = sin(qJ(2));
t226 = cos(pkin(6));
t232 = cos(qJ(2));
t259 = t226 * t232;
t214 = -t223 * t229 + t225 * t259;
t207 = t214 * qJD(2);
t260 = t226 * t229;
t215 = t223 * t232 + t225 * t260;
t228 = sin(qJ(3));
t224 = sin(pkin(6));
t263 = t224 * t228;
t248 = t225 * t263;
t231 = cos(qJ(3));
t256 = qJD(3) * t231;
t181 = -qJD(3) * t248 + t207 * t228 + t215 * t256;
t262 = t224 * t231;
t201 = t215 * t228 + t225 * t262;
t199 = t201 ^ 2;
t241 = -t226 * t231 + t229 * t263;
t212 = 0.1e1 / t241 ^ 2;
t195 = t199 * t212 + 0.1e1;
t193 = 0.1e1 / t195;
t219 = -t226 * t228 - t229 * t262;
t257 = qJD(2) * t232;
t247 = t224 * t257;
t205 = t219 * qJD(3) - t228 * t247;
t211 = 0.1e1 / t241;
t268 = t201 * t212;
t165 = (-t181 * t211 - t205 * t268) * t193;
t196 = atan2(t201, -t241);
t191 = sin(t196);
t192 = cos(t196);
t244 = t191 * t241 + t192 * t201;
t161 = t244 * t165 + t191 * t181 + t192 * t205;
t175 = t191 * t201 - t192 * t241;
t172 = 0.1e1 / t175;
t173 = 0.1e1 / t175 ^ 2;
t281 = t161 * t172 * t173;
t249 = t223 * t260;
t217 = t225 * t232 - t249;
t203 = -t217 * t228 + t223 * t262;
t280 = 0.2e1 * t203 * t281;
t266 = t205 * t211 * t212;
t279 = (t181 * t268 + t199 * t266) / t195 ^ 2;
t261 = t224 * t232;
t250 = t201 * t261;
t240 = -t211 * t214 + t212 * t250;
t278 = t228 * t240;
t204 = t217 * t231 + t223 * t263;
t216 = t223 * t259 + t225 * t229;
t227 = sin(qJ(4));
t230 = cos(qJ(4));
t190 = t204 * t230 + t216 * t227;
t186 = 0.1e1 / t190;
t187 = 0.1e1 / t190 ^ 2;
t277 = t173 * t203;
t209 = t216 * qJD(2);
t183 = -t204 * qJD(3) + t209 * t228;
t276 = t183 * t173;
t189 = -t204 * t227 + t216 * t230;
t185 = t189 ^ 2;
t275 = t185 * t187;
t188 = t186 * t187;
t274 = t185 * t188;
t273 = t186 * t227;
t272 = t187 * t189;
t271 = t189 * t230;
t270 = t191 * t203;
t269 = t192 * t203;
t267 = t201 * t219;
t265 = t216 * t228;
t264 = t216 * t231;
t258 = qJD(2) * t229;
t255 = qJD(4) * t189;
t200 = t203 ^ 2;
t171 = t173 * t200 + 0.1e1;
t254 = 0.2e1 * (-t200 * t281 + t203 * t276) / t171 ^ 2;
t184 = t203 * qJD(3) - t209 * t231;
t210 = -qJD(2) * t249 + t225 * t257;
t177 = t184 * t230 + t210 * t227 + t255;
t180 = 0.1e1 + t275;
t176 = -t190 * qJD(4) - t184 * t227 + t210 * t230;
t251 = t176 * t272;
t253 = 0.2e1 * (-t177 * t274 + t251) / t180 ^ 2;
t246 = 0.2e1 * t189 * t188 * t177;
t245 = qJD(4) * t264 - t209;
t243 = t187 * t271 + t273;
t202 = t215 * t231 - t248;
t242 = t202 * t211 + t212 * t267;
t239 = qJD(3) * t265 + qJD(4) * t217 - t210 * t231;
t208 = t215 * qJD(2);
t206 = t241 * qJD(3) - t231 * t247;
t198 = t217 * t227 - t230 * t264;
t197 = t217 * t230 + t227 * t264;
t182 = -t201 * qJD(3) + t207 * t231;
t178 = 0.1e1 / t180;
t168 = 0.1e1 / t171;
t167 = t193 * t278;
t166 = t242 * t193;
t163 = (t191 * t214 - t192 * t261) * t228 + t244 * t167;
t162 = -t244 * t166 + t191 * t202 + t192 * t219;
t160 = 0.2e1 * t242 * t279 + (-0.2e1 * t266 * t267 - t182 * t211 + (-t181 * t219 - t201 * t206 - t202 * t205) * t212) * t193;
t158 = -0.2e1 * t278 * t279 + (t240 * t256 + (0.2e1 * t250 * t266 + t208 * t211 + (-t205 * t214 + (t181 * t232 - t201 * t258) * t224) * t212) * t228) * t193;
t1 = [0, t158, t160, 0, 0, 0; 0 (t163 * t277 - t172 * t265) * t254 + ((t210 * t228 + t216 * t256) * t172 + (-t276 + t280) * t163 + (-t265 * t161 - (t158 * t201 + t167 * t181 + (t228 * t258 - t232 * t256) * t224 + (t167 * t241 + t214 * t228) * t165) * t269 - (t214 * t256 + t158 * t241 - t167 * t205 - t208 * t228 + (-t167 * t201 + t228 * t261) * t165) * t270) * t173) * t168 (t162 * t277 + t172 * t204) * t254 + (t162 * t280 - t184 * t172 + (t204 * t161 - t162 * t183 - (t160 * t201 - t166 * t181 + t206 + (-t166 * t241 + t202) * t165) * t269 - (t160 * t241 + t166 * t205 + t182 + (t166 * t201 - t219) * t165) * t270) * t173) * t168, 0, 0, 0; 0 (-t186 * t197 + t198 * t272) * t253 + (t198 * t246 + t245 * t186 * t230 - t239 * t273 + (-t245 * t189 * t227 - t198 * t176 - t197 * t177 - t239 * t271) * t187) * t178, t243 * t203 * t253 + (-t243 * t183 + ((-qJD(4) * t186 + t246) * t230 + (-t176 * t230 + (t177 + t255) * t227) * t187) * t203) * t178 (t186 * t190 + t275) * t253 + (-0.2e1 * t251 + (t187 * t190 - t186 + 0.2e1 * t274) * t177) * t178, 0, 0;];
JaD_rot  = t1;
