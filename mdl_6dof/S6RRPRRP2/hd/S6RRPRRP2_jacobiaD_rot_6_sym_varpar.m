% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:46:35
% EndTime: 2019-02-26 21:46:36
% DurationCPUTime: 1.12s
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
t1 = [t282 * t168 + t226 * t281, t143, 0, t143, t141, 0; t178 * t229 + (-t162 * t156 + (t144 * t178 + t147 * t160) * t157) * t152 + (t147 * t228 + (0.2e1 * t147 * t279 + (t160 * t168 - t160 - (-t148 * t168 * t243 + t250) * t181) * t157 * t164 + (-(-0.2e1 * t178 * t226 - t148) * t277 + (-(t148 + t245) * t181 + t282 * t178) * t157 * t168) * t165) * t152) * t181, t139, 0, t139 (t145 * t277 - t156 * t182) * t252 + (t145 * t276 + t161 * t156 + (t145 * t249 - t157 * t182) * t144 - (-t195 * t255 + t206 * t263 - t141 * t178 - t151 * t162 + (-t151 * t265 - t180) * t148) * t157 * t270 - (-t163 + (-t141 * t205 - t148 * t206) * t195 - t215 * t151) * t247) * t152, 0; t219 * t195 * t251 + (-t219 * t263 + ((qJD(1) * t175 + 0.2e1 * t180 * t275) * t207 + (-t280 * t161 - t163 * t207 + t180 * t231) * t176) * t195) * t166, t140, 0, t140, t176 * t227 * t283 + (t227 * t275 + (t160 * t264 + (t195 * t231 - t196 * t261) * t181) * t176) * t166, 0;];
JaD_rot  = t1;
