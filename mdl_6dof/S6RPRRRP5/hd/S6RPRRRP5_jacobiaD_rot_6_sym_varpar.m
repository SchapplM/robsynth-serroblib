% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:19
% EndTime: 2019-02-26 21:10:20
% DurationCPUTime: 1.10s
% Computational Cost: add. (8326->125), mult. (8382->272), div. (1515->15), fcn. (10508->9), ass. (0->118)
t197 = pkin(10) + qJ(3) + qJ(4);
t194 = cos(t197);
t203 = sin(qJ(5));
t277 = sin(qJ(1));
t234 = t277 * t203;
t204 = cos(qJ(5));
t205 = cos(qJ(1));
t255 = t205 * t204;
t180 = t194 * t255 + t234;
t174 = 0.1e1 / t180 ^ 2;
t193 = sin(t197);
t189 = t193 ^ 2;
t202 = t205 ^ 2;
t265 = t189 * t202;
t242 = t174 * t265;
t169 = 0.1e1 + t242;
t229 = qJD(1) * t277;
t198 = qJD(3) + qJD(4);
t258 = t198 * t205;
t237 = t193 * t258;
t215 = t194 * t229 + t237;
t228 = t277 * qJD(5);
t256 = t205 * t203;
t159 = (-qJD(5) * t194 + qJD(1)) * t256 + (t228 - t215) * t204;
t173 = 0.1e1 / t180;
t272 = t159 * t173 * t174;
t223 = t265 * t272;
t238 = t193 * t198 * t202;
t280 = (-t223 + (-t189 * t205 * t229 + t194 * t238) * t174) / t169 ^ 2;
t261 = t193 * t205;
t176 = t194 * t234 + t255;
t220 = t203 * t228;
t251 = qJD(5) * t205;
t231 = t204 * t251;
t158 = t176 * qJD(1) - t194 * t231 + t203 * t237 - t220;
t233 = t277 * t204;
t179 = t194 * t256 - t233;
t190 = 0.1e1 / t193;
t191 = 0.1e1 / t193 ^ 2;
t200 = 0.1e1 / t203 ^ 2;
t252 = qJD(5) * t204;
t232 = t200 * t252;
t199 = 0.1e1 / t203;
t259 = t198 * t199;
t236 = t194 * t259;
t264 = t190 * t199;
t279 = (t190 * t232 + t191 * t236) * t179 + t158 * t264;
t262 = t193 * t203;
t168 = atan2(-t176, t262);
t163 = cos(t168);
t162 = sin(t168);
t271 = t162 * t176;
t157 = t163 * t262 - t271;
t154 = 0.1e1 / t157;
t155 = 0.1e1 / t157 ^ 2;
t278 = 0.2e1 * t179;
t171 = t176 ^ 2;
t263 = t191 * t200;
t170 = t171 * t263 + 0.1e1;
t166 = 0.1e1 / t170;
t260 = t194 * t198;
t216 = t193 * t252 + t203 * t260;
t240 = t176 * t263;
t221 = t204 * t229;
t235 = t277 * t193;
t222 = t198 * t235;
t254 = qJD(1) * t205;
t160 = t204 * t228 * t194 - t221 + (t254 * t194 - t222 - t251) * t203;
t243 = t160 * t264;
t146 = (t216 * t240 - t243) * t166;
t213 = -t146 * t176 + t216;
t142 = (-t146 * t262 - t160) * t162 + t213 * t163;
t156 = t154 * t155;
t276 = t142 * t156;
t192 = t190 / t189;
t201 = t199 * t200;
t275 = (t160 * t240 + (-t191 * t201 * t252 - t192 * t200 * t260) * t171) / t170 ^ 2;
t274 = t155 * t179;
t273 = t158 * t155;
t270 = t162 * t179;
t269 = t162 * t193;
t268 = t163 * t176;
t267 = t163 * t179;
t266 = t163 * t194;
t257 = t200 * t204;
t253 = qJD(5) * t203;
t172 = t179 ^ 2;
t152 = t155 * t172 + 0.1e1;
t250 = 0.2e1 * (-t172 * t276 - t179 * t273) / t152 ^ 2;
t249 = 0.2e1 * t280;
t248 = -0.2e1 * t275;
t247 = t156 * t278;
t246 = t190 * t275;
t245 = t155 * t270;
t241 = t176 * t264;
t239 = t191 * t194 * t199;
t218 = t176 * t239 + t277;
t153 = t218 * t166;
t230 = t277 - t153;
t227 = t154 * t250;
t226 = t155 * t250;
t225 = t261 * t278;
t224 = t199 * t246;
t178 = t194 * t233 - t256;
t219 = t176 * t257 - t178 * t199;
t217 = t174 * t178 * t205 - t277 * t173;
t164 = 0.1e1 / t169;
t161 = t180 * qJD(1) - t194 * t220 - t204 * t222 - t231;
t150 = 0.1e1 / t152;
t149 = t219 * t190 * t166;
t145 = (-t162 + (t163 * t241 + t162) * t166) * t179;
t144 = -t153 * t268 + (t230 * t269 + t266) * t203;
t143 = t163 * t193 * t204 - t162 * t178 + (-t162 * t262 - t268) * t149;
t141 = t218 * t248 + (t160 * t239 + t254 + (-t190 * t259 + (-t191 * t232 - 0.2e1 * t192 * t236) * t194) * t176) * t166;
t139 = -0.2e1 * t219 * t246 + (-t219 * t191 * t260 + (t160 * t257 - t161 * t199 + (t178 * t257 + (-0.2e1 * t201 * t204 ^ 2 - t199) * t176) * qJD(5)) * t190) * t166;
t138 = (t173 * t194 * t205 + t204 * t242) * t249 + (0.2e1 * t204 * t223 + t215 * t173 + ((t159 * t205 - 0.2e1 * t204 * t238) * t194 + (t202 * t253 + 0.2e1 * t205 * t221) * t189) * t174) * t164;
t137 = t144 * t179 * t226 + (-(-t141 * t268 + (t146 * t271 - t160 * t163) * t153) * t274 + (t142 * t247 + t273) * t144 + (-t154 * t261 - (-t153 * t269 + t162 * t235 + t266) * t274) * t252) * t150 + (t227 * t261 + ((-t154 * t258 - (t230 * t198 - t146) * t245) * t194 + (t154 * t229 + (t205 * t142 - (-t141 + t254) * t270 - (t230 * t146 - t198) * t267) * t155) * t193) * t150) * t203;
t1 = [t279 * t166 + t224 * t278, 0, t141, t141, t139, 0; t176 * t227 + (-t160 * t154 + (t142 * t176 + t145 * t158) * t155) * t150 + (t145 * t226 + (0.2e1 * t145 * t276 + (t158 * t166 - t158 - (-t146 * t166 * t241 + t248) * t179) * t155 * t162 + (-(-0.2e1 * t176 * t224 - t146) * t274 + (-(t146 + t243) * t179 + t279 * t176) * t155 * t166) * t163) * t150) * t179, 0, t137, t137 (t143 * t274 - t154 * t180) * t250 + (t143 * t273 + t159 * t154 + (t143 * t247 - t155 * t180) * t142 - (-t193 * t253 + t204 * t260 - t139 * t176 - t149 * t160 + (-t149 * t262 - t178) * t146) * t155 * t267 - (-t161 + (-t139 * t203 - t146 * t204) * t193 - t213 * t149) * t245) * t150, 0; t217 * t193 * t249 + (-t217 * t260 + ((qJD(1) * t173 + 0.2e1 * t178 * t272) * t205 + (-t277 * t159 - t161 * t205 + t178 * t229) * t174) * t193) * t164, 0, t138, t138, t174 * t225 * t280 + (t225 * t272 + (t158 * t261 + (t193 * t229 - t194 * t258) * t179) * t174) * t164, 0;];
JaD_rot  = t1;
