% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR5
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
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:32:53
% EndTime: 2019-02-26 22:32:54
% DurationCPUTime: 1.09s
% Computational Cost: add. (5529->125), mult. (8382->273), div. (1515->15), fcn. (10508->9), ass. (0->119)
t201 = qJ(2) + qJ(3);
t195 = cos(t201);
t202 = sin(qJ(4));
t277 = sin(qJ(1));
t233 = t277 * t202;
t203 = cos(qJ(4));
t204 = cos(qJ(1));
t254 = t204 * t203;
t179 = t195 * t254 + t233;
t173 = 0.1e1 / t179 ^ 2;
t194 = sin(t201);
t190 = t194 ^ 2;
t200 = t204 ^ 2;
t265 = t190 * t200;
t235 = t173 * t265;
t168 = 0.1e1 + t235;
t228 = qJD(1) * t277;
t196 = qJD(2) + qJD(3);
t258 = t196 * t204;
t237 = t194 * t258;
t214 = t195 * t228 + t237;
t227 = t277 * qJD(4);
t255 = t204 * t202;
t158 = (-qJD(4) * t195 + qJD(1)) * t255 + (t227 - t214) * t203;
t172 = 0.1e1 / t179;
t272 = t158 * t172 * t173;
t222 = t265 * t272;
t238 = t194 * t196 * t200;
t280 = (-t222 + (-t190 * t204 * t228 + t195 * t238) * t173) / t168 ^ 2;
t261 = t194 * t204;
t175 = t195 * t233 + t254;
t219 = t202 * t227;
t250 = qJD(4) * t204;
t230 = t203 * t250;
t157 = t175 * qJD(1) - t195 * t230 + t202 * t237 - t219;
t232 = t277 * t203;
t178 = t195 * t255 - t232;
t191 = 0.1e1 / t194;
t192 = 0.1e1 / t194 ^ 2;
t198 = 0.1e1 / t202 ^ 2;
t251 = qJD(4) * t203;
t231 = t198 * t251;
t197 = 0.1e1 / t202;
t259 = t196 * t197;
t236 = t195 * t259;
t264 = t191 * t197;
t279 = (t191 * t231 + t192 * t236) * t178 + t157 * t264;
t262 = t194 * t202;
t167 = atan2(-t175, t262);
t162 = cos(t167);
t161 = sin(t167);
t271 = t161 * t175;
t156 = t162 * t262 - t271;
t153 = 0.1e1 / t156;
t154 = 0.1e1 / t156 ^ 2;
t278 = 0.2e1 * t178;
t170 = t175 ^ 2;
t263 = t192 * t198;
t169 = t170 * t263 + 0.1e1;
t165 = 0.1e1 / t169;
t260 = t195 * t196;
t215 = t194 * t251 + t202 * t260;
t240 = t175 * t263;
t220 = t203 * t228;
t234 = t277 * t194;
t221 = t196 * t234;
t253 = qJD(1) * t204;
t159 = t203 * t227 * t195 - t220 + (t253 * t195 - t221 - t250) * t202;
t242 = t159 * t264;
t145 = (t215 * t240 - t242) * t165;
t212 = -t145 * t175 + t215;
t141 = (-t145 * t262 - t159) * t161 + t212 * t162;
t155 = t153 * t154;
t276 = t141 * t155;
t193 = t191 / t190;
t199 = t197 * t198;
t275 = (t159 * t240 + (-t192 * t199 * t251 - t193 * t198 * t260) * t170) / t169 ^ 2;
t274 = t154 * t178;
t273 = t157 * t154;
t270 = t161 * t178;
t269 = t161 * t194;
t268 = t162 * t175;
t267 = t162 * t178;
t266 = t162 * t195;
t257 = t198 * t203;
t256 = t204 * t153;
t252 = qJD(4) * t202;
t171 = t178 ^ 2;
t151 = t171 * t154 + 0.1e1;
t249 = 0.2e1 * (-t171 * t276 - t178 * t273) / t151 ^ 2;
t248 = 0.2e1 * t280;
t247 = -0.2e1 * t275;
t246 = t155 * t278;
t245 = t191 * t275;
t244 = t154 * t270;
t241 = t175 * t264;
t239 = t192 * t195 * t197;
t217 = t175 * t239 + t277;
t152 = t217 * t165;
t229 = t277 - t152;
t226 = t153 * t249;
t225 = t154 * t249;
t224 = t261 * t278;
t223 = t197 * t245;
t177 = t195 * t232 - t255;
t218 = t175 * t257 - t177 * t197;
t216 = t173 * t177 * t204 - t277 * t172;
t163 = 0.1e1 / t168;
t160 = t179 * qJD(1) - t195 * t219 - t203 * t221 - t230;
t149 = 0.1e1 / t151;
t148 = t218 * t191 * t165;
t144 = (-t161 + (t162 * t241 + t161) * t165) * t178;
t143 = -t152 * t268 + (t229 * t269 + t266) * t202;
t142 = t162 * t194 * t203 - t161 * t177 + (-t161 * t262 - t268) * t148;
t140 = t217 * t247 + (t159 * t239 + t253 + (-t191 * t259 + (-t192 * t231 - 0.2e1 * t193 * t236) * t195) * t175) * t165;
t138 = -0.2e1 * t218 * t245 + (-t218 * t192 * t260 + (t159 * t257 - t160 * t197 + (t177 * t257 + (-0.2e1 * t199 * t203 ^ 2 - t197) * t175) * qJD(4)) * t191) * t165;
t137 = (t172 * t195 * t204 + t203 * t235) * t248 + (0.2e1 * t203 * t222 + t214 * t172 + ((t158 * t204 - 0.2e1 * t203 * t238) * t195 + (t200 * t252 + 0.2e1 * t204 * t220) * t190) * t173) * t163;
t136 = t143 * t178 * t225 + (-(-t140 * t268 + (t145 * t271 - t159 * t162) * t152) * t274 + (t141 * t246 + t273) * t143 + (-t194 * t256 - (-t152 * t269 + t161 * t234 + t266) * t274) * t251) * t149 + (t226 * t261 + ((-t196 * t256 - (t229 * t196 - t145) * t244) * t195 + (t153 * t228 + (t204 * t141 - (-t140 + t253) * t270 - (t229 * t145 - t196) * t267) * t154) * t194) * t149) * t202;
t1 = [t279 * t165 + t223 * t278, t140, t140, t138, 0, 0; t175 * t226 + (-t159 * t153 + (t141 * t175 + t144 * t157) * t154) * t149 + (t144 * t225 + (0.2e1 * t144 * t276 + (t157 * t165 - t157 - (-t145 * t165 * t241 + t247) * t178) * t154 * t161 + (-(-0.2e1 * t175 * t223 - t145) * t274 + (-(t145 + t242) * t178 + t279 * t175) * t154 * t165) * t162) * t149) * t178, t136, t136 (t142 * t274 - t153 * t179) * t249 + (t142 * t273 + t158 * t153 + (t142 * t246 - t179 * t154) * t141 - (-t194 * t252 + t203 * t260 - t138 * t175 - t148 * t159 + (-t148 * t262 - t177) * t145) * t154 * t267 - (-t160 + (-t138 * t202 - t145 * t203) * t194 - t212 * t148) * t244) * t149, 0, 0; t216 * t194 * t248 + (-t216 * t260 + ((qJD(1) * t172 + 0.2e1 * t177 * t272) * t204 + (-t277 * t158 - t160 * t204 + t177 * t228) * t173) * t194) * t163, t137, t137, t173 * t224 * t280 + (t224 * t272 + (t157 * t261 + (t194 * t228 - t195 * t258) * t178) * t173) * t163, 0, 0;];
JaD_rot  = t1;
