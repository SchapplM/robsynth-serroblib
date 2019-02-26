% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR9_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiaD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:24
% EndTime: 2019-02-26 20:53:25
% DurationCPUTime: 0.61s
% Computational Cost: add. (1382->69), mult. (4220->158), div. (139->12), fcn. (5452->15), ass. (0->87)
t217 = sin(pkin(13));
t221 = cos(pkin(13));
t225 = sin(qJ(3));
t227 = cos(qJ(3));
t209 = t225 * t217 - t227 * t221;
t267 = t209 * qJD(3);
t219 = sin(pkin(7));
t191 = t209 * t219;
t223 = cos(pkin(7));
t193 = t209 * t223;
t224 = cos(pkin(6));
t222 = cos(pkin(12));
t228 = cos(qJ(1));
t248 = t228 * t222;
t218 = sin(pkin(12));
t226 = sin(qJ(1));
t252 = t226 * t218;
t237 = t224 * t252 - t248;
t249 = t228 * t218;
t251 = t226 * t222;
t238 = t224 * t251 + t249;
t240 = t227 * t217 + t225 * t221;
t220 = sin(pkin(6));
t255 = t220 * t226;
t171 = -t191 * t255 + t193 * t238 + t237 * t240;
t164 = t171 ^ 2;
t192 = t240 * t219;
t194 = t240 * t223;
t235 = t192 * t255 - t194 * t238 + t209 * t237;
t166 = 0.1e1 / t235 ^ 2;
t266 = t164 * t166;
t202 = -t220 * t222 * t219 + t224 * t223;
t203 = -t224 * t248 + t252;
t254 = t220 * t228;
t239 = -t203 * t219 + t223 * t254;
t179 = atan2(t239, t202);
t174 = sin(t179);
t175 = cos(t179);
t163 = t174 * t239 + t175 * t202;
t160 = 0.1e1 / t163;
t165 = 0.1e1 / t235;
t264 = t239 ^ 2;
t199 = 0.1e1 / t202;
t161 = 0.1e1 / t163 ^ 2;
t200 = 0.1e1 / t202 ^ 2;
t263 = -0.2e1 * t199 * t200;
t197 = t238 * qJD(1);
t243 = qJD(1) * t220 * t223;
t181 = -t197 * t219 - t226 * t243;
t178 = t200 * t264 + 0.1e1;
t176 = 0.1e1 / t178;
t236 = t174 + (t175 * t199 * t239 - t174) * t176;
t150 = t236 * t181;
t262 = t150 * t160 * t161;
t187 = t267 * t219;
t189 = t267 * t223;
t195 = t203 * qJD(1);
t204 = -t224 * t249 - t251;
t196 = t204 * qJD(1);
t208 = t240 * qJD(3);
t246 = qJD(1) * t228;
t153 = t238 * t189 + t195 * t194 - t196 * t209 + t237 * t208 + (-t187 * t226 + t192 * t246) * t220;
t167 = t165 * t166;
t261 = t167 * t153;
t169 = t192 * t254 + t203 * t194 - t204 * t209;
t260 = t169 * t171;
t259 = t176 * t199;
t177 = 0.1e1 / t178 ^ 2;
t258 = t177 * t239;
t180 = t195 * t219 - t228 * t243;
t257 = t180 * t161;
t185 = -t219 * t238 - t223 * t255;
t256 = t181 * t185;
t247 = qJD(1) * t226;
t156 = 0.1e1 + t266;
t188 = t219 * t208;
t190 = t223 * t208;
t152 = t238 * t190 - t195 * t193 - t196 * t240 - t237 * t267 + (-t188 * t226 - t191 * t246) * t220;
t244 = t171 * t166 * t152;
t245 = 0.2e1 * (-t164 * t261 + t244) / t156 ^ 2;
t198 = t237 * qJD(1);
t182 = t185 ^ 2;
t168 = t191 * t254 + t203 * t193 + t204 * t240;
t159 = t182 * t161 + 0.1e1;
t154 = 0.1e1 / t156;
t151 = t236 * t185;
t1 = [t256 * t258 * t263 + t180 * t259, 0, 0, 0, 0, 0; 0.2e1 * (-t151 * t161 * t185 - t160 * t239) / t159 ^ 2 * (-t182 * t262 + t185 * t257) + (t181 * t160 + (-t150 * t239 + t151 * t180) * t161 + (-0.2e1 * t151 * t262 + t236 * t257 + (t174 * t200 * t258 + (0.2e1 * t259 + (t263 * t264 - t199) * t177) * t175) * t161 * t256) * t185) / t159, 0, 0, 0, 0, 0; (-t165 * t168 - t166 * t260) * t245 + ((t203 * t190 + t197 * t193 + t198 * t240 - t204 * t267 + (t188 * t228 - t191 * t247) * t220) * t165 - 0.2e1 * t260 * t261 + (-t168 * t153 + (-t203 * t189 + t197 * t194 - t198 * t209 - t204 * t208 + (-t187 * t228 - t192 * t247) * t220) * t171 + t169 * t152) * t166) * t154, 0 (-t165 * t235 - t266) * t245 + (0.2e1 * t244 + (-0.2e1 * t164 * t167 - t166 * t235 + t165) * t153) * t154, 0, 0, 0;];
JaD_rot  = t1;
