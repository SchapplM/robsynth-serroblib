% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:18
% EndTime: 2019-02-26 22:19:19
% DurationCPUTime: 0.70s
% Computational Cost: add. (2443->93), mult. (4665->202), div. (686->14), fcn. (5929->11), ass. (0->95)
t205 = cos(qJ(2));
t206 = cos(qJ(1));
t256 = cos(pkin(6));
t227 = t206 * t256;
t225 = t205 * t227;
t203 = sin(qJ(2));
t204 = sin(qJ(1));
t241 = t204 * t203;
t182 = -t225 + t241;
t202 = sin(pkin(6));
t243 = t202 * t205;
t176 = atan2(-t182, -t243);
t174 = sin(t176);
t175 = cos(t176);
t180 = t182 ^ 2;
t197 = 0.1e1 / t202 ^ 2;
t200 = 0.1e1 / t205 ^ 2;
t179 = t180 * t197 * t200 + 0.1e1;
t177 = 0.1e1 / t179;
t196 = 0.1e1 / t202;
t199 = 0.1e1 / t205;
t230 = t182 * t196 * t199;
t257 = (t175 * t230 - t174) * t177 + t174;
t158 = -t174 * t182 - t175 * t243;
t155 = 0.1e1 / t158;
t228 = t204 * t256;
t226 = t203 * t228;
t240 = t206 * t205;
t186 = -t226 + t240;
t195 = qJ(3) + pkin(12) + qJ(5);
t193 = sin(t195);
t194 = cos(t195);
t244 = t202 * t204;
t173 = t186 * t194 + t193 * t244;
t163 = 0.1e1 / t173;
t156 = 0.1e1 / t158 ^ 2;
t164 = 0.1e1 / t173 ^ 2;
t213 = -t203 * t227 - t204 * t205;
t214 = -t206 * t203 - t205 * t228;
t168 = -t214 * qJD(1) - t213 * qJD(2);
t238 = qJD(2) * t203;
t229 = t200 * t238;
t215 = t168 * t199 + t182 * t229;
t246 = t177 * t196;
t147 = t215 * t246;
t219 = t174 * t243 - t175 * t182;
t231 = t175 * t202 * t203;
t143 = qJD(2) * t231 + t219 * t147 - t174 * t168;
t255 = t143 * t155 * t156;
t198 = qJD(3) + qJD(5);
t239 = qJD(1) * t202;
t216 = -t186 * t198 + t206 * t239;
t167 = t213 * qJD(1) + t214 * qJD(2);
t224 = t198 * t244 + t167;
t148 = t224 * t193 - t216 * t194;
t172 = t186 * t193 - t194 * t244;
t162 = t172 ^ 2;
t161 = t162 * t164 + 0.1e1;
t250 = t164 * t172;
t149 = t216 * t193 + t224 * t194;
t252 = t149 * t163 * t164;
t254 = (t148 * t250 - t162 * t252) / t161 ^ 2;
t245 = t200 * t203;
t218 = t182 * t245 - t199 * t213;
t150 = t218 * t246;
t145 = t219 * t150 + t174 * t213 + t231;
t253 = t145 * t214;
t201 = t199 * t200;
t251 = (t168 * t182 * t200 + t180 * t201 * t238) * t197 / t179 ^ 2;
t222 = qJD(2) * t256 + qJD(1);
t237 = qJD(2) * t205;
t166 = -qJD(1) * t225 - t206 * t237 + t222 * t241;
t249 = t166 * t156;
t248 = t174 * t214;
t247 = t175 * t214;
t242 = t202 * t206;
t181 = t214 ^ 2;
t153 = t181 * t156 + 0.1e1;
t236 = 0.2e1 * (-t181 * t255 + t214 * t249) / t153 ^ 2;
t235 = 0.2e1 * t255;
t234 = 0.2e1 * t254;
t233 = -0.2e1 * t251;
t232 = t172 * t252;
t169 = -qJD(1) * t226 - t204 * t238 + t222 * t240;
t223 = t198 * t242 - t169;
t220 = t193 * t163 - t194 * t250;
t217 = t198 * t213 + t204 * t239;
t171 = t193 * t242 + t194 * t213;
t170 = t193 * t213 - t194 * t242;
t159 = 0.1e1 / t161;
t151 = 0.1e1 / t153;
t146 = t257 * t214;
t142 = (t218 * t233 + (t168 * t245 + t169 * t199 + (-t213 * t245 + (0.2e1 * t201 * t203 ^ 2 + t199) * t182) * qJD(2)) * t177) * t196;
t140 = -0.2e1 * t254 + 0.2e1 * (t148 * t164 * t159 + (-t159 * t252 - t164 * t254) * t172) * t172;
t1 = [(-t214 * t199 * t233 + (-t166 * t199 - t214 * t229) * t177) * t196, t142, 0, 0, 0, 0; t182 * t155 * t236 + (-t168 * t155 + (t143 * t182 + t146 * t166) * t156) * t151 - ((t146 * t235 - t257 * t249) * t151 + (t146 * t236 + ((t147 * t177 * t230 + t233) * t248 + (0.2e1 * t230 * t251 - t147 + (-t215 * t196 + t147) * t177) * t247) * t151) * t156) * t214 (-t155 * t186 - t156 * t253) * t236 + (-t235 * t253 + t167 * t155 + (-t186 * t143 + t145 * t166 + (t202 * t237 - t142 * t182 - t150 * t168 + (t150 * t243 + t213) * t147) * t247 + (t147 * t150 * t182 - t169 + (t142 * t205 + (-qJD(2) * t150 - t147) * t203) * t202) * t248) * t156) * t151, 0, 0, 0, 0; (-t163 * t170 + t171 * t250) * t234 + ((t223 * t193 + t217 * t194) * t163 + 0.2e1 * t171 * t232 + (-t170 * t149 - (-t217 * t193 + t223 * t194) * t172 - t171 * t148) * t164) * t159, -t220 * t214 * t234 + (t220 * t166 - ((-t163 * t198 - 0.2e1 * t232) * t194 + (t148 * t194 + (-t172 * t198 + t149) * t193) * t164) * t214) * t159, t140, 0, t140, 0;];
JaD_rot  = t1;
