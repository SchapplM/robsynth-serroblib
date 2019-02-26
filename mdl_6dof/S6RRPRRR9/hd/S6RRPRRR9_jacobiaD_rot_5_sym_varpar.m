% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:50
% EndTime: 2019-02-26 21:58:51
% DurationCPUTime: 0.61s
% Computational Cost: add. (2443->93), mult. (4665->202), div. (686->14), fcn. (5929->11), ass. (0->95)
t204 = cos(qJ(2));
t205 = cos(qJ(1));
t255 = cos(pkin(6));
t226 = t205 * t255;
t224 = t204 * t226;
t202 = sin(qJ(2));
t203 = sin(qJ(1));
t240 = t203 * t202;
t181 = -t224 + t240;
t201 = sin(pkin(6));
t242 = t201 * t204;
t175 = atan2(-t181, -t242);
t173 = sin(t175);
t174 = cos(t175);
t179 = t181 ^ 2;
t196 = 0.1e1 / t201 ^ 2;
t199 = 0.1e1 / t204 ^ 2;
t178 = t179 * t196 * t199 + 0.1e1;
t176 = 0.1e1 / t178;
t195 = 0.1e1 / t201;
t198 = 0.1e1 / t204;
t229 = t181 * t195 * t198;
t256 = (t174 * t229 - t173) * t176 + t173;
t157 = -t173 * t181 - t174 * t242;
t154 = 0.1e1 / t157;
t227 = t203 * t255;
t225 = t202 * t227;
t239 = t205 * t204;
t185 = -t225 + t239;
t194 = pkin(12) + qJ(4) + qJ(5);
t192 = sin(t194);
t193 = cos(t194);
t243 = t201 * t203;
t172 = t185 * t193 + t192 * t243;
t162 = 0.1e1 / t172;
t155 = 0.1e1 / t157 ^ 2;
t163 = 0.1e1 / t172 ^ 2;
t212 = -t202 * t226 - t203 * t204;
t213 = -t205 * t202 - t204 * t227;
t167 = -t213 * qJD(1) - t212 * qJD(2);
t237 = qJD(2) * t202;
t228 = t199 * t237;
t214 = t167 * t198 + t181 * t228;
t245 = t176 * t195;
t146 = t214 * t245;
t218 = t173 * t242 - t174 * t181;
t230 = t174 * t201 * t202;
t142 = qJD(2) * t230 + t218 * t146 - t173 * t167;
t254 = t142 * t154 * t155;
t197 = qJD(4) + qJD(5);
t238 = qJD(1) * t201;
t215 = -t185 * t197 + t205 * t238;
t166 = t212 * qJD(1) + t213 * qJD(2);
t223 = t197 * t243 + t166;
t147 = t223 * t192 - t215 * t193;
t171 = t185 * t192 - t193 * t243;
t161 = t171 ^ 2;
t160 = t161 * t163 + 0.1e1;
t249 = t163 * t171;
t148 = t215 * t192 + t223 * t193;
t251 = t148 * t162 * t163;
t253 = (t147 * t249 - t161 * t251) / t160 ^ 2;
t244 = t199 * t202;
t217 = t181 * t244 - t198 * t212;
t149 = t217 * t245;
t144 = t218 * t149 + t173 * t212 + t230;
t252 = t144 * t213;
t200 = t198 * t199;
t250 = (t167 * t181 * t199 + t179 * t200 * t237) * t196 / t178 ^ 2;
t221 = qJD(2) * t255 + qJD(1);
t236 = qJD(2) * t204;
t165 = -qJD(1) * t224 - t205 * t236 + t221 * t240;
t248 = t165 * t155;
t247 = t173 * t213;
t246 = t174 * t213;
t241 = t201 * t205;
t180 = t213 ^ 2;
t152 = t180 * t155 + 0.1e1;
t235 = 0.2e1 * (-t180 * t254 + t213 * t248) / t152 ^ 2;
t234 = 0.2e1 * t254;
t233 = 0.2e1 * t253;
t232 = -0.2e1 * t250;
t231 = t171 * t251;
t168 = -qJD(1) * t225 - t203 * t237 + t221 * t239;
t222 = t197 * t241 - t168;
t219 = t192 * t162 - t193 * t249;
t216 = t197 * t212 + t203 * t238;
t170 = t192 * t241 + t193 * t212;
t169 = t192 * t212 - t193 * t241;
t158 = 0.1e1 / t160;
t150 = 0.1e1 / t152;
t145 = t256 * t213;
t141 = (t217 * t232 + (t167 * t244 + t168 * t198 + (-t212 * t244 + (0.2e1 * t200 * t202 ^ 2 + t198) * t181) * qJD(2)) * t176) * t195;
t139 = -0.2e1 * t253 + 0.2e1 * (t147 * t163 * t158 + (-t158 * t251 - t163 * t253) * t171) * t171;
t1 = [(-t213 * t198 * t232 + (-t165 * t198 - t213 * t228) * t176) * t195, t141, 0, 0, 0, 0; t181 * t154 * t235 + (-t167 * t154 + (t142 * t181 + t145 * t165) * t155) * t150 - ((t145 * t234 - t256 * t248) * t150 + (t145 * t235 + ((t146 * t176 * t229 + t232) * t247 + (0.2e1 * t229 * t250 - t146 + (-t214 * t195 + t146) * t176) * t246) * t150) * t155) * t213 (-t154 * t185 - t155 * t252) * t235 + (-t234 * t252 + t166 * t154 + (-t185 * t142 + t144 * t165 + (t201 * t236 - t141 * t181 - t149 * t167 + (t149 * t242 + t212) * t146) * t246 + (t146 * t149 * t181 - t168 + (t141 * t204 + (-qJD(2) * t149 - t146) * t202) * t201) * t247) * t155) * t150, 0, 0, 0, 0; (-t162 * t169 + t170 * t249) * t233 + ((t222 * t192 + t216 * t193) * t162 + 0.2e1 * t170 * t231 + (-t169 * t148 - (-t216 * t192 + t222 * t193) * t171 - t170 * t147) * t163) * t158, -t219 * t213 * t233 + (t219 * t165 - ((-t162 * t197 - 0.2e1 * t231) * t193 + (t147 * t193 + (-t171 * t197 + t148) * t192) * t163) * t213) * t158, 0, t139, t139, 0;];
JaD_rot  = t1;
