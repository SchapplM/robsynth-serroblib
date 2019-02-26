% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR12_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:00:24
% EndTime: 2019-02-26 22:00:25
% DurationCPUTime: 0.60s
% Computational Cost: add. (1905->92), mult. (4665->200), div. (686->14), fcn. (5929->11), ass. (0->95)
t204 = sin(qJ(2));
t205 = sin(qJ(1));
t206 = cos(qJ(2));
t207 = cos(qJ(1));
t255 = cos(pkin(6));
t227 = t207 * t255;
t184 = t204 * t227 + t205 * t206;
t203 = sin(pkin(6));
t243 = t203 * t204;
t178 = atan2(-t184, t243);
t174 = sin(t178);
t175 = cos(t178);
t181 = t184 ^ 2;
t197 = 0.1e1 / t203 ^ 2;
t200 = 0.1e1 / t204 ^ 2;
t179 = t181 * t197 * t200 + 0.1e1;
t176 = 0.1e1 / t179;
t196 = 0.1e1 / t203;
t199 = 0.1e1 / t204;
t230 = t184 * t196 * t199;
t256 = t176 * (t175 * t230 + t174) - t174;
t158 = -t174 * t184 + t175 * t243;
t155 = 0.1e1 / t158;
t228 = t205 * t255;
t186 = t207 * t204 + t206 * t228;
t202 = qJ(4) + qJ(5);
t194 = sin(t202);
t195 = cos(t202);
t242 = t203 * t205;
t171 = t186 * t194 + t195 * t242;
t167 = 0.1e1 / t171;
t156 = 0.1e1 / t158 ^ 2;
t168 = 0.1e1 / t171 ^ 2;
t221 = qJD(2) * t255 + qJD(1);
t225 = t204 * t228;
t237 = qJD(2) * t204;
t239 = t207 * t206;
t165 = -qJD(1) * t225 - t205 * t237 + t221 * t239;
t236 = qJD(2) * t206;
t229 = t200 * t236;
t214 = -t165 * t199 + t184 * t229;
t245 = t176 * t196;
t147 = t214 * t245;
t218 = -t174 * t243 - t175 * t184;
t231 = t175 * t203 * t206;
t143 = qJD(2) * t231 + t147 * t218 - t174 * t165;
t254 = t143 * t155 * t156;
t224 = t206 * t227;
t240 = t205 * t204;
t183 = -t224 + t240;
t244 = t200 * t206;
t217 = t183 * t199 + t184 * t244;
t148 = t217 * t245;
t144 = t148 * t218 + t174 * t183 + t231;
t187 = -t225 + t239;
t253 = t144 * t187;
t198 = qJD(4) + qJD(5);
t238 = qJD(1) * t203;
t215 = t186 * t198 + t207 * t238;
t162 = -qJD(1) * t224 - t207 * t236 + t221 * t240;
t223 = t198 * t242 + t162;
t149 = t194 * t215 + t195 * t223;
t170 = -t186 * t195 + t194 * t242;
t166 = t170 ^ 2;
t161 = t166 * t168 + 0.1e1;
t248 = t168 * t170;
t150 = -t194 * t223 + t195 * t215;
t251 = t150 * t167 * t168;
t252 = (t149 * t248 - t166 * t251) / t161 ^ 2;
t201 = t199 * t200;
t250 = (t165 * t184 * t200 - t181 * t201 * t236) * t197 / t179 ^ 2;
t163 = qJD(1) * t184 + qJD(2) * t186;
t249 = t163 * t156;
t247 = t174 * t187;
t246 = t175 * t187;
t241 = t203 * t207;
t182 = t187 ^ 2;
t153 = t182 * t156 + 0.1e1;
t235 = 0.2e1 * (-t182 * t254 - t187 * t249) / t153 ^ 2;
t234 = 0.2e1 * t254;
t233 = 0.2e1 * t252;
t232 = -0.2e1 * t250;
t226 = 0.2e1 * t170 * t251;
t164 = qJD(1) * t186 + t184 * qJD(2);
t222 = t198 * t241 + t164;
t219 = t195 * t167 + t194 * t248;
t216 = -t183 * t198 - t205 * t238;
t173 = -t183 * t194 + t195 * t241;
t172 = t183 * t195 + t194 * t241;
t159 = 0.1e1 / t161;
t151 = 0.1e1 / t153;
t146 = t256 * t187;
t142 = (t217 * t232 + (t165 * t244 + t164 * t199 + (-t183 * t244 + (-0.2e1 * t201 * t206 ^ 2 - t199) * t184) * qJD(2)) * t176) * t196;
t140 = -0.2e1 * t252 + 0.2e1 * (t149 * t168 * t159 + (-t159 * t251 - t168 * t252) * t170) * t170;
t1 = [(0.2e1 * t187 * t199 * t250 + (t163 * t199 + t187 * t229) * t176) * t196, t142, 0, 0, 0, 0; t184 * t155 * t235 + (-t165 * t155 + (t143 * t184 + t146 * t163) * t156) * t151 + ((t146 * t234 + t256 * t249) * t151 + (t146 * t235 + (-(-t147 * t176 * t230 + t232) * t247 - (t230 * t232 - t147 + (-t196 * t214 + t147) * t176) * t246) * t151) * t156) * t187 (t155 * t186 + t156 * t253) * t235 + (t234 * t253 + t162 * t155 + (t186 * t143 + t144 * t163 - (-t203 * t237 - t142 * t184 - t148 * t165 + (-t148 * t243 + t183) * t147) * t246 - (t147 * t148 * t184 + t164 + (-t142 * t204 + (-qJD(2) * t148 - t147) * t206) * t203) * t247) * t156) * t151, 0, 0, 0, 0; (-t167 * t172 + t173 * t248) * t233 + ((t194 * t216 + t195 * t222) * t167 + t173 * t226 + (-t172 * t150 - (-t194 * t222 + t195 * t216) * t170 - t173 * t149) * t168) * t159, t219 * t187 * t233 + (t219 * t163 + ((t167 * t198 + t226) * t194 + (-t149 * t194 + (-t170 * t198 + t150) * t195) * t168) * t187) * t159, 0, t140, t140, 0;];
JaD_rot  = t1;
