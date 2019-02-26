% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:41
% EndTime: 2019-02-26 21:54:41
% DurationCPUTime: 0.74s
% Computational Cost: add. (6414->98), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->98)
t199 = sin(qJ(1));
t260 = 0.2e1 * t199;
t196 = t199 ^ 2;
t191 = qJ(2) + pkin(11) + qJ(4);
t189 = sin(t191);
t185 = t189 ^ 2;
t190 = cos(t191);
t187 = 0.1e1 / t190 ^ 2;
t245 = t185 * t187;
t181 = t196 * t245 + 0.1e1;
t178 = 0.1e1 / t181;
t186 = 0.1e1 / t190;
t200 = cos(qJ(1));
t231 = qJD(1) * t200;
t221 = t189 * t231;
t195 = qJD(2) + qJD(4);
t237 = t195 * t199;
t224 = t187 * t237;
t152 = (-(-t190 * t237 - t221) * t186 + t185 * t224) * t178;
t259 = t152 - t237;
t198 = qJ(5) + qJ(6);
t192 = sin(t198);
t234 = t199 * t192;
t193 = cos(t198);
t238 = t193 * t200;
t174 = t190 * t238 + t234;
t235 = t199 * t189;
t177 = atan2(-t235, -t190);
t176 = cos(t177);
t175 = sin(t177);
t225 = t175 * t235;
t162 = -t176 * t190 - t225;
t159 = 0.1e1 / t162;
t168 = 0.1e1 / t174;
t160 = 0.1e1 / t162 ^ 2;
t169 = 0.1e1 / t174 ^ 2;
t258 = t178 - 0.1e1;
t247 = t176 * t189;
t147 = (-t152 * t199 + t195) * t247 + (t259 * t190 - t221) * t175;
t257 = t147 * t159 * t160;
t194 = qJD(5) + qJD(6);
t215 = -qJD(1) * t190 + t194;
t216 = t190 * t194 - qJD(1);
t236 = t195 * t200;
t222 = t189 * t236;
t239 = t192 * t200;
t154 = -t216 * t239 + (t215 * t199 - t222) * t193;
t256 = t154 * t168 * t169;
t184 = t189 * t185;
t242 = t186 * t189;
t209 = t195 * (t184 * t186 * t187 + t242);
t243 = t185 * t199;
t213 = t231 * t243;
t255 = (t187 * t213 + t196 * t209) / t181 ^ 2;
t254 = t160 * t189;
t253 = t160 * t200;
t210 = t190 * t234 + t238;
t153 = t210 * qJD(1) - t174 * t194 + t192 * t222;
t233 = t199 * t193;
t173 = t190 * t239 - t233;
t167 = t173 ^ 2;
t166 = t167 * t169 + 0.1e1;
t250 = t169 * t173;
t252 = 0.1e1 / t166 ^ 2 * (-t153 * t250 - t167 * t256);
t251 = t168 * t192;
t249 = t173 * t193;
t248 = t175 * t199;
t246 = t185 * t186;
t197 = t200 ^ 2;
t244 = t185 * t197;
t241 = t189 * t200;
t240 = t190 * t195;
t232 = qJD(1) * t199;
t157 = t160 * t244 + 0.1e1;
t230 = 0.2e1 * (-t244 * t257 + (t189 * t197 * t240 - t213) * t160) / t157 ^ 2;
t229 = 0.2e1 * t257;
t228 = -0.2e1 * t252;
t227 = t173 * t256;
t226 = t160 * t241;
t220 = 0.1e1 + t245;
t219 = t189 * t230;
t218 = -0.2e1 * t189 * t255;
t217 = t255 * t260;
t214 = t176 * t178 * t246;
t212 = t220 * t200;
t211 = t169 * t249 - t251;
t208 = t195 * t235 + t215 * t200;
t172 = -t190 * t233 + t239;
t164 = 0.1e1 / t166;
t163 = t220 * t199 * t178;
t155 = 0.1e1 / t157;
t151 = (t258 * t189 * t175 - t199 * t214) * t200;
t150 = -t190 * t248 + t247 + (t175 * t190 - t176 * t235) * t163;
t148 = -t220 * t217 + (qJD(1) * t212 + t209 * t260) * t178;
t145 = t228 + 0.2e1 * (-t153 * t169 * t164 + (-t164 * t256 - t169 * t252) * t173) * t173;
t144 = t211 * t228 * t241 + (t211 * t190 * t236 + (-t211 * t232 + ((-t168 * t194 - 0.2e1 * t227) * t193 + (-t153 * t193 + (-t173 * t194 + t154) * t192) * t169) * t200) * t189) * t164;
t143 = (t150 * t254 - t159 * t190) * t200 * t230 + ((-t159 * t232 + (-t150 * t195 - t147) * t253) * t190 + (-t159 * t236 - (-t148 * t176 * t199 - t259 * t175 + (t152 * t248 - t175 * t195 - t176 * t231) * t163) * t226 + (t160 * t232 + t200 * t229) * t150 - ((t148 - t231) * t175 + ((-t163 * t199 + 0.1e1) * t195 + (t163 - t199) * t152) * t176) * t190 * t253) * t189) * t155;
t1 = [t200 * t186 * t218 + (t195 * t212 - t232 * t242) * t178, t148, 0, t148, 0, 0; (t159 * t219 + (-t159 * t240 + (qJD(1) * t151 + t147) * t254) * t155) * t199 + (t160 * t219 * t151 + (-((t218 - t240 + (t152 * t186 * t243 + t240) * t178) * t175 + (t217 * t246 - t152 * t189 + (-t184 * t224 + (t152 - 0.2e1 * t237) * t189) * t178) * t176) * t226 + (-t160 * t240 + t189 * t229) * t151 + (-t159 + ((-t196 + t197) * t214 + t258 * t225) * t160) * t189 * qJD(1)) * t155) * t200, t143, 0, t143, 0, 0; 0.2e1 * (t168 * t210 + t172 * t250) * t252 + (0.2e1 * t172 * t227 - t216 * t168 * t233 + t208 * t251 + (-t216 * t173 * t234 + t172 * t153 + t154 * t210 - t208 * t249) * t169) * t164, t144, 0, t144, t145, t145;];
JaD_rot  = t1;
