% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:36:28
% EndTime: 2019-02-26 20:36:29
% DurationCPUTime: 0.77s
% Computational Cost: add. (2448->97), mult. (4991->216), div. (498->12), fcn. (6402->11), ass. (0->95)
t249 = sin(pkin(10));
t250 = cos(pkin(10));
t251 = sin(qJ(1));
t252 = cos(qJ(1));
t178 = t252 * t249 - t251 * t250;
t175 = t178 ^ 2;
t194 = sin(qJ(4));
t189 = t194 ^ 2;
t195 = cos(qJ(4));
t191 = 0.1e1 / t195 ^ 2;
t227 = t189 * t191;
t171 = t175 * t227 + 0.1e1;
t177 = -t249 * t251 - t250 * t252;
t173 = t177 * qJD(1);
t188 = t194 * t189;
t190 = 0.1e1 / t195;
t226 = t190 * t194;
t206 = qJD(4) * (t188 * t190 * t191 + t226);
t242 = (t173 * t178 * t227 + t175 * t206) / t171 ^ 2;
t254 = -0.2e1 * t242;
t215 = 0.1e1 + t227;
t253 = t178 * t215;
t232 = t178 * t194;
t168 = atan2(t232, t195);
t166 = sin(t168);
t167 = cos(t168);
t157 = t166 * t232 + t167 * t195;
t154 = 0.1e1 / t157;
t193 = qJ(5) + qJ(6);
t185 = sin(t193);
t186 = cos(t193);
t230 = t186 * t195;
t165 = -t177 * t230 + t178 * t185;
t159 = 0.1e1 / t165;
t155 = 0.1e1 / t157 ^ 2;
t160 = 0.1e1 / t165 ^ 2;
t176 = t177 ^ 2;
t233 = t176 * t189;
t148 = t155 * t233 + 0.1e1;
t174 = t178 * qJD(1);
t222 = qJD(4) * t195;
t169 = 0.1e1 / t171;
t217 = t178 * t222;
t224 = qJD(4) * t178;
t218 = t191 * t224;
t234 = t173 * t194;
t145 = ((t217 + t234) * t190 + t189 * t218) * t169;
t235 = t167 * t194;
t140 = (t145 * t178 - qJD(4)) * t235 + (t234 + (-t145 + t224) * t195) * t166;
t246 = t140 * t154 * t155;
t248 = (-t233 * t246 + (-t174 * t177 * t189 + t176 * t194 * t222) * t155) / t148 ^ 2;
t187 = qJD(5) + qJD(6);
t223 = qJD(4) * t194;
t203 = t174 * t195 + t177 * t223 + t178 * t187;
t229 = t187 * t195;
t209 = t177 * t229 + t173;
t143 = t185 * t203 - t186 * t209;
t231 = t185 * t195;
t164 = -t177 * t231 - t178 * t186;
t158 = t164 ^ 2;
t151 = t158 * t160 + 0.1e1;
t239 = t160 * t164;
t144 = t185 * t209 + t186 * t203;
t243 = t144 * t159 * t160;
t247 = (t143 * t239 - t158 * t243) / t151 ^ 2;
t153 = t169 * t253;
t236 = t166 * t195;
t141 = t178 * t236 - t235 + (t167 * t232 - t236) * t153;
t245 = t141 * t155;
t228 = t189 * t190;
t219 = t178 * t228;
t212 = t167 * t219;
t237 = t166 * t194;
t142 = (t237 + (t212 - t237) * t169) * t177;
t244 = t142 * t155;
t241 = t155 * t177;
t240 = t159 * t185;
t238 = t164 * t186;
t225 = t153 - t178;
t221 = 0.2e1 * t247;
t220 = -0.2e1 * t246;
t216 = t153 * t178 - 0.1e1;
t214 = -0.2e1 * t194 * t248;
t213 = 0.2e1 * t164 * t243;
t208 = t178 * t229 + t174;
t207 = t160 * t238 - t240;
t205 = t207 * t194;
t204 = t173 * t195 + t177 * t187 - t178 * t223;
t163 = t177 * t185 + t178 * t230;
t162 = -t177 * t186 + t178 * t231;
t149 = 0.1e1 / t151;
t146 = 0.1e1 / t148;
t138 = t253 * t254 + (t173 * t215 + 0.2e1 * t178 * t206) * t169;
t136 = -0.2e1 * t247 + 0.2e1 * (t143 * t160 * t149 + (-t149 * t243 - t160 * t247) * t164) * t164;
t1 = [t177 * t226 * t254 + (qJD(4) * t177 * t215 - t174 * t226) * t169, 0, 0, t138, 0, 0; t178 * t154 * t214 + (t154 * t217 + (t154 * t173 + (-t140 * t178 - t142 * t174) * t155) * t194) * t146 + (t214 * t244 + (t222 * t244 + (t142 * t220 + ((t166 * t222 + t212 * t254) * t177 + (-t166 * t174 + (t145 * t167 + 0.2e1 * t166 * t242) * t177) * t194 + (((-t145 * t219 - t222) * t177 + t174 * t194) * t166 + (-t174 * t219 + (t188 * t218 + t173 * t228 + (-t145 + 0.2e1 * t224) * t194) * t177) * t167) * t169) * t155) * t194) * t146) * t177, 0, 0, 0.2e1 * (t154 * t195 - t194 * t245) * t177 * t248 + ((t174 * t154 + (qJD(4) * t141 + t140) * t241) * t195 + (-t174 * t245 + (qJD(4) * t154 + t141 * t220 + ((t138 * t178 + t153 * t173) * t235 + (qJD(4) * t225 - t145 * t216) * t237) * t155) * t177 + ((-t138 + t173) * t166 + (qJD(4) * t216 - t145 * t225) * t167) * t195 * t241) * t194) * t146, 0, 0; (-t159 * t162 + t163 * t239) * t221 + (t163 * t213 + t208 * t159 * t186 + t204 * t240 + (t164 * t185 * t208 - t163 * t143 - t162 * t144 - t204 * t238) * t160) * t149, 0, 0, t177 * t205 * t221 + (t174 * t205 + (-t207 * t222 + ((t159 * t187 + t213) * t186 + (-t143 * t186 + (t164 * t187 - t144) * t185) * t160) * t194) * t177) * t149, t136, t136;];
JaD_rot  = t1;
