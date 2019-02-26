% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:33
% EndTime: 2019-02-26 22:16:34
% DurationCPUTime: 0.76s
% Computational Cost: add. (6414->98), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->98)
t202 = sin(qJ(1));
t263 = 0.2e1 * t202;
t199 = t202 ^ 2;
t194 = qJ(2) + qJ(3) + pkin(11);
t192 = sin(t194);
t188 = t192 ^ 2;
t193 = cos(t194);
t190 = 0.1e1 / t193 ^ 2;
t248 = t188 * t190;
t184 = t199 * t248 + 0.1e1;
t181 = 0.1e1 / t184;
t189 = 0.1e1 / t193;
t203 = cos(qJ(1));
t234 = qJD(1) * t203;
t224 = t192 * t234;
t198 = qJD(2) + qJD(3);
t240 = t198 * t202;
t227 = t190 * t240;
t155 = (-(-t193 * t240 - t224) * t189 + t188 * t227) * t181;
t262 = t155 - t240;
t201 = qJ(5) + qJ(6);
t195 = sin(t201);
t237 = t202 * t195;
t196 = cos(t201);
t241 = t196 * t203;
t177 = t193 * t241 + t237;
t238 = t202 * t192;
t180 = atan2(-t238, -t193);
t179 = cos(t180);
t178 = sin(t180);
t228 = t178 * t238;
t165 = -t179 * t193 - t228;
t162 = 0.1e1 / t165;
t171 = 0.1e1 / t177;
t163 = 0.1e1 / t165 ^ 2;
t172 = 0.1e1 / t177 ^ 2;
t261 = t181 - 0.1e1;
t250 = t179 * t192;
t150 = (-t155 * t202 + t198) * t250 + (t262 * t193 - t224) * t178;
t260 = t150 * t162 * t163;
t197 = qJD(5) + qJD(6);
t218 = -qJD(1) * t193 + t197;
t219 = t193 * t197 - qJD(1);
t239 = t198 * t203;
t225 = t192 * t239;
t242 = t195 * t203;
t157 = -t219 * t242 + (t218 * t202 - t225) * t196;
t259 = t157 * t171 * t172;
t187 = t192 * t188;
t245 = t189 * t192;
t212 = t198 * (t187 * t189 * t190 + t245);
t246 = t188 * t202;
t216 = t234 * t246;
t258 = (t190 * t216 + t199 * t212) / t184 ^ 2;
t257 = t163 * t192;
t256 = t163 * t203;
t213 = t193 * t237 + t241;
t156 = t213 * qJD(1) - t177 * t197 + t195 * t225;
t236 = t202 * t196;
t176 = t193 * t242 - t236;
t170 = t176 ^ 2;
t169 = t170 * t172 + 0.1e1;
t253 = t172 * t176;
t255 = 0.1e1 / t169 ^ 2 * (-t156 * t253 - t170 * t259);
t254 = t171 * t195;
t252 = t176 * t196;
t251 = t178 * t202;
t249 = t188 * t189;
t200 = t203 ^ 2;
t247 = t188 * t200;
t244 = t192 * t203;
t243 = t193 * t198;
t235 = qJD(1) * t202;
t160 = t163 * t247 + 0.1e1;
t233 = 0.2e1 * (-t247 * t260 + (t192 * t200 * t243 - t216) * t163) / t160 ^ 2;
t232 = 0.2e1 * t260;
t231 = -0.2e1 * t255;
t230 = t176 * t259;
t229 = t163 * t244;
t223 = 0.1e1 + t248;
t222 = t192 * t233;
t221 = -0.2e1 * t192 * t258;
t220 = t258 * t263;
t217 = t179 * t181 * t249;
t215 = t223 * t203;
t214 = t172 * t252 - t254;
t211 = t198 * t238 + t218 * t203;
t175 = -t193 * t236 + t242;
t167 = 0.1e1 / t169;
t166 = t223 * t202 * t181;
t158 = 0.1e1 / t160;
t154 = (t261 * t192 * t178 - t202 * t217) * t203;
t153 = -t193 * t251 + t250 + (t178 * t193 - t179 * t238) * t166;
t151 = -t223 * t220 + (qJD(1) * t215 + t212 * t263) * t181;
t148 = t231 + 0.2e1 * (-t156 * t172 * t167 + (-t167 * t259 - t172 * t255) * t176) * t176;
t147 = t214 * t231 * t244 + (t214 * t193 * t239 + (-t214 * t235 + ((-t171 * t197 - 0.2e1 * t230) * t196 + (-t156 * t196 + (-t176 * t197 + t157) * t195) * t172) * t203) * t192) * t167;
t146 = (t153 * t257 - t162 * t193) * t203 * t233 + ((-t162 * t235 + (-t153 * t198 - t150) * t256) * t193 + (-t162 * t239 - (-t151 * t179 * t202 - t262 * t178 + (t155 * t251 - t178 * t198 - t179 * t234) * t166) * t229 + (t163 * t235 + t203 * t232) * t153 - ((t151 - t234) * t178 + ((-t166 * t202 + 0.1e1) * t198 + (t166 - t202) * t155) * t179) * t193 * t256) * t192) * t158;
t1 = [t203 * t189 * t221 + (t198 * t215 - t235 * t245) * t181, t151, t151, 0, 0, 0; (t162 * t222 + (-t162 * t243 + (qJD(1) * t154 + t150) * t257) * t158) * t202 + (t163 * t222 * t154 + (-((t221 - t243 + (t155 * t189 * t246 + t243) * t181) * t178 + (t220 * t249 - t155 * t192 + (-t187 * t227 + (t155 - 0.2e1 * t240) * t192) * t181) * t179) * t229 + (-t163 * t243 + t192 * t232) * t154 + (-t162 + ((-t199 + t200) * t217 + t261 * t228) * t163) * t192 * qJD(1)) * t158) * t203, t146, t146, 0, 0, 0; 0.2e1 * (t171 * t213 + t175 * t253) * t255 + (0.2e1 * t175 * t230 - t219 * t171 * t236 + t211 * t254 + (-t219 * t176 * t237 + t175 * t156 + t157 * t213 - t211 * t252) * t172) * t167, t147, t147, 0, t148, t148;];
JaD_rot  = t1;
