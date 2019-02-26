% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:59
% EndTime: 2019-02-26 22:27:01
% DurationCPUTime: 1.12s
% Computational Cost: add. (11090->123), mult. (8378->264), div. (1558->15), fcn. (10537->9), ass. (0->116)
t191 = qJ(3) + qJ(4) + pkin(10);
t190 = cos(t191);
t274 = 0.2e1 * t190;
t189 = sin(t191);
t199 = cos(qJ(2));
t200 = cos(qJ(1));
t247 = t200 * t190;
t230 = t199 * t247;
t268 = sin(qJ(1));
t175 = t189 * t268 + t230;
t169 = 0.1e1 / t175 ^ 2;
t198 = sin(qJ(2));
t193 = t198 ^ 2;
t197 = t200 ^ 2;
t252 = t193 * t197;
t231 = t169 * t252;
t165 = 0.1e1 + t231;
t223 = qJD(1) * t268;
t245 = qJD(2) * t199;
t208 = t193 * t200 * t223 - t197 * t198 * t245;
t192 = qJD(3) + qJD(4);
t244 = qJD(2) * t200;
t225 = t198 * t244;
t211 = t199 * t223 + t225;
t229 = t268 * t192;
t248 = t200 * t189;
t154 = (-t192 * t199 + qJD(1)) * t248 + (t229 - t211) * t190;
t168 = 0.1e1 / t175;
t263 = t154 * t168 * t169;
t218 = t252 * t263;
t273 = (-t169 * t208 - t218) / t165 ^ 2;
t227 = t268 * t199;
t171 = t189 * t227 + t247;
t272 = t171 * t192;
t250 = t198 * t200;
t153 = t171 * qJD(1) - t192 * t230 + (t225 - t229) * t189;
t174 = -t190 * t268 + t199 * t248;
t186 = 0.1e1 / t189;
t187 = 0.1e1 / t189 ^ 2;
t194 = 0.1e1 / t198;
t195 = 0.1e1 / t198 ^ 2;
t226 = t195 * t245;
t254 = t190 * t192;
t256 = t186 * t194;
t271 = t174 * (t187 * t194 * t254 + t186 * t226) + t153 * t256;
t251 = t198 * t189;
t161 = atan2(-t171, t251);
t158 = cos(t161);
t157 = sin(t161);
t262 = t157 * t171;
t152 = t158 * t251 - t262;
t149 = 0.1e1 / t152;
t150 = 0.1e1 / t152 ^ 2;
t270 = -0.2e1 * t171;
t269 = 0.2e1 * t174;
t166 = t171 ^ 2;
t255 = t187 * t195;
t162 = t166 * t255 + 0.1e1;
t159 = 0.1e1 / t162;
t253 = t190 * t198;
t212 = t189 * t245 + t192 * t253;
t234 = t171 * t255;
t228 = t268 * t198;
t216 = qJD(2) * t228;
t246 = qJD(1) * t200;
t155 = -t189 * t216 - t192 * t248 - t190 * t223 + (t189 * t246 + t190 * t229) * t199;
t236 = t155 * t256;
t141 = (t212 * t234 - t236) * t159;
t209 = -t141 * t171 + t212;
t136 = (-t141 * t251 - t155) * t157 + t209 * t158;
t151 = t149 * t150;
t267 = t136 * t151;
t188 = t186 * t187;
t196 = t194 / t193;
t232 = t195 * t254;
t266 = (t155 * t234 + (-t187 * t196 * t245 - t188 * t232) * t166) / t162 ^ 2;
t265 = t150 * t174;
t264 = t153 * t150;
t261 = t157 * t174;
t260 = t157 * t198;
t259 = t158 * t171;
t258 = t158 * t174;
t257 = t158 * t199;
t249 = t199 * t200;
t167 = t174 ^ 2;
t147 = t167 * t150 + 0.1e1;
t243 = 0.2e1 * (-t167 * t267 - t174 * t264) / t147 ^ 2;
t242 = -0.2e1 * t266;
t241 = 0.2e1 * t273;
t240 = t151 * t269;
t239 = t194 * t266;
t238 = t150 * t261;
t235 = t171 * t256;
t233 = t186 * t195 * t199;
t214 = t171 * t233 + t268;
t148 = t214 * t159;
t224 = t268 - t148;
t222 = t149 * t243;
t221 = t150 * t243;
t220 = t250 * t269;
t219 = t186 * t239;
t173 = t190 * t227 - t248;
t215 = t171 * t187 * t190 - t173 * t186;
t213 = t169 * t173 * t200 - t168 * t268;
t163 = 0.1e1 / t165;
t156 = t175 * qJD(1) - t190 * t216 - t272;
t145 = 0.1e1 / t147;
t144 = t215 * t194 * t159;
t140 = (-t157 + (t158 * t235 + t157) * t159) * t174;
t139 = -t148 * t259 + (t224 * t260 + t257) * t189;
t138 = t158 * t253 - t157 * t173 + (-t157 * t251 - t259) * t144;
t137 = t169 * t220 * t273 + (t220 * t263 + (t153 * t250 + (t198 * t223 - t199 * t244) * t174) * t169) * t163;
t135 = t214 * t242 + (t155 * t233 + t246 + (-t187 * t199 * t232 + (-0.2e1 * t196 * t199 ^ 2 - t194) * t186 * qJD(2)) * t171) * t159;
t133 = -0.2e1 * t215 * t239 + (-t215 * t226 + ((-t156 - t272) * t186 + (t188 * t254 * t270 + (t173 * t192 + t155) * t187) * t190) * t194) * t159;
t132 = (t138 * t265 - t149 * t175) * t243 + (t138 * t264 + t154 * t149 + (t138 * t240 - t175 * t150) * t136 - (t190 * t245 - t192 * t251 - t133 * t171 - t144 * t155 + (-t144 * t251 - t173) * t141) * t150 * t258 - (-t156 + (-t133 * t189 - t141 * t190) * t198 - t209 * t144) * t238) * t145;
t1 = [t271 * t159 + t219 * t269, t135, t133, t133, 0, 0; t171 * t222 + (-t155 * t149 + (t136 * t171 + t140 * t153) * t150) * t145 + (t140 * t221 + (0.2e1 * t140 * t267 + (t153 * t159 - t153 - (-t141 * t159 * t235 + t242) * t174) * t150 * t157 + (-(t219 * t270 - t141) * t265 + (-(t141 + t236) * t174 + t271 * t171) * t150 * t159) * t158) * t145) * t174, t139 * t174 * t221 + (-(-t135 * t259 + (t141 * t262 - t155 * t158) * t148) * t265 + (-t149 * t250 - (-t148 * t260 + t157 * t228 + t257) * t265) * t254 + (t136 * t240 + t264) * t139) * t145 + (t222 * t250 + ((-t149 * t244 - (qJD(2) * t224 - t141) * t238) * t199 + (t149 * t223 + (t200 * t136 - (-t135 + t246) * t261 - (t141 * t224 - qJD(2)) * t258) * t150) * t198) * t145) * t189, t132, t132, 0, 0; t213 * t198 * t241 + (-t213 * t245 + ((qJD(1) * t168 + 0.2e1 * t173 * t263) * t200 + (-t154 * t268 - t156 * t200 + t173 * t223) * t169) * t198) * t163 (t168 * t249 + t190 * t231) * t241 + (t218 * t274 + t211 * t168 + (t189 * t192 * t252 + t154 * t249 + t208 * t274) * t169) * t163, t137, t137, 0, 0;];
JaD_rot  = t1;
