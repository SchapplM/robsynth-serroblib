% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:25
% EndTime: 2019-02-26 21:49:26
% DurationCPUTime: 1.22s
% Computational Cost: add. (3233->93), mult. (9702->187), div. (711->12), fcn. (12089->11), ass. (0->93)
t279 = sin(qJ(4));
t280 = sin(qJ(2));
t282 = cos(qJ(4));
t283 = cos(qJ(2));
t201 = -t283 * t279 + t280 * t282;
t281 = sin(qJ(1));
t192 = t201 * t281;
t200 = t280 * t279 + t283 * t282;
t179 = atan2(t192, t200);
t174 = sin(t179);
t175 = cos(t179);
t194 = t200 * t281;
t244 = -t174 * t200 + t175 * t192;
t190 = t192 ^ 2;
t198 = 0.1e1 / t200 ^ 2;
t178 = t190 * t198 + 0.1e1;
t176 = 0.1e1 / t178;
t197 = 0.1e1 / t200;
t263 = t192 * t201;
t240 = -t194 * t197 - t198 * t263;
t291 = t240 * t176;
t151 = -t174 * t194 + t175 * t201 + t244 * t291;
t166 = t174 * t192 + t175 * t200;
t163 = 0.1e1 / t166;
t284 = cos(qJ(1));
t195 = t200 * t284;
t164 = 0.1e1 / t166 ^ 2;
t196 = t201 * t284;
t191 = t196 ^ 2;
t160 = t191 * t164 + 0.1e1;
t299 = qJD(2) - qJD(4);
t168 = -t192 * qJD(1) + t195 * t299;
t272 = t164 * t196;
t169 = -t196 * qJD(1) - t194 * t299;
t181 = t299 * t201;
t264 = t192 * t198;
t242 = -t169 * t197 + t181 * t264;
t154 = t242 * t176;
t149 = t244 * t154 - t174 * t169 - t175 * t181;
t278 = t149 * t163 * t164;
t261 = 0.2e1 * (t168 * t272 - t191 * t278) / t160 ^ 2;
t309 = (t151 * t272 + t163 * t195) * t261;
t167 = t194 * qJD(1) + t196 * t299;
t219 = sin(qJ(5));
t220 = cos(qJ(5));
t189 = t195 * t220 - t281 * t219;
t254 = qJD(1) * t284;
t161 = t189 * qJD(5) - t167 * t219 + t220 * t254;
t238 = -t195 * t219 - t281 * t220;
t300 = t238 * qJD(5);
t162 = -t167 * t220 - t219 * t254 + t300;
t182 = t238 ^ 2;
t184 = 0.1e1 / t189 ^ 2;
t173 = t182 * t184 + 0.1e1;
t171 = 0.1e1 / t173;
t183 = 0.1e1 / t189;
t266 = t184 * t238;
t241 = -t219 * t183 - t220 * t266;
t274 = t162 * t183 * t184;
t286 = -0.2e1 * t238;
t251 = t274 * t286;
t308 = (t241 * t168 - t196 * (((-t162 - t300) * t219 - t161 * t220) * t184 + (qJD(5) * t183 + t251) * t220)) * t171;
t307 = -t195 * t149 + t151 * t168;
t260 = 0.2e1 * t278;
t306 = -t151 * t196 * t260 - t167 * t163;
t180 = t299 * t200;
t304 = (t200 * t291 + t194) * t154 + t291 * t169 - t180;
t170 = t195 * qJD(1) - t192 * t299;
t303 = -(-t192 * t291 - t201) * t154 - t291 * t181 + t170;
t293 = t181 * t198;
t267 = t197 * t293;
t301 = ((t169 * t201 - t180 * t192 - t181 * t194) * t198 - t170 * t197 - 0.2e1 * t263 * t267) * t176;
t265 = t192 * t197;
t289 = (-t175 * t265 + t174) * t176 - t174;
t285 = -0.2e1 * t196;
t277 = (-t161 * t266 - t182 * t274) / t173 ^ 2;
t276 = (-t169 * t264 + t190 * t267) / t178 ^ 2;
t270 = t171 * t184;
t269 = t174 * t196;
t268 = t175 * t196;
t259 = -0.2e1 * t277;
t258 = -0.2e1 * t276;
t257 = t184 * t277;
t256 = t197 * t276;
t255 = t161 * t270;
t253 = qJD(1) * t281;
t187 = -t194 * t220 - t284 * t219;
t239 = t194 * t219 - t284 * t220;
t158 = 0.1e1 / t160;
t153 = t289 * t196;
t147 = t240 * t258 + t301;
t146 = 0.2e1 * t240 * t276 - t301;
t1 = [t256 * t285 + (t168 * t197 + t196 * t293) * t176, t146, 0, t147, 0, 0; -t192 * t163 * t261 + (-t169 * t163 + (-t149 * t192 - t153 * t168) * t164) * t158 - (-t153 * t260 * t158 + (-t153 * t261 + ((t154 * t176 * t265 + t258) * t269 + (0.2e1 * t192 * t256 - t154 + (t154 - t242) * t176) * t268 + t289 * t168) * t158) * t164) * t196, t309 + (((t146 * t192 + t304) * t268 + (-t146 * t200 + t303) * t269 - t307) * t164 - t306) * t158, 0, -t309 + (((t147 * t192 - t304) * t268 + (-t147 * t200 - t303) * t269 + t307) * t164 + t306) * t158, 0, 0; (t257 * t286 - t255) * t187 - (-t162 * t270 + t183 * t259) * t239 + ((t187 * qJD(5) - t170 * t219 - t220 * t253) * t183 + (t239 * qJD(5) - t170 * t220 + t219 * t253) * t266 + t187 * t251) * t171, t241 * t277 * t285 + t308, 0, -t241 * t196 * t259 - t308, t259 + (t255 - (-t171 * t274 - t257) * t238) * t286, 0;];
JaD_rot  = t1;
