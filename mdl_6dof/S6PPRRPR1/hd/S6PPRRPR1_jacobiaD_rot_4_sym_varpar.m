% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRPR1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiaD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:19
% EndTime: 2019-02-26 19:40:20
% DurationCPUTime: 0.68s
% Computational Cost: add. (2591->65), mult. (8150->143), div. (281->12), fcn. (10563->15), ass. (0->78)
t208 = sin(pkin(7));
t209 = sin(pkin(6));
t211 = cos(pkin(6));
t207 = sin(pkin(12));
t250 = sin(pkin(11));
t234 = t250 * t207;
t210 = cos(pkin(12));
t251 = cos(pkin(11));
t235 = t251 * t210;
t252 = cos(pkin(7));
t256 = (t211 * t235 - t234) * t252 - t208 * t209 * t251;
t233 = t250 * t210;
t236 = t251 * t207;
t224 = t211 * t233 + t236;
t238 = t209 * t250;
t255 = -t208 * t238 + t224 * t252;
t202 = t211 * t236 + t233;
t213 = sin(qJ(3));
t253 = cos(qJ(3));
t188 = t202 * t253 + t256 * t213;
t237 = t210 * t252;
t241 = t208 * t211;
t254 = (-t207 * t213 + t253 * t237) * t209 + t253 * t241;
t186 = t202 * t213 - t256 * t253;
t179 = atan2(-t186, -t254);
t174 = sin(t179);
t175 = cos(t179);
t162 = -t174 * t186 - t175 * t254;
t159 = 0.1e1 / t162;
t203 = -t211 * t234 + t235;
t190 = t203 * t253 - t255 * t213;
t199 = t224 * t208 + t252 * t238;
t212 = sin(qJ(4));
t214 = cos(qJ(4));
t173 = t190 * t214 + t199 * t212;
t169 = 0.1e1 / t173;
t194 = 0.1e1 / t254;
t160 = 0.1e1 / t162 ^ 2;
t170 = 0.1e1 / t173 ^ 2;
t195 = 0.1e1 / t254 ^ 2;
t184 = t186 ^ 2;
t178 = t184 * t195 + 0.1e1;
t176 = 0.1e1 / t178;
t181 = t188 * qJD(3);
t198 = t213 * t241 + (t253 * t207 + t213 * t237) * t209;
t192 = t198 * qJD(3);
t244 = t186 * t195;
t153 = (t181 * t194 + t192 * t244) * t176;
t229 = t174 * t254 - t175 * t186;
t150 = t229 * t153 - t174 * t181 + t175 * t192;
t249 = t150 * t159 * t160;
t172 = t190 * t212 - t199 * t214;
t168 = t172 ^ 2;
t165 = t168 * t170 + 0.1e1;
t189 = t203 * t213 + t255 * t253;
t182 = t189 * qJD(3);
t166 = t173 * qJD(4) - t182 * t212;
t245 = t170 * t172;
t240 = qJD(4) * t172;
t167 = -t182 * t214 - t240;
t246 = t167 * t169 * t170;
t248 = (t166 * t245 - t168 * t246) / t165 ^ 2;
t247 = t160 * t189;
t243 = t186 * t198;
t242 = t192 * t194 * t195;
t239 = -0.2e1 * t248;
t227 = -t169 * t212 + t214 * t245;
t226 = t188 * t194 + t195 * t243;
t191 = t254 * qJD(3);
t185 = t189 ^ 2;
t183 = t190 * qJD(3);
t180 = t186 * qJD(3);
t163 = 0.1e1 / t165;
t157 = t185 * t160 + 0.1e1;
t154 = t226 * t176;
t151 = t229 * t154 - t174 * t188 + t175 * t198;
t149 = -0.2e1 * t226 / t178 ^ 2 * (t181 * t244 + t184 * t242) + (0.2e1 * t242 * t243 - t180 * t194 + (t181 * t198 + t186 * t191 + t188 * t192) * t195) * t176;
t1 = [0, 0, t149, 0, 0, 0; 0, 0, 0.2e1 * (t151 * t247 - t159 * t190) / t157 ^ 2 * (t183 * t247 - t185 * t249) + (-t182 * t159 + (-t190 * t150 - t151 * t183) * t160 + (0.2e1 * t151 * t249 + (-(-t149 * t186 - t154 * t181 + t191 + (t154 * t254 - t188) * t153) * t175 - (t149 * t254 - t154 * t192 + t180 + (t154 * t186 - t198) * t153) * t174) * t160) * t189) / t157, 0, 0, 0; 0, 0, t227 * t189 * t239 + (t227 * t183 + ((-qJD(4) * t169 - 0.2e1 * t172 * t246) * t214 + (t166 * t214 + (t167 - t240) * t212) * t170) * t189) * t163, t239 + 0.2e1 * (t163 * t166 * t170 + (-t163 * t246 - t170 * t248) * t172) * t172, 0, 0;];
JaD_rot  = t1;
