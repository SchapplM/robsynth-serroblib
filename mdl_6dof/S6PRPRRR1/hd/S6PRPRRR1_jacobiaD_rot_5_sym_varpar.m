% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:43
% EndTime: 2019-02-26 19:53:44
% DurationCPUTime: 0.51s
% Computational Cost: add. (2232->61), mult. (5863->138), div. (299->12), fcn. (7536->13), ass. (0->75)
t217 = cos(pkin(6));
t212 = sin(pkin(12));
t215 = cos(pkin(12));
t218 = sin(qJ(2));
t219 = cos(qJ(2));
t230 = t219 * t212 + t218 * t215;
t200 = t230 * t217;
t203 = t218 * t212 - t219 * t215;
t201 = t203 * qJD(2);
t213 = sin(pkin(11));
t216 = cos(pkin(11));
t227 = t203 * t217;
t185 = -t213 * t230 - t216 * t227;
t214 = sin(pkin(6));
t198 = t203 * t214;
t171 = atan2(t185, t198);
t166 = sin(t171);
t167 = cos(t171);
t160 = t166 * t185 + t167 * t198;
t157 = 0.1e1 / t160;
t211 = qJ(4) + qJ(5);
t208 = sin(t211);
t209 = cos(t211);
t231 = -t213 * t200 - t216 * t203;
t236 = t213 * t214;
t177 = t208 * t236 + t209 * t231;
t173 = 0.1e1 / t177;
t195 = 0.1e1 / t198;
t158 = 0.1e1 / t160 ^ 2;
t174 = 0.1e1 / t177 ^ 2;
t196 = 0.1e1 / t198 ^ 2;
t182 = t185 ^ 2;
t170 = t182 * t196 + 0.1e1;
t168 = 0.1e1 / t170;
t226 = qJD(2) * t200;
t178 = t213 * t201 - t216 * t226;
t199 = t230 * t214;
t192 = qJD(2) * t199;
t240 = t185 * t196;
t151 = (t178 * t195 - t192 * t240) * t168;
t232 = -t166 * t198 + t167 * t185;
t148 = t232 * t151 + t166 * t178 + t167 * t192;
t245 = t148 * t157 * t158;
t176 = t208 * t231 - t209 * t236;
t172 = t176 ^ 2;
t163 = t172 * t174 + 0.1e1;
t194 = t217 * t201;
t202 = t230 * qJD(2);
t181 = t213 * t194 - t216 * t202;
t210 = qJD(4) + qJD(5);
t233 = t210 * t236 + t181;
t238 = t231 * t210;
t164 = t233 * t208 + t209 * t238;
t241 = t174 * t176;
t165 = -t208 * t238 + t233 * t209;
t242 = t165 * t173 * t174;
t244 = (t164 * t241 - t172 * t242) / t163 ^ 2;
t187 = t213 * t227 - t216 * t230;
t243 = t158 * t187;
t239 = t185 * t199;
t237 = t192 * t195 * t196;
t229 = -t173 * t208 + t209 * t241;
t184 = -t216 * t200 + t213 * t203;
t228 = -t184 * t195 + t196 * t239;
t193 = t214 * t201;
t183 = t187 ^ 2;
t180 = t216 * t201 + t213 * t226;
t179 = t216 * t194 + t213 * t202;
t161 = 0.1e1 / t163;
t155 = t183 * t158 + 0.1e1;
t152 = t228 * t168;
t149 = -t232 * t152 + t166 * t184 + t167 * t199;
t147 = 0.2e1 * t228 / t170 ^ 2 * (t178 * t240 - t182 * t237) + (0.2e1 * t237 * t239 + t179 * t195 + (-t178 * t199 - t184 * t192 + t185 * t193) * t196) * t168;
t145 = -0.2e1 * t244 + 0.2e1 * (t161 * t164 * t174 + (-t161 * t242 - t174 * t244) * t176) * t176;
t1 = [0, t147, 0, 0, 0, 0; 0, 0.2e1 * (-t149 * t243 - t157 * t231) / t155 ^ 2 * (t180 * t243 - t183 * t245) + (t181 * t157 + (-t148 * t231 + t149 * t180) * t158 + (-0.2e1 * t149 * t245 + ((t147 * t185 - t152 * t178 - t193 + (t152 * t198 + t184) * t151) * t167 + (-t147 * t198 + t152 * t192 + t179 + (t152 * t185 - t199) * t151) * t166) * t158) * t187) / t155, 0, 0, 0, 0; 0, 0.2e1 * t229 * t187 * t244 + (-t229 * t180 + ((t173 * t210 + 0.2e1 * t176 * t242) * t209 + (-t164 * t209 + (t176 * t210 - t165) * t208) * t174) * t187) * t161, 0, t145, t145, 0;];
JaD_rot  = t1;
