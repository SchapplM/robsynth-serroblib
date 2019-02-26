% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:40
% EndTime: 2019-02-26 21:18:41
% DurationCPUTime: 0.69s
% Computational Cost: add. (4138->96), mult. (4025->205), div. (771->12), fcn. (4686->9), ass. (0->99)
t196 = cos(qJ(1));
t255 = 0.2e1 * t196;
t192 = t196 ^ 2;
t194 = qJ(3) + qJ(4);
t186 = sin(t194);
t181 = 0.1e1 / t186 ^ 2;
t188 = cos(t194);
t184 = t188 ^ 2;
t238 = t181 * t184;
t178 = t192 * t238 + 0.1e1;
t176 = 0.1e1 / t178;
t180 = 0.1e1 / t186;
t195 = sin(qJ(1));
t228 = qJD(1) * t195;
t220 = t188 * t228;
t190 = qJD(3) + qJD(4);
t234 = t190 * t196;
t221 = t181 * t234;
t150 = ((t186 * t234 + t220) * t180 + t184 * t221) * t176;
t254 = -t150 + t234;
t229 = t196 * t188;
t175 = atan2(-t229, t186);
t173 = sin(t175);
t174 = cos(t175);
t160 = -t173 * t229 + t174 * t186;
t157 = 0.1e1 / t160;
t193 = qJ(5) + qJ(6);
t185 = sin(t193);
t231 = t196 * t185;
t187 = cos(t193);
t232 = t195 * t187;
t170 = t186 * t232 + t231;
t166 = 0.1e1 / t170;
t158 = 0.1e1 / t160 ^ 2;
t167 = 0.1e1 / t170 ^ 2;
t191 = t195 ^ 2;
t237 = t184 * t191;
t155 = t158 * t237 + 0.1e1;
t227 = qJD(1) * t196;
t212 = t184 * t195 * t227;
t236 = t186 * t190;
t241 = t174 * t188;
t250 = t150 * t196;
t145 = (t190 - t250) * t241 + (t254 * t186 + t220) * t173;
t252 = t145 * t157 * t158;
t253 = (-t237 * t252 + (-t188 * t191 * t236 + t212) * t158) / t155 ^ 2;
t189 = qJD(5) + qJD(6);
t215 = qJD(1) * t186 + t189;
t235 = t190 * t195;
t204 = t188 * t235 + t215 * t196;
t216 = t186 * t189 + qJD(1);
t209 = t187 * t216;
t151 = t204 * t185 + t195 * t209;
t230 = t196 * t187;
t233 = t195 * t185;
t169 = t186 * t233 - t230;
t165 = t169 ^ 2;
t164 = t165 * t167 + 0.1e1;
t244 = t167 * t169;
t210 = t185 * t216;
t152 = t204 * t187 - t195 * t210;
t249 = t152 * t166 * t167;
t251 = (t151 * t244 - t165 * t249) / t164 ^ 2;
t183 = t188 * t184;
t239 = t180 * t188;
t207 = t190 * (-t180 * t181 * t183 - t239);
t248 = (-t181 * t212 + t192 * t207) / t178 ^ 2;
t247 = t158 * t188;
t246 = t158 * t195;
t245 = t166 * t185;
t243 = t169 * t187;
t242 = t173 * t196;
t240 = t180 * t184;
t226 = -0.2e1 * t252;
t225 = 0.2e1 * t251;
t224 = t188 * t253;
t223 = t188 * t248;
t222 = t188 * t246;
t219 = 0.1e1 + t238;
t218 = t248 * t255;
t217 = 0.2e1 * t169 * t249;
t214 = t174 * t176 * t240;
t213 = (-t176 + 0.1e1) * t188 * t173;
t211 = t219 * t195;
t208 = t167 * t243 - t245;
t206 = t208 * t195;
t205 = t190 * t229 - t215 * t195;
t172 = t186 * t230 - t233;
t171 = t186 * t231 + t232;
t163 = t219 * t196 * t176;
t161 = 0.1e1 / t164;
t153 = 0.1e1 / t155;
t149 = (-t196 * t214 + t213) * t195;
t148 = t186 * t242 + t241 + (-t173 * t186 - t174 * t229) * t163;
t146 = -t219 * t218 + (-qJD(1) * t211 + t207 * t255) * t176;
t143 = -0.2e1 * t251 + 0.2e1 * (t151 * t167 * t161 + (-t161 * t249 - t167 * t251) * t169) * t169;
t142 = t188 * t206 * t225 + (t206 * t236 + (-t208 * t227 + ((t166 * t189 + t217) * t187 + (-t151 * t187 + (t169 * t189 - t152) * t185) * t167) * t195) * t188) * t161;
t141 = 0.2e1 * (-t148 * t247 - t157 * t186) * t195 * t253 + ((t157 * t227 + (-t148 * t190 - t145) * t246) * t186 + (t157 * t235 + (-t146 * t174 * t196 + t254 * t173 + (t150 * t242 - t173 * t190 + t174 * t228) * t163) * t222 + (t158 * t227 + t195 * t226) * t148 + ((-t146 - t228) * t173 + ((t163 * t196 - 0.1e1) * t190 + (-t163 + t196) * t150) * t174) * t186 * t246) * t188) * t153;
t1 = [-0.2e1 * t195 * t180 * t223 + (-t190 * t211 + t227 * t239) * t176, 0, t146, t146, 0, 0; (0.2e1 * t157 * t224 + (t157 * t236 + (qJD(1) * t149 + t145) * t247) * t153) * t196 + (-0.2e1 * t158 * t224 * t149 + (((0.2e1 * t223 - t236 + (t240 * t250 + t236) * t176) * t173 + (t218 * t240 + t150 * t188 + (t183 * t221 + (-t150 + 0.2e1 * t234) * t188) * t176) * t174) * t222 + (-t158 * t236 + t188 * t226) * t149 + (t157 + ((t191 - t192) * t214 + t196 * t213) * t158) * t188 * qJD(1)) * t153) * t195, 0, t141, t141, 0, 0; (-t166 * t171 + t172 * t244) * t225 + (t172 * t217 + t196 * t166 * t209 + t205 * t245 + (t196 * t169 * t210 - t172 * t151 - t171 * t152 - t205 * t243) * t167) * t161, 0, t142, t142, t143, t143;];
JaD_rot  = t1;
