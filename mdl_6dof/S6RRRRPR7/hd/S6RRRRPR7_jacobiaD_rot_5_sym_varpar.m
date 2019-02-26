% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:07
% EndTime: 2019-02-26 22:34:07
% DurationCPUTime: 0.68s
% Computational Cost: add. (2443->93), mult. (4665->202), div. (686->14), fcn. (5929->11), ass. (0->95)
t207 = cos(qJ(2));
t208 = cos(qJ(1));
t258 = cos(pkin(6));
t229 = t208 * t258;
t227 = t207 * t229;
t205 = sin(qJ(2));
t206 = sin(qJ(1));
t243 = t206 * t205;
t184 = -t227 + t243;
t204 = sin(pkin(6));
t245 = t204 * t207;
t178 = atan2(-t184, -t245);
t176 = sin(t178);
t177 = cos(t178);
t182 = t184 ^ 2;
t199 = 0.1e1 / t204 ^ 2;
t202 = 0.1e1 / t207 ^ 2;
t181 = t182 * t199 * t202 + 0.1e1;
t179 = 0.1e1 / t181;
t198 = 0.1e1 / t204;
t201 = 0.1e1 / t207;
t232 = t184 * t198 * t201;
t259 = (t177 * t232 - t176) * t179 + t176;
t160 = -t176 * t184 - t177 * t245;
t157 = 0.1e1 / t160;
t230 = t206 * t258;
t228 = t205 * t230;
t242 = t208 * t207;
t188 = -t228 + t242;
t197 = qJ(3) + qJ(4) + pkin(12);
t195 = sin(t197);
t196 = cos(t197);
t246 = t204 * t206;
t175 = t188 * t196 + t195 * t246;
t165 = 0.1e1 / t175;
t158 = 0.1e1 / t160 ^ 2;
t166 = 0.1e1 / t175 ^ 2;
t215 = -t205 * t229 - t206 * t207;
t216 = -t208 * t205 - t207 * t230;
t170 = -t216 * qJD(1) - t215 * qJD(2);
t240 = qJD(2) * t205;
t231 = t202 * t240;
t217 = t170 * t201 + t184 * t231;
t248 = t179 * t198;
t149 = t217 * t248;
t221 = t176 * t245 - t177 * t184;
t233 = t177 * t204 * t205;
t145 = qJD(2) * t233 + t221 * t149 - t176 * t170;
t257 = t145 * t157 * t158;
t200 = qJD(3) + qJD(4);
t241 = qJD(1) * t204;
t218 = -t188 * t200 + t208 * t241;
t169 = t215 * qJD(1) + t216 * qJD(2);
t226 = t200 * t246 + t169;
t150 = t226 * t195 - t218 * t196;
t174 = t188 * t195 - t196 * t246;
t164 = t174 ^ 2;
t163 = t164 * t166 + 0.1e1;
t252 = t166 * t174;
t151 = t218 * t195 + t226 * t196;
t254 = t151 * t165 * t166;
t256 = (t150 * t252 - t164 * t254) / t163 ^ 2;
t247 = t202 * t205;
t220 = t184 * t247 - t201 * t215;
t152 = t220 * t248;
t147 = t221 * t152 + t176 * t215 + t233;
t255 = t147 * t216;
t203 = t201 * t202;
t253 = (t170 * t184 * t202 + t182 * t203 * t240) * t199 / t181 ^ 2;
t224 = qJD(2) * t258 + qJD(1);
t239 = qJD(2) * t207;
t168 = -qJD(1) * t227 - t208 * t239 + t224 * t243;
t251 = t168 * t158;
t250 = t176 * t216;
t249 = t177 * t216;
t244 = t204 * t208;
t183 = t216 ^ 2;
t155 = t183 * t158 + 0.1e1;
t238 = 0.2e1 * (-t183 * t257 + t216 * t251) / t155 ^ 2;
t237 = 0.2e1 * t257;
t236 = 0.2e1 * t256;
t235 = -0.2e1 * t253;
t234 = t174 * t254;
t171 = -qJD(1) * t228 - t206 * t240 + t224 * t242;
t225 = t200 * t244 - t171;
t222 = t195 * t165 - t196 * t252;
t219 = t200 * t215 + t206 * t241;
t173 = t195 * t244 + t196 * t215;
t172 = t195 * t215 - t196 * t244;
t161 = 0.1e1 / t163;
t153 = 0.1e1 / t155;
t148 = t259 * t216;
t144 = (t220 * t235 + (t170 * t247 + t171 * t201 + (-t215 * t247 + (0.2e1 * t203 * t205 ^ 2 + t201) * t184) * qJD(2)) * t179) * t198;
t142 = -0.2e1 * t256 + 0.2e1 * (t150 * t166 * t161 + (-t161 * t254 - t166 * t256) * t174) * t174;
t1 = [(-t216 * t201 * t235 + (-t168 * t201 - t216 * t231) * t179) * t198, t144, 0, 0, 0, 0; t184 * t157 * t238 + (-t170 * t157 + (t145 * t184 + t148 * t168) * t158) * t153 - ((t148 * t237 - t259 * t251) * t153 + (t148 * t238 + ((t149 * t179 * t232 + t235) * t250 + (0.2e1 * t232 * t253 - t149 + (-t217 * t198 + t149) * t179) * t249) * t153) * t158) * t216 (-t157 * t188 - t158 * t255) * t238 + (-t237 * t255 + t169 * t157 + (-t188 * t145 + t147 * t168 + (t204 * t239 - t144 * t184 - t152 * t170 + (t152 * t245 + t215) * t149) * t249 + (t149 * t152 * t184 - t171 + (t144 * t207 + (-qJD(2) * t152 - t149) * t205) * t204) * t250) * t158) * t153, 0, 0, 0, 0; (-t165 * t172 + t173 * t252) * t236 + ((t225 * t195 + t219 * t196) * t165 + 0.2e1 * t173 * t234 + (-t172 * t151 - (-t219 * t195 + t225 * t196) * t174 - t173 * t150) * t166) * t161, -t222 * t216 * t236 + (t222 * t168 - ((-t165 * t200 - 0.2e1 * t234) * t196 + (t150 * t196 + (-t174 * t200 + t151) * t195) * t166) * t216) * t161, t142, t142, 0, 0;];
JaD_rot  = t1;
