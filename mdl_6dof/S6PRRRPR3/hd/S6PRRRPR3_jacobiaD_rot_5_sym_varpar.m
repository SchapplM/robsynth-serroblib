% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:43
% EndTime: 2019-02-26 20:11:44
% DurationCPUTime: 0.98s
% Computational Cost: add. (7756->94), mult. (10920->196), div. (767->12), fcn. (14140->11), ass. (0->96)
t204 = sin(pkin(11));
t206 = cos(pkin(11));
t208 = sin(qJ(2));
t207 = cos(pkin(6));
t209 = cos(qJ(2));
t231 = t207 * t209;
t192 = -t204 * t208 + t206 * t231;
t185 = t192 * qJD(2);
t232 = t207 * t208;
t193 = t204 * t209 + t206 * t232;
t203 = qJ(3) + qJ(4);
t200 = sin(t203);
t202 = qJD(3) + qJD(4);
t205 = sin(pkin(6));
t235 = t205 * t206;
t222 = t200 * t235;
t201 = cos(t203);
t237 = t201 * t202;
t157 = t185 * t200 + t193 * t237 - t202 * t222;
t175 = t193 * t200 + t201 * t235;
t172 = t175 ^ 2;
t234 = t205 * t208;
t224 = t200 * t234;
t183 = -t207 * t201 + t224;
t181 = 0.1e1 / t183 ^ 2;
t165 = t172 * t181 + 0.1e1;
t163 = 0.1e1 / t165;
t229 = qJD(2) * t209;
t217 = t202 * t207 + t205 * t229;
t223 = t201 * t234;
t170 = t200 * t217 + t202 * t223;
t180 = 0.1e1 / t183;
t242 = t175 * t181;
t145 = (-t157 * t180 + t170 * t242) * t163;
t166 = atan2(-t175, t183);
t161 = sin(t166);
t162 = cos(t166);
t219 = -t161 * t183 - t162 * t175;
t141 = t145 * t219 - t161 * t157 + t162 * t170;
t156 = -t161 * t175 + t162 * t183;
t153 = 0.1e1 / t156;
t154 = 0.1e1 / t156 ^ 2;
t251 = t141 * t153 * t154;
t225 = t204 * t232;
t195 = t206 * t209 - t225;
t236 = t204 * t205;
t178 = t195 * t200 - t201 * t236;
t250 = 0.2e1 * t178 * t251;
t233 = t205 * t209;
t216 = -t180 * t192 + t233 * t242;
t249 = t200 * t216;
t243 = t170 * t180 * t181;
t248 = -0.2e1 * (t157 * t242 - t172 * t243) / t165 ^ 2;
t194 = t204 * t231 + t206 * t208;
t189 = 0.1e1 / t194;
t190 = 0.1e1 / t194 ^ 2;
t247 = t154 * t178;
t187 = t194 * qJD(2);
t220 = t202 * t236 - t187;
t159 = t195 * t237 + t200 * t220;
t246 = t159 * t154;
t245 = t161 * t178;
t244 = t162 * t178;
t179 = t195 * t201 + t200 * t236;
t241 = t179 * t195;
t240 = t189 * t194;
t239 = t194 * t200;
t238 = t200 * t202;
t230 = qJD(2) * t208;
t173 = t178 ^ 2;
t151 = t173 * t154 + 0.1e1;
t228 = 0.2e1 * (-t173 * t251 + t178 * t246) / t151 ^ 2;
t160 = -t195 * t238 + t201 * t220;
t174 = t179 ^ 2;
t169 = t174 * t190 + 0.1e1;
t188 = -qJD(2) * t225 + t206 * t229;
t191 = t189 * t190;
t227 = 0.2e1 * (t179 * t190 * t160 - t174 * t191 * t188) / t169 ^ 2;
t221 = -0.2e1 * t175 * t243;
t177 = t193 * t201 - t222;
t184 = t207 * t200 + t223;
t218 = -t177 * t180 + t184 * t242;
t186 = t193 * qJD(2);
t171 = t201 * t217 - t202 * t224;
t167 = 0.1e1 / t169;
t158 = -t193 * t238 + (-t202 * t235 + t185) * t201;
t149 = 0.1e1 / t151;
t147 = t163 * t249;
t146 = t218 * t163;
t144 = t178 * t189 * t227 + (t178 * t188 * t190 - t159 * t189) * t167;
t143 = (-t161 * t192 + t162 * t233) * t200 + t219 * t147;
t142 = t146 * t219 - t161 * t177 + t162 * t184;
t140 = t218 * t248 + (t184 * t221 - t158 * t180 + (t157 * t184 + t170 * t177 + t171 * t175) * t181) * t163;
t138 = t248 * t249 + (t216 * t237 + (t221 * t233 + t180 * t186 + (t170 * t192 + (t157 * t209 - t175 * t230) * t205) * t181) * t200) * t163;
t137 = (t142 * t247 - t153 * t179) * t228 + (t142 * t250 + t160 * t153 + (-t179 * t141 - t142 * t159 - (-t140 * t175 - t146 * t157 + t171 + (-t146 * t183 - t177) * t145) * t244 - (-t140 * t183 - t146 * t170 - t158 + (t146 * t175 - t184) * t145) * t245) * t154) * t149;
t1 = [0, t138, t140, t140, 0, 0; 0 (t143 * t247 + t153 * t239) * t228 + ((-t188 * t200 - t194 * t237) * t153 + (-t246 + t250) * t143 + (t239 * t141 - (-t138 * t175 - t147 * t157 + (-t200 * t230 + t209 * t237) * t205 + (-t147 * t183 - t192 * t200) * t145) * t244 - (-t192 * t237 - t138 * t183 - t147 * t170 + t186 * t200 + (t147 * t175 - t200 * t233) * t145) * t245) * t154) * t149, t137, t137, 0, 0; 0 (t190 * t241 + t201 * t240) * t227 + (t238 * t240 + (-t160 * t195 + t179 * t187) * t190 + (0.2e1 * t191 * t241 + (t190 * t194 - t189) * t201) * t188) * t167, t144, t144, 0, 0;];
JaD_rot  = t1;
