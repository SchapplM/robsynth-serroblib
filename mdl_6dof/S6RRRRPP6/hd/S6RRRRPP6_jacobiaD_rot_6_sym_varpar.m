% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:14
% EndTime: 2019-02-26 22:28:16
% DurationCPUTime: 1.17s
% Computational Cost: add. (6849->119), mult. (8141->258), div. (1608->14), fcn. (10341->9), ass. (0->116)
t183 = qJ(3) + qJ(4);
t173 = sin(t183);
t174 = cos(t183);
t187 = cos(qJ(1));
t235 = t187 * t174;
t185 = sin(qJ(1));
t186 = cos(qJ(2));
t237 = t185 * t186;
t159 = t173 * t237 + t235;
t175 = qJD(3) + qJD(4);
t264 = t159 * t175;
t184 = sin(qJ(2));
t263 = t184 * t187;
t236 = t187 * t173;
t238 = t185 * t174;
t262 = t186 * t236 - t238;
t176 = 0.1e1 / t184;
t177 = 0.1e1 / t184 ^ 2;
t178 = t176 * t177;
t261 = qJD(2) * (0.2e1 * t178 * t186 ^ 2 + t176);
t229 = qJD(2) * t187;
t209 = t184 * t229;
t232 = qJD(1) * t186;
t213 = t185 * t232;
t231 = qJD(1) * t187;
t144 = -t173 * t231 + (t209 + t213) * t174 + t262 * t175;
t163 = t185 * t173 + t186 * t235;
t170 = 0.1e1 / t174;
t171 = 0.1e1 / t174 ^ 2;
t230 = qJD(2) * t186;
t212 = t177 * t230;
t245 = t173 * t175;
t247 = t170 * t176;
t260 = (-t171 * t176 * t245 + t170 * t212) * t163 + t144 * t247;
t160 = t174 * t237 - t236;
t240 = t184 * t174;
t151 = atan2(-t160, t240);
t148 = cos(t151);
t147 = sin(t151);
t253 = t147 * t160;
t142 = t148 * t240 - t253;
t139 = 0.1e1 / t142;
t180 = 0.1e1 / t187;
t140 = 0.1e1 / t142 ^ 2;
t181 = 0.1e1 / t187 ^ 2;
t259 = -0.2e1 * t160;
t258 = 0.2e1 * t163;
t156 = t160 ^ 2;
t246 = t171 * t177;
t152 = t156 * t246 + 0.1e1;
t149 = 0.1e1 / t152;
t244 = t173 * t184;
t198 = t174 * t230 - t175 * t244;
t219 = t160 * t246;
t239 = t184 * t185;
t210 = qJD(2) * t239;
t146 = t163 * qJD(1) - t174 * t210 - t264;
t221 = t146 * t247;
t131 = (t198 * t219 - t221) * t149;
t196 = -t131 * t160 + t198;
t126 = (-t131 * t240 - t146) * t147 + t196 * t148;
t141 = t139 * t140;
t257 = t126 * t141;
t172 = t170 * t171;
t211 = t178 * t230;
t217 = t177 * t245;
t256 = (t146 * t219 + (-t171 * t211 + t172 * t217) * t156) / t152 ^ 2;
t255 = t140 * t163;
t254 = t144 * t140;
t252 = t147 * t163;
t251 = t147 * t184;
t250 = t148 * t160;
t249 = t148 * t163;
t248 = t148 * t186;
t243 = t177 * t181;
t242 = t177 * t186;
t241 = t181 * t185;
t218 = t170 * t242;
t201 = t160 * t218 + t185;
t138 = t201 * t149;
t234 = -t138 + t185;
t233 = qJD(1) * t185;
t158 = t163 ^ 2;
t137 = t158 * t140 + 0.1e1;
t228 = 0.2e1 * (-t158 * t257 - t163 * t254) / t137 ^ 2;
t227 = -0.2e1 * t256;
t203 = -t175 + t232;
t204 = t175 * t186 - qJD(1);
t143 = -t204 * t235 + (t203 * t185 + t209) * t173;
t157 = t262 ^ 2;
t155 = t157 * t243 + 0.1e1;
t182 = t180 * t181;
t226 = 0.2e1 * (-t262 * t143 * t243 + (t177 * t182 * t233 - t181 * t211) * t157) / t155 ^ 2;
t225 = t141 * t258;
t224 = t176 * t256;
t223 = t140 * t252;
t220 = t160 * t247;
t216 = t180 * t242;
t208 = t139 * t228;
t207 = t140 * t228;
t206 = t176 * t226;
t202 = t170 * t224;
t200 = t160 * t171 * t173 - t159 * t170;
t199 = -t159 * t180 + t241 * t262;
t153 = 0.1e1 / t155;
t145 = t204 * t238 + (t203 * t187 - t210) * t173;
t135 = 0.1e1 / t137;
t134 = t200 * t176 * t149;
t130 = (-t147 + (t148 * t220 + t147) * t149) * t163;
t129 = -t138 * t250 + (t234 * t251 + t248) * t174;
t128 = t163 * t180 * t206 + (t144 * t176 * t180 + (-t176 * t181 * t233 + t180 * t212) * t163) * t153;
t127 = -t148 * t244 + t147 * t159 - (-t147 * t240 - t250) * t134;
t125 = t201 * t227 + (t146 * t218 + t231 + (t171 * t186 * t217 - t170 * t261) * t160) * t149;
t123 = 0.2e1 * t200 * t224 + (t200 * t212 + ((-t160 * t175 + t145) * t170 + (t172 * t245 * t259 + (-t146 + t264) * t171) * t173) * t176) * t149;
t122 = (t127 * t255 + t139 * t262) * t228 + (t127 * t254 + t143 * t139 + (t127 * t225 + t140 * t262) * t126 - (-t173 * t230 - t175 * t240 - t123 * t160 + t134 * t146 + (t134 * t240 + t159) * t131) * t140 * t249 - (t145 + (-t123 * t174 + t131 * t173) * t184 + t196 * t134) * t223) * t135;
t1 = [t260 * t149 + t202 * t258, t125, t123, t123, 0, 0; t160 * t208 + (-t146 * t139 + (t126 * t160 + t130 * t144) * t140) * t135 + (t130 * t207 + (0.2e1 * t130 * t257 + (t144 * t149 - t144 - (-t131 * t149 * t220 + t227) * t163) * t140 * t147 + (-(t202 * t259 - t131) * t255 + (-(t131 + t221) * t163 + t260 * t160) * t140 * t149) * t148) * t135) * t163, t129 * t163 * t207 + (-(-t125 * t250 + (t131 * t253 - t146 * t148) * t138) * t255 + (t139 * t263 - (t138 * t251 - t147 * t239 - t248) * t255) * t245 + (t126 * t225 + t254) * t129) * t135 + (t208 * t263 + ((-t139 * t229 - (t234 * qJD(2) - t131) * t223) * t186 + (t139 * t233 + (t187 * t126 - (-t125 + t231) * t252 - (t234 * t131 - qJD(2)) * t249) * t140) * t184) * t135) * t174, t122, t122, 0, 0; t199 * t206 + (t199 * t212 + (t143 * t241 + t145 * t180 + (t159 * t241 - (0.2e1 * t182 * t185 ^ 2 + t180) * t262) * qJD(1)) * t176) * t153 (-t216 * t262 - t173) * t226 + (-t143 * t216 + t174 * t175 - (t180 * t261 - t213 * t243) * t262) * t153, t128, t128, 0, 0;];
JaD_rot  = t1;
