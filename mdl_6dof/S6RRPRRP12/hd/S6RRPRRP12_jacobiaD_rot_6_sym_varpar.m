% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP12
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
% Datum: 2019-02-26 21:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:13
% EndTime: 2019-02-26 21:52:14
% DurationCPUTime: 1.15s
% Computational Cost: add. (6985->119), mult. (8378->264), div. (1558->15), fcn. (10537->9), ass. (0->114)
t180 = qJ(4) + qJ(5);
t172 = sin(t180);
t181 = sin(qJ(2));
t184 = cos(qJ(1));
t173 = cos(t180);
t182 = sin(qJ(1));
t230 = t182 * t173;
t161 = t172 * t184 + t181 * t230;
t183 = cos(qJ(2));
t229 = t183 * t173;
t149 = atan2(t161, t229);
t145 = sin(t149);
t146 = cos(t149);
t158 = t161 ^ 2;
t170 = 0.1e1 / t173 ^ 2;
t177 = 0.1e1 / t183 ^ 2;
t237 = t170 * t177;
t150 = t158 * t237 + 0.1e1;
t147 = 0.1e1 / t150;
t169 = 0.1e1 / t173;
t211 = t169 * t177 * t181;
t198 = t161 * t211 + t182;
t136 = t198 * t147;
t227 = t136 - t182;
t255 = t227 * t183 * t145 + t146 * t181;
t231 = t181 * t184;
t160 = t172 * t231 + t230;
t156 = 0.1e1 / t160 ^ 2;
t175 = t183 ^ 2;
t179 = t184 ^ 2;
t232 = t175 * t179;
t214 = t156 * t232;
t153 = 0.1e1 + t214;
t223 = qJD(2) * t183;
t225 = qJD(1) * t184;
t194 = -t175 * t182 * t225 - t179 * t181 * t223;
t174 = qJD(4) + qJD(5);
t202 = qJD(1) * t181 + t174;
t203 = t174 * t181 + qJD(1);
t222 = qJD(2) * t184;
t206 = t183 * t222;
t233 = t173 * t184;
t144 = t203 * t233 + (-t182 * t202 + t206) * t172;
t155 = 0.1e1 / t160;
t245 = t144 * t155 * t156;
t201 = t232 * t245;
t254 = (t156 * t194 - t201) / t153 ^ 2;
t176 = 0.1e1 / t183;
t235 = t172 * t182;
t162 = -t181 * t235 + t233;
t196 = t161 * t170 * t172 + t162 * t169;
t253 = t176 * t196;
t228 = t183 * t184;
t243 = t145 * t161;
t140 = t146 * t229 + t243;
t137 = 0.1e1 / t140;
t138 = 0.1e1 / t140 ^ 2;
t209 = t173 * t231;
t159 = -t209 + t235;
t251 = 0.2e1 * t159;
t154 = t159 ^ 2;
t135 = t138 * t154 + 0.1e1;
t143 = qJD(1) * t161 + t160 * t174 - t173 * t206;
t246 = t143 * t138;
t224 = qJD(2) * t181;
t234 = t172 * t183;
t195 = -t173 * t224 - t174 * t234;
t212 = t161 * t237;
t199 = t203 * t182;
t207 = t182 * t223;
t141 = -qJD(1) * t209 + t172 * t199 - t173 * t207 - t174 * t233;
t238 = t169 * t176;
t215 = t141 * t238;
t129 = (-t195 * t212 - t215) * t147;
t193 = -t129 * t161 - t195;
t124 = (-t129 * t229 - t141) * t145 - t193 * t146;
t139 = t137 * t138;
t249 = t124 * t139;
t250 = (-t154 * t249 + t159 * t246) / t135 ^ 2;
t171 = t169 * t170;
t178 = t176 / t175;
t236 = t172 * t174;
t210 = t177 * t236;
t248 = (-t141 * t212 + (t170 * t178 * t224 + t171 * t210) * t158) / t150 ^ 2;
t247 = t138 * t159;
t244 = t145 * t159;
t241 = t146 * t159;
t240 = t146 * t161;
t226 = qJD(1) * t182;
t221 = 0.2e1 * t250;
t220 = -0.2e1 * t248;
t219 = 0.2e1 * t254;
t218 = t139 * t251;
t217 = t137 * t250;
t216 = t138 * t244;
t213 = t161 * t238;
t208 = t177 * t224;
t205 = t138 * t221;
t204 = t228 * t251;
t200 = 0.2e1 * t238 * t248;
t197 = -t156 * t162 * t184 - t155 * t182;
t192 = t143 * t238 - (-t170 * t176 * t236 - t169 * t208) * t159;
t151 = 0.1e1 / t153;
t142 = t173 * t199 + (t184 * t202 + t207) * t172;
t133 = 0.1e1 / t135;
t132 = t147 * t253;
t128 = (-t145 + (-t146 * t213 + t145) * t147) * t159;
t127 = t136 * t240 - t255 * t173;
t126 = -t146 * t234 + t145 * t162 + (-t145 * t229 + t240) * t132;
t125 = t156 * t204 * t254 + (t204 * t245 + (-t143 * t228 + (t181 * t222 + t183 * t226) * t159) * t156) * t151;
t123 = t198 * t220 + (-t141 * t211 + t225 + (t170 * t181 * t210 + (0.2e1 * t178 * t181 ^ 2 + t176) * t169 * qJD(2)) * t161) * t147;
t121 = t220 * t253 + (t196 * t208 + ((t161 * t174 - t142) * t169 + (0.2e1 * t161 * t171 * t236 + (t162 * t174 - t141) * t170) * t172) * t176) * t147;
t120 = (t126 * t247 - t137 * t160) * t221 + (-t126 * t246 + t144 * t137 + (t126 * t218 - t138 * t160) * t124 - (t172 * t224 - t174 * t229 + t121 * t161 - t132 * t141 + (-t132 * t229 + t162) * t129) * t138 * t241 - (-t142 + (-t121 * t173 + t129 * t172) * t183 + t193 * t132) * t216) * t133;
t1 = [-t147 * t192 + t159 * t200, t123, 0, t121, t121, 0; -0.2e1 * t161 * t217 + (-t141 * t137 + (-t124 * t161 - t128 * t143) * t138) * t133 + (t128 * t205 + (0.2e1 * t128 * t249 + (-t143 * t147 + t143 - (t129 * t147 * t213 + t220) * t159) * t138 * t145 + (-(t161 * t200 - t129) * t247 + (-(t129 + t215) * t159 + t192 * t161) * t138 * t147) * t146) * t133) * t159, t127 * t159 * t205 + (-(t123 * t240 + (-t129 * t243 - t141 * t146) * t136) * t247 + (t137 * t228 - t255 * t247) * t236 + (t124 * t218 - t246) * t127) * t133 + (0.2e1 * t217 * t228 + ((t137 * t222 - (t227 * qJD(2) + t129) * t216) * t181 + (t137 * t226 + (t184 * t124 - (-t123 + t225) * t244 - (-t227 * t129 - qJD(2)) * t241) * t138) * t183) * t133) * t173, 0, t120, t120, 0; t197 * t183 * t219 + (t197 * t224 + ((qJD(1) * t155 - 0.2e1 * t162 * t245) * t184 + (-t142 * t184 + (-qJD(1) * t162 - t144) * t182) * t156) * t183) * t151 (-t155 * t231 - t172 * t214) * t219 + (-0.2e1 * t172 * t201 + (-t181 * t226 + t206) * t155 + (t173 * t174 * t232 - t144 * t231 + 0.2e1 * t172 * t194) * t156) * t151, 0, t125, t125, 0;];
JaD_rot  = t1;
