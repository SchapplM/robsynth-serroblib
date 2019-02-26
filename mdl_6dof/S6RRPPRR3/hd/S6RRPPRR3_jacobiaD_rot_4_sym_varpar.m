% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiaD_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:33
% EndTime: 2019-02-26 21:29:34
% DurationCPUTime: 0.82s
% Computational Cost: add. (2974->87), mult. (8937->192), div. (424->12), fcn. (11581->13), ass. (0->95)
t184 = sin(pkin(6));
t237 = sin(pkin(11));
t238 = cos(pkin(11));
t240 = sin(qJ(2));
t241 = cos(qJ(2));
t198 = t241 * t237 + t240 * t238;
t174 = t198 * t184;
t165 = qJD(2) * t174;
t177 = t240 * t237 - t241 * t238;
t173 = t177 * t184;
t170 = 0.1e1 / t173 ^ 2;
t225 = t165 * t170;
t239 = cos(pkin(6));
t206 = t239 * t237;
t207 = t239 * t238;
t243 = t240 * t206 - t241 * t207;
t197 = t177 * qJD(2);
t186 = sin(qJ(1));
t187 = cos(qJ(1));
t157 = -t186 * t198 - t187 * t243;
t145 = atan2(t157, t173);
t140 = sin(t145);
t141 = cos(t145);
t154 = t157 ^ 2;
t144 = t154 * t170 + 0.1e1;
t142 = 0.1e1 / t144;
t169 = 0.1e1 / t173;
t227 = t157 * t169;
t242 = t142 * (t141 * t227 - t140) + t140;
t129 = t140 * t157 + t141 * t173;
t126 = 0.1e1 / t129;
t183 = sin(pkin(12));
t185 = cos(pkin(12));
t175 = t241 * t206 + t240 * t207;
t203 = t186 * t175 + t187 * t177;
t222 = t184 * t186;
t153 = t183 * t222 - t185 * t203;
t147 = 0.1e1 / t153;
t127 = 0.1e1 / t129 ^ 2;
t148 = 0.1e1 / t153 ^ 2;
t194 = t186 * t243;
t160 = -t187 * t198 + t194;
t155 = t160 ^ 2;
t125 = t155 * t127 + 0.1e1;
t168 = t175 * qJD(2);
t135 = qJD(1) * t157 - t186 * t168 - t187 * t197;
t231 = t135 * t127;
t219 = qJD(1) * t187;
t138 = qJD(1) * t194 - t187 * t168 + t186 * t197 - t198 * t219;
t202 = t138 * t169 - t157 * t225;
t120 = t202 * t142;
t205 = -t140 * t173 + t141 * t157;
t116 = t120 * t205 + t140 * t138 + t141 * t165;
t235 = t116 * t126 * t127;
t236 = (-t155 * t235 - t160 * t231) / t125 ^ 2;
t204 = -t187 * t175 + t186 * t177;
t226 = t157 * t174;
t201 = -t169 * t204 + t170 * t226;
t121 = t201 * t142;
t117 = -t121 * t205 + t140 * t204 + t141 * t174;
t234 = t117 * t160;
t224 = t169 * t225;
t233 = (t157 * t170 * t138 - t154 * t224) / t144 ^ 2;
t167 = t243 * qJD(2);
t176 = t198 * qJD(2);
t136 = qJD(1) * t204 + t186 * t167 - t187 * t176;
t213 = t184 * t219;
t134 = t136 * t185 + t183 * t213;
t232 = t134 * t147 * t148;
t230 = t140 * t160;
t229 = t141 * t160;
t152 = -t183 * t203 - t185 * t222;
t228 = t148 * t152;
t223 = t183 * t147;
t221 = t184 * t187;
t220 = t185 * t152;
t218 = -0.2e1 * t236;
t217 = -0.2e1 * t235;
t146 = t152 ^ 2;
t132 = t146 * t148 + 0.1e1;
t133 = t136 * t183 - t185 * t213;
t216 = 0.2e1 * (t133 * t228 - t146 * t232) / t132 ^ 2;
t215 = 0.2e1 * t233;
t214 = qJD(1) * t222;
t212 = -0.2e1 * t169 * t233;
t211 = 0.2e1 * t152 * t232;
t196 = qJD(1) * t203 + t187 * t167 + t186 * t176;
t166 = t184 * t197;
t151 = t183 * t221 + t185 * t204;
t150 = t183 * t204 - t185 * t221;
t130 = 0.1e1 / t132;
t123 = 0.1e1 / t125;
t118 = t242 * t160;
t115 = t201 * t215 + (0.2e1 * t224 * t226 + t196 * t169 + (-t138 * t174 + t157 * t166 - t165 * t204) * t170) * t142;
t1 = [t160 * t212 + (-t135 * t169 - t160 * t225) * t142, t115, 0, 0, 0, 0; t157 * t126 * t218 + (t138 * t126 + (-t116 * t157 - t118 * t135) * t127) * t123 + ((t118 * t217 - t242 * t231) * t123 + (t118 * t218 + ((-t120 * t142 * t227 + t215) * t230 + (t157 * t212 + t120 + (-t120 + t202) * t142) * t229) * t123) * t127) * t160, 0.2e1 * (t126 * t203 - t127 * t234) * t236 + (t136 * t126 + t217 * t234 + (t203 * t116 - t117 * t135 + (t115 * t157 - t121 * t138 - t166 + (t121 * t173 + t204) * t120) * t229 + (-t115 * t173 + t121 * t165 + t196 + (t121 * t157 - t174) * t120) * t230) * t127) * t123, 0, 0, 0, 0; (-t147 * t150 + t151 * t228) * t216 + ((t183 * t196 + t185 * t214) * t147 + t151 * t211 + (-t150 * t134 - (-t183 * t214 + t185 * t196) * t152 - t151 * t133) * t148) * t130 (t148 * t220 - t223) * t160 * t216 + (t160 * t185 * t211 - t135 * t223 + (t135 * t220 + (-t133 * t185 - t134 * t183) * t160) * t148) * t130, 0, 0, 0, 0;];
JaD_rot  = t1;
