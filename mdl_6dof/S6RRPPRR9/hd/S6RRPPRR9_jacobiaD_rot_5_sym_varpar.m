% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:13
% EndTime: 2019-02-26 21:33:14
% DurationCPUTime: 0.60s
% Computational Cost: add. (1334->89), mult. (4303->202), div. (668->14), fcn. (5516->11), ass. (0->93)
t175 = sin(qJ(2));
t176 = sin(qJ(1));
t178 = cos(qJ(2));
t179 = cos(qJ(1));
t226 = cos(pkin(6));
t195 = t179 * t226;
t158 = t175 * t195 + t176 * t178;
t196 = t176 * t226;
t159 = t179 * t175 + t178 * t196;
t138 = t159 * qJD(1) + t158 * qJD(2);
t192 = t178 * t195;
t211 = t176 * t175;
t157 = -t192 + t211;
t155 = t157 ^ 2;
t173 = sin(pkin(6));
t169 = 0.1e1 / t173 ^ 2;
t171 = 0.1e1 / t178 ^ 2;
t153 = t155 * t169 * t171 + 0.1e1;
t170 = 0.1e1 / t178;
t172 = t170 * t171;
t208 = qJD(2) * t175;
t221 = (t138 * t157 * t171 + t155 * t172 * t208) * t169 / t153 ^ 2;
t229 = -0.2e1 * t221;
t168 = 0.1e1 / t173;
t228 = t157 * t168;
t213 = t173 * t178;
t152 = atan2(t157, t213);
t148 = sin(t152);
t149 = cos(t152);
t150 = 0.1e1 / t153;
t200 = t170 * t228;
t227 = (t149 * t200 - t148) * t150 + t148;
t132 = t148 * t157 + t149 * t213;
t129 = 0.1e1 / t132;
t193 = t175 * t196;
t210 = t179 * t178;
t161 = -t193 + t210;
t174 = sin(qJ(5));
t177 = cos(qJ(5));
t214 = t173 * t176;
t145 = t161 * t174 + t177 * t214;
t141 = 0.1e1 / t145;
t130 = 0.1e1 / t132 ^ 2;
t142 = 0.1e1 / t145 ^ 2;
t156 = t159 ^ 2;
t125 = t156 * t130 + 0.1e1;
t191 = qJD(2) * t226 + qJD(1);
t207 = qJD(2) * t178;
t136 = -qJD(1) * t192 - t179 * t207 + t191 * t211;
t219 = t136 * t130;
t197 = t171 * t208;
t186 = (t138 * t170 + t157 * t197) * t168;
t121 = t150 * t186;
t188 = -t148 * t213 + t149 * t157;
t201 = t149 * t173 * t175;
t117 = -qJD(2) * t201 + t188 * t121 + t148 * t138;
t224 = t117 * t129 * t130;
t225 = (-t156 * t224 - t159 * t219) / t125 ^ 2;
t215 = t171 * t175;
t187 = t157 * t215 + t158 * t170;
t122 = t187 * t168 * t150;
t118 = t188 * t122 + t148 * t158 - t201;
t223 = t118 * t159;
t137 = t158 * qJD(1) + t159 * qJD(2);
t209 = qJD(1) * t173;
t198 = t179 * t209;
t127 = t145 * qJD(5) + t137 * t177 + t174 * t198;
t144 = -t161 * t177 + t174 * t214;
t140 = t144 ^ 2;
t135 = t140 * t142 + 0.1e1;
t218 = t142 * t144;
t206 = qJD(5) * t144;
t128 = -t137 * t174 + t177 * t198 - t206;
t220 = t128 * t141 * t142;
t222 = (t127 * t218 - t140 * t220) / t135 ^ 2;
t217 = t148 * t159;
t216 = t149 * t159;
t212 = t173 * t179;
t205 = -0.2e1 * t225;
t204 = -0.2e1 * t224;
t203 = 0.2e1 * t222;
t202 = t144 * t220;
t199 = t176 * t209;
t194 = t170 * t229;
t189 = -t177 * t141 - t174 * t218;
t147 = -t158 * t174 + t177 * t212;
t146 = t158 * t177 + t174 * t212;
t139 = -qJD(1) * t193 - t176 * t208 + t191 * t210;
t133 = 0.1e1 / t135;
t123 = 0.1e1 / t125;
t120 = t227 * t159;
t116 = (t187 * t229 + (t138 * t215 + t139 * t170 + (t158 * t215 + (0.2e1 * t172 * t175 ^ 2 + t170) * t157) * qJD(2)) * t150) * t168;
t1 = [(t159 * t194 + (-t136 * t170 + t159 * t197) * t150) * t168, t116, 0, 0, 0, 0; t157 * t129 * t205 + (t138 * t129 + (-t117 * t157 - t120 * t136) * t130) * t123 + ((t120 * t204 - t227 * t219) * t123 + (t120 * t205 + ((-t121 * t150 * t200 + 0.2e1 * t221) * t217 + (t194 * t228 + t121 + (-t121 + t186) * t150) * t216) * t123) * t130) * t159, 0.2e1 * (t129 * t161 - t130 * t223) * t225 + (t204 * t223 + t137 * t129 + (t161 * t117 - t118 * t136 + (-t173 * t207 + t116 * t157 + t122 * t138 + (-t122 * t213 + t158) * t121) * t216 + (-t121 * t122 * t157 + t139 + (-t116 * t178 + (qJD(2) * t122 + t121) * t175) * t173) * t217) * t130) * t123, 0, 0, 0, 0; (-t141 * t146 + t147 * t218) * t203 + ((t147 * qJD(5) + t139 * t177 - t174 * t199) * t141 + 0.2e1 * t147 * t202 + (-t146 * t128 - (-t146 * qJD(5) - t139 * t174 - t177 * t199) * t144 - t147 * t127) * t142) * t133, t189 * t159 * t203 + (t189 * t136 + ((-qJD(5) * t141 - 0.2e1 * t202) * t174 + (t127 * t174 + (-t128 + t206) * t177) * t142) * t159) * t133, 0, 0, -0.2e1 * t222 + 0.2e1 * (t127 * t142 * t133 + (-t133 * t220 - t142 * t222) * t144) * t144, 0;];
JaD_rot  = t1;
