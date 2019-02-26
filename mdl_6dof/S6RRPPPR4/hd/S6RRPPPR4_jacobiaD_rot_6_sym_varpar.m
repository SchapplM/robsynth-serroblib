% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:45
% EndTime: 2019-02-26 21:23:46
% DurationCPUTime: 0.83s
% Computational Cost: add. (1115->105), mult. (3345->234), div. (468->12), fcn. (4023->11), ass. (0->103)
t158 = sin(qJ(2));
t149 = 0.1e1 / t158 ^ 2;
t161 = cos(qJ(2));
t153 = t161 ^ 2;
t211 = t149 * t153;
t227 = t161 * t211;
t159 = sin(qJ(1));
t151 = t159 ^ 2;
t147 = t151 * t211 + 0.1e1;
t145 = 0.1e1 / t147;
t148 = 0.1e1 / t158;
t162 = cos(qJ(1));
t202 = qJD(1) * t162;
t187 = t161 * t202;
t200 = qJD(2) * t159;
t113 = (-(-t158 * t200 + t187) * t148 + t200 * t211) * t145;
t226 = t113 - t200;
t157 = sin(qJ(6));
t160 = cos(qJ(6));
t197 = qJD(6) * t162;
t204 = qJD(1) * t159;
t225 = -t157 * t204 + t160 * t197;
t224 = -t157 * t197 - t160 * t204;
t155 = sin(pkin(9));
t156 = cos(pkin(9));
t207 = t159 * t156;
t209 = t158 * t162;
t139 = t155 * t209 + t207;
t208 = t159 * t155;
t171 = t156 * t209 - t208;
t122 = t139 * t160 - t157 * t171;
t118 = 0.1e1 / t122;
t206 = t159 * t161;
t144 = atan2(t206, -t158);
t143 = cos(t144);
t142 = sin(t144);
t191 = t142 * t206;
t128 = -t143 * t158 + t191;
t125 = 0.1e1 / t128;
t119 = 0.1e1 / t122 ^ 2;
t126 = 0.1e1 / t128 ^ 2;
t223 = t145 - 0.1e1;
t154 = t162 ^ 2;
t210 = t153 * t154;
t116 = t126 * t210 + 0.1e1;
t177 = t153 * t159 * t202;
t199 = qJD(2) * t161;
t213 = t143 * t161;
t104 = (t113 * t159 - qJD(2)) * t213 + (t226 * t158 + t187) * t142;
t221 = t104 * t125 * t126;
t222 = (-t210 * t221 + (-t154 * t158 * t199 - t177) * t126) / t116 ^ 2;
t140 = t155 * t162 + t158 * t207;
t198 = qJD(2) * t162;
t185 = t161 * t198;
t133 = t140 * qJD(1) - t156 * t185;
t141 = t156 * t162 - t158 * t208;
t134 = t141 * qJD(1) + t155 * t185;
t175 = -t139 * t157 - t160 * t171;
t108 = t175 * qJD(6) + t133 * t157 + t134 * t160;
t120 = t118 * t119;
t220 = t108 * t120;
t107 = t122 * qJD(6) - t133 * t160 + t134 * t157;
t117 = t175 ^ 2;
t112 = t117 * t119 + 0.1e1;
t217 = t119 * t175;
t219 = 0.1e1 / t112 ^ 2 * (-t107 * t217 - t117 * t220);
t218 = t113 * t161;
t173 = t155 * t160 - t156 * t157;
t205 = t161 * t162;
t136 = t173 * t205;
t216 = t119 * t136;
t215 = t126 * t161;
t214 = t142 * t159;
t212 = t148 * t153;
t203 = qJD(1) * t161;
t201 = qJD(2) * t158;
t196 = -0.2e1 * t221;
t195 = 0.2e1 * t219;
t194 = -0.2e1 * t120 * t175;
t170 = qJD(2) * (-t161 - t227) * t148;
t193 = 0.2e1 * (t149 * t177 + t151 * t170) / t147 ^ 2;
t192 = t126 * t205;
t190 = t145 * t212;
t186 = t159 * t199;
t182 = 0.1e1 + t211;
t181 = -0.2e1 * t161 * t222;
t180 = t159 * t193;
t179 = t161 * t193;
t178 = t159 * t190;
t176 = t182 * t162;
t174 = t140 * t160 - t141 * t157;
t124 = t140 * t157 + t141 * t160;
t172 = t155 * t157 + t156 * t160;
t135 = t172 * t205;
t132 = -t139 * qJD(1) - t155 * t186;
t131 = t171 * qJD(1) + t156 * t186;
t130 = t182 * t159 * t145;
t114 = 0.1e1 / t116;
t111 = (-t223 * t161 * t142 - t143 * t178) * t162;
t109 = 0.1e1 / t112;
t106 = -t158 * t214 - t213 + (t142 * t158 + t143 * t206) * t130;
t105 = -t182 * t180 + (qJD(1) * t176 + 0.2e1 * t159 * t170) * t145;
t1 = [t162 * t148 * t179 + (t148 * t159 * t203 + qJD(2) * t176) * t145, t105, 0, 0, 0, 0; (t125 * t181 + (-t125 * t201 + (-qJD(1) * t111 - t104) * t215) * t114) * t159 + (t126 * t181 * t111 + (((t113 * t178 + t223 * t201 + t179) * t142 + (t180 * t212 + t218 + (-t218 + (0.2e1 * t161 + t227) * t200) * t145) * t143) * t192 + (-t126 * t201 + t161 * t196) * t111 + (t125 + ((t151 - t154) * t143 * t190 + t223 * t191) * t126) * t203) * t114) * t162, 0.2e1 * (-t106 * t215 - t125 * t158) * t162 * t222 + ((-t125 * t204 + (-qJD(2) * t106 - t104) * t162 * t126) * t158 + (t125 * t198 + (t105 * t143 * t159 + t226 * t142 + (qJD(2) * t142 - t113 * t214 + t143 * t202) * t130) * t192 + (-t126 * t204 + t162 * t196) * t106 + ((t105 - t202) * t142 + ((-t130 * t159 + 0.1e1) * qJD(2) + (t130 - t159) * t113) * t143) * t126 * t209) * t161) * t114, 0, 0, 0, 0; (t118 * t174 - t124 * t217) * t195 + ((t124 * qJD(6) - t131 * t160 + t132 * t157) * t118 + t124 * t108 * t194 + (t174 * t108 + (t174 * qJD(6) + t131 * t157 + t132 * t160) * t175 - t124 * t107) * t119) * t109 (-t118 * t135 - t175 * t216) * t195 + (-t107 * t216 + (-t119 * t135 + t136 * t194) * t108 + (-t172 * t118 - t173 * t217) * t158 * t198 + ((t224 * t118 - t225 * t217) * t156 + (t225 * t118 + t224 * t217) * t155) * t161) * t109, 0, 0, 0, -0.2e1 * t219 - 0.2e1 * (t107 * t119 * t109 - (-t109 * t220 - t119 * t219) * t175) * t175;];
JaD_rot  = t1;
