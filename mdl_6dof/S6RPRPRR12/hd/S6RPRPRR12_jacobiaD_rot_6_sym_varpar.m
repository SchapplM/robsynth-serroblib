% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:02
% EndTime: 2019-02-26 20:55:03
% DurationCPUTime: 0.71s
% Computational Cost: add. (1192->94), mult. (2734->207), div. (498->12), fcn. (3199->9), ass. (0->95)
t155 = sin(qJ(3));
t148 = t155 ^ 2;
t157 = cos(qJ(3));
t151 = 0.1e1 / t157 ^ 2;
t199 = t148 * t151;
t158 = cos(qJ(1));
t219 = 0.2e1 * t158;
t218 = t155 * t199;
t153 = t158 ^ 2;
t143 = t153 * t199 + 0.1e1;
t140 = 0.1e1 / t143;
t150 = 0.1e1 / t157;
t156 = sin(qJ(1));
t191 = qJD(1) * t156;
t179 = t155 * t191;
t187 = qJD(3) * t158;
t115 = ((t157 * t187 - t179) * t150 + t187 * t199) * t140;
t217 = t115 - t187;
t146 = qJD(5) + qJD(6);
t173 = qJD(1) * t157 + t146;
t216 = t155 * t187 + t173 * t156;
t192 = t158 * t155;
t142 = atan2(t192, t157);
t137 = sin(t142);
t138 = cos(t142);
t125 = t137 * t192 + t138 * t157;
t122 = 0.1e1 / t125;
t154 = qJ(5) + qJ(6);
t144 = sin(t154);
t145 = cos(t154);
t193 = t158 * t145;
t195 = t156 * t157;
t135 = -t144 * t195 + t193;
t129 = 0.1e1 / t135;
t123 = 0.1e1 / t125 ^ 2;
t130 = 0.1e1 / t135 ^ 2;
t215 = 0.2e1 * t156;
t214 = t140 - 0.1e1;
t149 = t156 ^ 2;
t201 = t148 * t149;
t120 = t123 * t201 + 0.1e1;
t190 = qJD(1) * t158;
t170 = t148 * t156 * t190;
t188 = qJD(3) * t157;
t202 = t138 * t155;
t109 = (t115 * t158 - qJD(3)) * t202 + (-t217 * t157 - t179) * t137;
t212 = t109 * t122 * t123;
t213 = (-t201 * t212 + (t149 * t155 * t188 + t170) * t123) / t120 ^ 2;
t174 = t146 * t157 + qJD(1);
t189 = qJD(3) * t156;
t196 = t156 * t145;
t114 = -t174 * t196 + (t155 * t189 - t173 * t158) * t144;
t211 = t114 * t129 * t130;
t210 = t115 * t155;
t180 = t157 * t193;
t113 = -qJD(1) * t180 - t146 * t193 + (qJD(3) * t145 * t155 + t174 * t144) * t156;
t194 = t158 * t144;
t134 = t145 * t195 + t194;
t128 = t134 ^ 2;
t121 = t128 * t130 + 0.1e1;
t205 = t130 * t134;
t209 = 0.1e1 / t121 ^ 2 * (-t113 * t205 - t128 * t211);
t208 = t123 * t155;
t198 = t150 * t155;
t166 = qJD(3) * (t150 * t218 + t198);
t207 = (-t151 * t170 + t153 * t166) / t143 ^ 2;
t206 = t129 * t145;
t204 = t134 * t144;
t203 = t137 * t158;
t200 = t148 * t150;
t197 = t155 * t156;
t186 = 0.2e1 * t212;
t185 = 0.2e1 * t209;
t184 = t155 * t213;
t183 = t123 * t197;
t182 = t155 * t207;
t181 = t140 * t200;
t177 = 0.1e1 + t199;
t176 = 0.2e1 * t134 * t211;
t175 = t207 * t219;
t172 = t158 * t181;
t171 = t214 * t155 * t137;
t169 = t177 * t156;
t168 = t174 * t158;
t167 = t130 * t204 + t206;
t133 = -t157 * t194 - t196;
t132 = -t156 * t144 + t180;
t127 = t177 * t158 * t140;
t118 = 0.1e1 / t121;
t116 = 0.1e1 / t120;
t112 = (-t138 * t172 + t171) * t156;
t111 = t157 * t203 - t202 + (-t137 * t157 + t138 * t192) * t127;
t110 = -t177 * t175 + (-qJD(1) * t169 + t166 * t219) * t140;
t106 = -0.2e1 * t209 + 0.2e1 * (-t113 * t130 * t118 + (-t118 * t211 - t130 * t209) * t134) * t134;
t1 = [t150 * t182 * t215 + (-qJD(3) * t169 - t190 * t198) * t140, 0, t110, 0, 0, 0; (-0.2e1 * t122 * t184 + (t122 * t188 + (-qJD(1) * t112 - t109) * t208) * t116) * t158 + (0.2e1 * t123 * t184 * t112 + (-((t115 * t172 + t214 * t188 - 0.2e1 * t182) * t137 + (t175 * t200 - t210 + (t210 + (-0.2e1 * t155 - t218) * t187) * t140) * t138) * t183 + (-t123 * t188 + t155 * t186) * t112 + (-t122 + ((-t149 + t153) * t138 * t181 - t158 * t171) * t123) * t155 * qJD(1)) * t116) * t156, 0 (t111 * t208 - t122 * t157) * t213 * t215 + ((t122 * t190 + (-qJD(3) * t111 - t109) * t156 * t123) * t157 + (-t122 * t189 - (t110 * t138 * t158 + t217 * t137 + (qJD(3) * t137 - t115 * t203 - t138 * t191) * t127) * t183 + (-t123 * t190 + t156 * t186) * t111 - ((-t110 - t191) * t137 + ((t127 * t158 - 0.1e1) * qJD(3) + (-t127 + t158) * t115) * t138) * t123 * t195) * t155) * t116, 0, 0, 0; (-t129 * t132 + t133 * t205) * t185 + (t133 * t176 - t129 * t144 * t168 - t216 * t206 + (t134 * t145 * t168 + t133 * t113 - t132 * t114 - t216 * t204) * t130) * t118, 0, t167 * t185 * t197 + (-t167 * t156 * t188 + (-t167 * t190 + ((t129 * t146 + t176) * t144 + (t113 * t144 + (-t134 * t146 + t114) * t145) * t130) * t156) * t155) * t118, 0, t106, t106;];
JaD_rot  = t1;
