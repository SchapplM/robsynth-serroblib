% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:35:57
% EndTime: 2019-02-26 20:35:58
% DurationCPUTime: 0.71s
% Computational Cost: add. (2446->95), mult. (2734->206), div. (498->12), fcn. (3199->9), ass. (0->97)
t155 = sin(qJ(4));
t150 = 0.1e1 / t155 ^ 2;
t156 = cos(qJ(4));
t153 = t156 ^ 2;
t193 = t150 * t153;
t148 = qJ(1) + pkin(10);
t144 = cos(t148);
t217 = 0.2e1 * t144;
t216 = t156 * t193;
t197 = t144 * t156;
t136 = atan2(-t197, t155);
t134 = sin(t136);
t135 = cos(t136);
t124 = -t134 * t197 + t135 * t155;
t121 = 0.1e1 / t124;
t143 = sin(t148);
t154 = qJ(5) + qJ(6);
t145 = sin(t154);
t146 = cos(t154);
t195 = t146 * t155;
t131 = t143 * t195 + t144 * t145;
t127 = 0.1e1 / t131;
t149 = 0.1e1 / t155;
t122 = 0.1e1 / t124 ^ 2;
t128 = 0.1e1 / t131 ^ 2;
t142 = t144 ^ 2;
t139 = t142 * t193 + 0.1e1;
t137 = 0.1e1 / t139;
t215 = t137 - 0.1e1;
t141 = t143 ^ 2;
t199 = t141 * t153;
t116 = t122 * t199 + 0.1e1;
t191 = qJD(1) * t144;
t171 = t143 * t153 * t191;
t187 = qJD(4) * t156;
t190 = qJD(1) * t156;
t180 = t143 * t190;
t188 = qJD(4) * t155;
t189 = qJD(4) * t144;
t113 = ((t144 * t188 + t180) * t149 + t189 * t193) * t137;
t200 = t135 * t156;
t107 = (-t113 * t144 + qJD(4)) * t200 + (t180 + (-t113 + t189) * t155) * t134;
t212 = t107 * t121 * t122;
t214 = (-t199 * t212 + (-t141 * t155 * t187 + t171) * t122) / t116 ^ 2;
t147 = qJD(5) + qJD(6);
t175 = t147 * t155 + qJD(1);
t165 = t145 * t187 + t175 * t146;
t174 = qJD(1) * t155 + t147;
t168 = t144 * t174;
t111 = t165 * t143 + t145 * t168;
t196 = t145 * t155;
t130 = t143 * t196 - t144 * t146;
t126 = t130 ^ 2;
t119 = t126 * t128 + 0.1e1;
t203 = t128 * t130;
t164 = -t175 * t145 + t146 * t187;
t112 = t164 * t143 + t146 * t168;
t211 = t112 * t127 * t128;
t213 = (t111 * t203 - t126 * t211) / t119 ^ 2;
t210 = t113 * t134;
t209 = t113 * t156;
t166 = qJD(4) * (-t156 - t216) * t149;
t208 = (t142 * t166 - t150 * t171) / t139 ^ 2;
t207 = t122 * t143;
t206 = t122 * t156;
t178 = 0.1e1 + t193;
t125 = t178 * t144 * t137;
t205 = t125 * t144;
t204 = t127 * t145;
t202 = t130 * t146;
t201 = t134 * t155;
t198 = t143 * t156;
t194 = t149 * t153;
t192 = qJD(1) * t143;
t186 = 0.2e1 * t213;
t185 = -0.2e1 * t212;
t184 = t156 * t214;
t183 = t156 * t208;
t182 = t122 * t198;
t181 = t137 * t194;
t179 = t144 * t190;
t177 = 0.2e1 * t130 * t211;
t176 = t208 * t217;
t173 = t144 * t181;
t172 = t215 * t156 * t134;
t170 = t178 * t143;
t169 = t143 * t174;
t167 = t128 * t202 - t204;
t133 = -t143 * t145 + t144 * t195;
t132 = t143 * t146 + t144 * t196;
t117 = 0.1e1 / t119;
t114 = 0.1e1 / t116;
t110 = (-t135 * t173 - t172) * t143;
t109 = t144 * t201 + t200 + (-t135 * t197 - t201) * t125;
t108 = -t178 * t176 + (-qJD(1) * t170 + t166 * t217) * t137;
t104 = -0.2e1 * t213 + 0.2e1 * (t111 * t117 * t128 + (-t117 * t211 - t128 * t213) * t130) * t130;
t1 = [-0.2e1 * t143 * t149 * t183 + (-qJD(4) * t170 + t149 * t179) * t137, 0, 0, t108, 0, 0; (0.2e1 * t121 * t184 + (t121 * t188 + (qJD(1) * t110 + t107) * t206) * t114) * t144 + (-0.2e1 * t122 * t184 * t110 + (((t113 * t173 + t215 * t188 + 0.2e1 * t183) * t134 + (t176 * t194 + t209 + (-t209 + (0.2e1 * t156 + t216) * t189) * t137) * t135) * t182 + (-t122 * t188 + t156 * t185) * t110 + (t121 + ((t141 - t142) * t135 * t181 - t144 * t172) * t122) * t190) * t114) * t143, 0, 0, 0.2e1 * (-t109 * t206 - t121 * t155) * t143 * t214 + ((t121 * t191 + (-qJD(4) * t109 - t107) * t207) * t155 + (t143 * qJD(4) * t121 + (-t108 * t135 * t144 + t134 * t189 + t205 * t210 - t210 + (-qJD(4) * t134 + t135 * t192) * t125) * t182 + (t122 * t191 + t143 * t185) * t109 + ((-t108 - t192) * t134 + ((-0.1e1 + t205) * qJD(4) + (-t125 + t144) * t113) * t135) * t155 * t207) * t156) * t114, 0, 0; (-t127 * t132 + t133 * t203) * t186 + (t133 * t177 - t169 * t204 + t165 * t127 * t144 + (-t164 * t130 * t144 - t133 * t111 - t132 * t112 + t169 * t202) * t128) * t117, 0, 0, t167 * t186 * t198 + (-t167 * t179 + (t167 * t188 + ((t127 * t147 + t177) * t146 + (-t111 * t146 + (t130 * t147 - t112) * t145) * t128) * t156) * t143) * t117, t104, t104;];
JaD_rot  = t1;
