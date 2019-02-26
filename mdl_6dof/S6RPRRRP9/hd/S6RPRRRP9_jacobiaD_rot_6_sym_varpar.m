% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:12:30
% EndTime: 2019-02-26 21:12:31
% DurationCPUTime: 0.72s
% Computational Cost: add. (1381->94), mult. (2734->205), div. (498->12), fcn. (3199->9), ass. (0->97)
t152 = sin(qJ(3));
t145 = 0.1e1 / t152 ^ 2;
t154 = cos(qJ(3));
t149 = t154 ^ 2;
t199 = t145 * t149;
t155 = cos(qJ(1));
t217 = 0.2e1 * t155;
t216 = t154 * t199;
t192 = t155 * t154;
t136 = atan2(-t192, t152);
t134 = sin(t136);
t135 = cos(t136);
t123 = -t134 * t192 + t135 * t152;
t120 = 0.1e1 / t123;
t151 = qJ(4) + qJ(5);
t141 = sin(t151);
t142 = cos(t151);
t153 = sin(qJ(1));
t195 = t153 * t142;
t131 = t141 * t155 + t152 * t195;
t127 = 0.1e1 / t131;
t144 = 0.1e1 / t152;
t121 = 0.1e1 / t123 ^ 2;
t128 = 0.1e1 / t131 ^ 2;
t150 = t155 ^ 2;
t139 = t150 * t199 + 0.1e1;
t137 = 0.1e1 / t139;
t215 = t137 - 0.1e1;
t147 = t153 ^ 2;
t198 = t147 * t149;
t116 = t121 * t198 + 0.1e1;
t189 = qJD(1) * t155;
t170 = t149 * t153 * t189;
t187 = qJD(3) * t154;
t190 = qJD(1) * t154;
t179 = t153 * t190;
t186 = qJD(3) * t155;
t113 = ((t152 * t186 + t179) * t144 + t186 * t199) * t137;
t201 = t135 * t154;
t107 = (-t113 * t155 + qJD(3)) * t201 + (t179 + (-t113 + t186) * t152) * t134;
t212 = t107 * t120 * t121;
t214 = (-t198 * t212 + (-t147 * t152 * t187 + t170) * t121) / t116 ^ 2;
t143 = qJD(4) + qJD(5);
t173 = qJD(1) * t152 + t143;
t163 = t153 * t187 + t173 * t155;
t174 = t143 * t152 + qJD(1);
t167 = t142 * t174;
t111 = t163 * t141 + t153 * t167;
t193 = t155 * t142;
t196 = t153 * t141;
t130 = t152 * t196 - t193;
t126 = t130 ^ 2;
t119 = t126 * t128 + 0.1e1;
t203 = t128 * t130;
t168 = t141 * t174;
t112 = t163 * t142 - t153 * t168;
t211 = t112 * t127 * t128;
t213 = (t111 * t203 - t126 * t211) / t119 ^ 2;
t210 = t113 * t134;
t209 = t113 * t154;
t208 = t121 * t153;
t207 = t121 * t154;
t165 = qJD(3) * (-t154 - t216) * t144;
t206 = (-t145 * t170 + t150 * t165) / t139 ^ 2;
t177 = 0.1e1 + t199;
t125 = t177 * t155 * t137;
t205 = t125 * t155;
t204 = t127 * t141;
t202 = t130 * t142;
t200 = t144 * t149;
t197 = t152 * t155;
t194 = t153 * t154;
t191 = qJD(1) * t153;
t188 = qJD(3) * t152;
t185 = 0.2e1 * t213;
t184 = -0.2e1 * t212;
t183 = t154 * t214;
t182 = t121 * t194;
t181 = t154 * t206;
t180 = t137 * t200;
t178 = t154 * t189;
t176 = 0.2e1 * t130 * t211;
t175 = t206 * t217;
t172 = t155 * t180;
t171 = t215 * t154 * t134;
t169 = t177 * t153;
t166 = t128 * t202 - t204;
t164 = -t173 * t153 + t154 * t186;
t133 = t152 * t193 - t196;
t132 = t141 * t197 + t195;
t117 = 0.1e1 / t119;
t114 = 0.1e1 / t116;
t110 = (-t135 * t172 - t171) * t153;
t109 = t134 * t197 + t201 + (-t134 * t152 - t135 * t192) * t125;
t108 = -t177 * t175 + (-qJD(1) * t169 + t165 * t217) * t137;
t104 = -0.2e1 * t213 + 0.2e1 * (t111 * t117 * t128 + (-t117 * t211 - t128 * t213) * t130) * t130;
t1 = [-0.2e1 * t153 * t144 * t181 + (-qJD(3) * t169 + t144 * t178) * t137, 0, t108, 0, 0, 0; (0.2e1 * t120 * t183 + (t120 * t188 + (qJD(1) * t110 + t107) * t207) * t114) * t155 + (-0.2e1 * t121 * t183 * t110 + (((t113 * t172 + t215 * t188 + 0.2e1 * t181) * t134 + (t175 * t200 + t209 + (-t209 + (0.2e1 * t154 + t216) * t186) * t137) * t135) * t182 + (-t121 * t188 + t154 * t184) * t110 + (t120 + ((t147 - t150) * t135 * t180 - t155 * t171) * t121) * t190) * t114) * t153, 0, 0.2e1 * (-t109 * t207 - t120 * t152) * t153 * t214 + ((t120 * t189 + (-qJD(3) * t109 - t107) * t208) * t152 + (t153 * qJD(3) * t120 + (-t108 * t135 * t155 + t134 * t186 + t205 * t210 - t210 + (-qJD(3) * t134 + t135 * t191) * t125) * t182 + (t121 * t189 + t153 * t184) * t109 + ((-t108 - t191) * t134 + ((-0.1e1 + t205) * qJD(3) + (-t125 + t155) * t113) * t135) * t152 * t208) * t154) * t114, 0, 0, 0; (-t127 * t132 + t133 * t203) * t185 + (t133 * t176 + t155 * t127 * t167 + t164 * t204 + (t155 * t130 * t168 - t133 * t111 - t132 * t112 - t164 * t202) * t128) * t117, 0, t166 * t185 * t194 + (-t166 * t178 + (t166 * t188 + ((t127 * t143 + t176) * t142 + (-t111 * t142 + (t130 * t143 - t112) * t141) * t128) * t154) * t153) * t117, t104, t104, 0;];
JaD_rot  = t1;
