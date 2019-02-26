% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:22:26
% EndTime: 2019-02-26 21:22:27
% DurationCPUTime: 0.70s
% Computational Cost: add. (2429->94), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
t146 = qJ(2) + pkin(9);
t142 = sin(t146);
t137 = 0.1e1 / t142 ^ 2;
t144 = cos(t146);
t140 = t144 ^ 2;
t194 = t137 * t140;
t213 = t144 * t194;
t149 = sin(qJ(1));
t170 = 0.1e1 + t194;
t212 = t149 * t170;
t166 = qJD(1) * t142 + qJD(6);
t150 = cos(qJ(1));
t182 = qJD(2) * t150;
t211 = -t144 * t182 + t166 * t149;
t183 = qJD(2) * t149;
t210 = t144 * t183 + t166 * t150;
t187 = t149 * t144;
t131 = atan2(-t187, t142);
t130 = cos(t131);
t129 = sin(t131);
t175 = t129 * t187;
t115 = t130 * t142 - t175;
t112 = 0.1e1 / t115;
t145 = pkin(10) + qJ(6);
t141 = sin(t145);
t143 = cos(t145);
t188 = t149 * t143;
t190 = t142 * t150;
t126 = t141 * t190 + t188;
t122 = 0.1e1 / t126;
t136 = 0.1e1 / t142;
t113 = 0.1e1 / t115 ^ 2;
t123 = 0.1e1 / t126 ^ 2;
t147 = t149 ^ 2;
t134 = t147 * t194 + 0.1e1;
t132 = 0.1e1 / t134;
t209 = t132 - 0.1e1;
t185 = qJD(1) * t150;
t173 = t144 * t185;
t106 = ((t142 * t183 - t173) * t136 + t183 * t194) * t132;
t196 = t130 * t144;
t101 = (-t106 * t149 + qJD(2)) * t196 + (-t173 + (-t106 + t183) * t142) * t129;
t208 = t101 * t112 * t113;
t167 = qJD(6) * t142 + qJD(1);
t161 = t167 * t150;
t107 = t141 * t161 + t143 * t211;
t189 = t143 * t150;
t125 = t141 * t149 - t142 * t189;
t121 = t125 ^ 2;
t120 = t121 * t123 + 0.1e1;
t198 = t123 * t125;
t108 = -t141 * t211 + t143 * t161;
t204 = t108 * t122 * t123;
t207 = (t107 * t198 - t121 * t204) / t120 ^ 2;
t206 = t106 * t129;
t205 = t106 * t144;
t203 = t113 * t144;
t202 = t113 * t150;
t195 = t136 * t144;
t159 = qJD(2) * (-t136 * t213 - t195);
t192 = t140 * t149;
t164 = t185 * t192;
t201 = (t137 * t164 + t147 * t159) / t134 ^ 2;
t119 = t132 * t212;
t200 = t119 * t149;
t199 = t122 * t143;
t197 = t125 * t141;
t148 = t150 ^ 2;
t193 = t140 * t148;
t191 = t142 * t149;
t186 = qJD(1) * t149;
t184 = qJD(2) * t142;
t111 = t113 * t193 + 0.1e1;
t181 = 0.2e1 * (-t193 * t208 + (-t144 * t148 * t184 - t164) * t113) / t111 ^ 2;
t180 = 0.2e1 * t208;
t179 = 0.2e1 * t207;
t178 = -0.2e1 * t201;
t177 = t144 * t202;
t176 = t144 * t201;
t174 = t136 * t192;
t169 = t144 * t181;
t168 = 0.2e1 * t125 * t204;
t165 = t132 * t174;
t163 = t170 * t150;
t162 = t167 * t149;
t160 = t123 * t197 + t199;
t158 = t160 * t150;
t128 = -t141 * t191 + t189;
t127 = t141 * t150 + t142 * t188;
t117 = 0.1e1 / t120;
t109 = 0.1e1 / t111;
t105 = (t209 * t144 * t129 + t130 * t165) * t150;
t104 = t129 * t191 + t196 + (-t129 * t142 - t130 * t187) * t119;
t102 = t178 * t212 + (qJD(1) * t163 + 0.2e1 * t149 * t159) * t132;
t1 = [0.2e1 * t136 * t150 * t176 + (qJD(2) * t163 + t186 * t195) * t132, t102, 0, 0, 0, 0; (t112 * t169 + (t112 * t184 + (qJD(1) * t105 + t101) * t203) * t109) * t149 + (t113 * t169 * t105 + (-((-t106 * t165 - t209 * t184 - 0.2e1 * t176) * t129 + (t174 * t178 - t205 + (t205 + (-0.2e1 * t144 - t213) * t183) * t132) * t130) * t177 + (t113 * t184 + t144 * t180) * t105 + (-t112 + ((t147 - t148) * t140 * t136 * t132 * t130 + t209 * t175) * t113) * t144 * qJD(1)) * t109) * t150 (t104 * t203 + t112 * t142) * t150 * t181 + ((t112 * t186 + (qJD(2) * t104 + t101) * t202) * t142 + (-t112 * t182 - (-t102 * t130 * t149 + t129 * t183 + t200 * t206 - t206 + (-qJD(2) * t129 - t130 * t185) * t119) * t177 + (t113 * t186 + t150 * t180) * t104 - ((-t102 + t185) * t129 + ((-0.1e1 + t200) * qJD(2) + (-t119 + t149) * t106) * t130) * t113 * t190) * t144) * t109, 0, 0, 0, 0; (-t122 * t127 + t128 * t198) * t179 + (t128 * t168 - t122 * t141 * t162 + t210 * t199 + (t125 * t143 * t162 - t128 * t107 - t127 * t108 + t197 * t210) * t123) * t117, t144 * t158 * t179 + (t158 * t184 + (t160 * t186 + ((qJD(6) * t122 + t168) * t141 + (-t107 * t141 + (-qJD(6) * t125 + t108) * t143) * t123) * t150) * t144) * t117, 0, 0, 0, -0.2e1 * t207 + 0.2e1 * (t107 * t117 * t123 + (-t117 * t204 - t123 * t207) * t125) * t125;];
JaD_rot  = t1;
