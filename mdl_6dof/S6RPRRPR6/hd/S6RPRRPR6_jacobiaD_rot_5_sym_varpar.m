% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:47
% EndTime: 2019-02-26 21:03:48
% DurationCPUTime: 0.74s
% Computational Cost: add. (2618->95), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->94)
t145 = pkin(10) + qJ(3);
t141 = sin(t145);
t137 = t141 ^ 2;
t143 = cos(t145);
t139 = 0.1e1 / t143 ^ 2;
t195 = t137 * t139;
t149 = sin(qJ(1));
t212 = 0.2e1 * t149;
t211 = t141 * t195;
t147 = t149 ^ 2;
t132 = t147 * t195 + 0.1e1;
t130 = 0.1e1 / t132;
t138 = 0.1e1 / t143;
t150 = cos(qJ(1));
t184 = qJD(1) * t150;
t172 = t141 * t184;
t182 = qJD(3) * t149;
t104 = (-(-t143 * t182 - t172) * t138 + t182 * t195) * t130;
t210 = t104 - t182;
t146 = qJ(4) + pkin(11);
t142 = sin(t146);
t187 = t149 * t142;
t144 = cos(t146);
t189 = t144 * t150;
t126 = t143 * t189 + t187;
t188 = t149 * t141;
t129 = atan2(-t188, -t143);
t128 = cos(t129);
t127 = sin(t129);
t175 = t127 * t188;
t113 = -t128 * t143 - t175;
t110 = 0.1e1 / t113;
t120 = 0.1e1 / t126;
t111 = 0.1e1 / t113 ^ 2;
t121 = 0.1e1 / t126 ^ 2;
t209 = -0.2e1 * t141;
t208 = t130 - 0.1e1;
t197 = t128 * t141;
t99 = (-t104 * t149 + qJD(3)) * t197 + (t210 * t143 - t172) * t127;
t207 = t110 * t111 * t99;
t160 = t143 * t187 + t189;
t181 = qJD(3) * t150;
t171 = t141 * t181;
t105 = t160 * qJD(1) - t126 * qJD(4) + t142 * t171;
t186 = t149 * t144;
t190 = t143 * t150;
t125 = t142 * t190 - t186;
t119 = t125 ^ 2;
t118 = t119 * t121 + 0.1e1;
t200 = t121 * t125;
t165 = -qJD(1) * t143 + qJD(4);
t166 = qJD(4) * t143 - qJD(1);
t191 = t142 * t150;
t106 = -t166 * t191 + (t165 * t149 - t171) * t144;
t204 = t106 * t120 * t121;
t206 = (-t105 * t200 - t119 * t204) / t118 ^ 2;
t205 = t104 * t141;
t203 = t111 * t141;
t193 = t138 * t141;
t159 = qJD(3) * (t138 * t211 + t193);
t163 = t137 * t149 * t184;
t202 = (t139 * t163 + t147 * t159) / t132 ^ 2;
t201 = t120 * t142;
t199 = t125 * t144;
t198 = t127 * t149;
t196 = t137 * t138;
t148 = t150 ^ 2;
t194 = t137 * t148;
t192 = t141 * t150;
t185 = qJD(1) * t149;
t183 = qJD(3) * t143;
t109 = t111 * t194 + 0.1e1;
t180 = 0.2e1 / t109 ^ 2 * (-t194 * t207 + (t141 * t148 * t183 - t163) * t111);
t179 = 0.2e1 * t207;
t178 = -0.2e1 * t206;
t177 = t125 * t204;
t176 = t111 * t192;
t174 = t130 * t196;
t170 = 0.1e1 + t195;
t169 = t141 * t180;
t168 = t202 * t209;
t167 = t202 * t212;
t164 = t149 * t174;
t162 = t170 * t150;
t161 = t121 * t199 - t201;
t158 = t141 * t182 + t165 * t150;
t124 = -t143 * t186 + t191;
t117 = t170 * t149 * t130;
t115 = 0.1e1 / t118;
t107 = 0.1e1 / t109;
t103 = (t208 * t141 * t127 - t128 * t164) * t150;
t102 = -t143 * t198 + t197 + (t127 * t143 - t128 * t188) * t117;
t100 = -t170 * t167 + (qJD(1) * t162 + t159 * t212) * t130;
t1 = [t138 * t150 * t168 + (qJD(3) * t162 - t185 * t193) * t130, 0, t100, 0, 0, 0; (t110 * t169 + (-t110 * t183 + (qJD(1) * t103 + t99) * t203) * t107) * t149 + (t111 * t169 * t103 + (-((t104 * t164 + t208 * t183 + t168) * t127 + (t167 * t196 - t205 + (t205 + (t209 - t211) * t182) * t130) * t128) * t176 + (-t111 * t183 + t141 * t179) * t103 + (-t110 + ((-t147 + t148) * t128 * t174 + t208 * t175) * t111) * t141 * qJD(1)) * t107) * t150, 0 (t102 * t203 - t110 * t143) * t150 * t180 + ((-t110 * t185 + (-qJD(3) * t102 - t99) * t150 * t111) * t143 + (-t110 * t181 - (-t100 * t128 * t149 - t210 * t127 + (-qJD(3) * t127 + t104 * t198 - t128 * t184) * t117) * t176 + (t111 * t185 + t150 * t179) * t102 - ((t100 - t184) * t127 + ((-t117 * t149 + 0.1e1) * qJD(3) + (t117 - t149) * t104) * t128) * t111 * t190) * t141) * t107, 0, 0, 0; 0.2e1 * (t120 * t160 + t124 * t200) * t206 + (0.2e1 * t124 * t177 - t166 * t120 * t186 + t158 * t201 + (-t166 * t125 * t187 + t124 * t105 + t106 * t160 - t158 * t199) * t121) * t115, 0, t161 * t178 * t192 + (t161 * t143 * t181 + (-t161 * t185 + ((-qJD(4) * t120 - 0.2e1 * t177) * t144 + (-t105 * t144 + (-qJD(4) * t125 + t106) * t142) * t121) * t150) * t141) * t115, t178 + 0.2e1 * (-t105 * t115 * t121 + (-t115 * t204 - t121 * t206) * t125) * t125, 0, 0;];
JaD_rot  = t1;
