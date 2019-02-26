% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP4_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_jacobiaD_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:36:31
% EndTime: 2019-02-26 21:36:31
% DurationCPUTime: 0.64s
% Computational Cost: add. (813->90), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->93)
t131 = cos(qJ(2));
t129 = sin(qJ(1));
t123 = t129 ^ 2;
t128 = sin(qJ(2));
t121 = 0.1e1 / t128 ^ 2;
t125 = t131 ^ 2;
t179 = t121 * t125;
t118 = t123 * t179 + 0.1e1;
t116 = 0.1e1 / t118;
t120 = 0.1e1 / t128;
t169 = qJD(2) * t129;
t157 = t121 * t169;
t132 = cos(qJ(1));
t171 = qJD(1) * t132;
t158 = t131 * t171;
t90 = ((t128 * t169 - t158) * t120 + t125 * t157) * t116;
t193 = t131 * t90;
t148 = -t90 + t169;
t146 = qJD(1) * t128 + qJD(4);
t167 = qJD(2) * t132;
t192 = t129 * t146 - t131 * t167;
t176 = t129 * t131;
t115 = atan2(-t176, t128);
t114 = cos(t115);
t113 = sin(t115);
t161 = t113 * t176;
t99 = t114 * t128 - t161;
t96 = 0.1e1 / t99;
t127 = sin(qJ(4));
t175 = t132 * t127;
t130 = cos(qJ(4));
t177 = t129 * t130;
t110 = t128 * t175 + t177;
t106 = 0.1e1 / t110;
t107 = 0.1e1 / t110 ^ 2;
t97 = 0.1e1 / t99 ^ 2;
t191 = t116 - 0.1e1;
t149 = -t129 * t90 + qJD(2);
t181 = t114 * t131;
t85 = t149 * t181 + (t128 * t148 - t158) * t113;
t190 = t85 * t96 * t97;
t126 = t132 ^ 2;
t95 = t126 * t125 * t97 + 0.1e1;
t93 = 0.1e1 / t95;
t189 = t93 * t97;
t188 = t96 * t93;
t174 = t132 * t130;
t178 = t129 * t127;
t109 = -t128 * t174 + t178;
t105 = t109 ^ 2;
t104 = t105 * t107 + 0.1e1;
t184 = t107 * t109;
t147 = qJD(4) * t128 + qJD(1);
t143 = t147 * t132;
t92 = -t192 * t127 + t130 * t143;
t186 = t106 * t107 * t92;
t91 = t127 * t143 + t192 * t130;
t187 = 0.1e1 / t104 ^ 2 * (-t105 * t186 + t184 * t91);
t185 = t132 * t97;
t183 = t109 * t127;
t182 = t113 * t129;
t180 = t120 * t125;
t173 = qJD(1) * t129;
t172 = qJD(1) * t131;
t170 = qJD(2) * t128;
t168 = qJD(2) * t131;
t159 = t129 * t171;
t162 = t97 * t170;
t166 = 0.2e1 * (-t126 * t131 * t162 + (-t126 * t190 - t159 * t97) * t125) / t95 ^ 2;
t165 = 0.2e1 * t190;
t164 = 0.2e1 * t187;
t124 = t131 * t125;
t141 = qJD(2) * (-t121 * t124 - t131) * t120;
t163 = 0.2e1 * (t123 * t141 + t159 * t179) / t118 ^ 2;
t160 = t116 * t180;
t155 = t96 * t166;
t154 = t97 * t166;
t153 = 0.1e1 + t179;
t152 = 0.2e1 * t109 * t186;
t151 = t129 * t163;
t150 = t131 * t163;
t145 = t129 * t160;
t144 = t153 * t132;
t142 = t106 * t130 + t107 * t183;
t140 = t142 * t132;
t112 = -t128 * t178 + t174;
t111 = t128 * t177 + t175;
t103 = t153 * t129 * t116;
t101 = 0.1e1 / t104;
t89 = (t113 * t131 * t191 + t114 * t145) * t132;
t88 = t128 * t182 + t181 + (-t113 * t128 - t114 * t176) * t103;
t86 = -t153 * t151 + (qJD(1) * t144 + 0.2e1 * t129 * t141) * t116;
t1 = [t132 * t120 * t150 + (t120 * t129 * t172 + qJD(2) * t144) * t116, t86, 0, 0, 0, 0; (t170 * t188 + (t155 + (qJD(1) * t89 + t85) * t189) * t131) * t129 + (t89 * t154 * t131 + (t89 * t162 + (t89 * t165 + ((t145 * t90 + t170 * t191 + t150) * t113 + (t151 * t180 + t193 + (t124 * t157 - (t90 - 0.2e1 * t169) * t131) * t116) * t114) * t185) * t131 + (-t96 + (-(-t123 + t126) * t114 * t160 + t191 * t161) * t97) * t172) * t93) * t132 (t173 * t188 + (t155 + (qJD(2) * t88 + t85) * t189) * t132) * t128 + (t88 * t132 * t154 + (-t96 * t167 + (t132 * t165 + t173 * t97) * t88 + (-t103 * t182 * t193 + (-(-t103 * t171 - t129 * t86) * t131 - (t103 * t148 - t149) * t128) * t114 + (-(-qJD(2) * t103 + t148) * t131 - (-t86 + t171) * t128) * t113) * t185) * t93) * t131, 0, 0, 0, 0; (-t106 * t111 + t112 * t184) * t164 + (t112 * t152 + (-t111 * t92 - t112 * t91 + t147 * t109 * t177 - (-t129 * t168 - t132 * t146) * t183) * t107 + (t146 * t174 + (-t127 * t147 + t130 * t168) * t129) * t106) * t101, t131 * t140 * t164 + (t140 * t170 + (t142 * t173 + ((qJD(4) * t106 + t152) * t127 + (-t127 * t91 + (-qJD(4) * t109 + t92) * t130) * t107) * t132) * t131) * t101, 0, -0.2e1 * t187 + 0.2e1 * (t91 * t107 * t101 + (-t101 * t186 - t107 * t187) * t109) * t109, 0, 0;];
JaD_rot  = t1;
