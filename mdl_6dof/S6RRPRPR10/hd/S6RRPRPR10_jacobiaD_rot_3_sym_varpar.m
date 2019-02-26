% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR10_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiaD_rot_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:00
% EndTime: 2019-02-26 21:43:01
% DurationCPUTime: 0.57s
% Computational Cost: add. (1295->85), mult. (3755->195), div. (644->14), fcn. (4884->11), ass. (0->88)
t141 = sin(qJ(2));
t143 = cos(qJ(2));
t142 = sin(qJ(1));
t181 = cos(pkin(6));
t159 = t142 * t181;
t189 = cos(qJ(1));
t151 = -t189 * t141 - t143 * t159;
t156 = t181 * t189;
t155 = t143 * t156;
t173 = t142 * t141;
t123 = -t155 + t173;
t139 = sin(pkin(6));
t174 = t139 * t143;
t117 = atan2(-t123, -t174);
t115 = sin(t117);
t116 = cos(t117);
t97 = -t115 * t123 - t116 * t174;
t95 = 0.1e1 / t97 ^ 2;
t182 = t151 * t95;
t121 = t123 ^ 2;
t134 = 0.1e1 / t139 ^ 2;
t136 = 0.1e1 / t143 ^ 2;
t120 = t121 * t134 * t136 + 0.1e1;
t118 = 0.1e1 / t120;
t133 = 0.1e1 / t139;
t135 = 0.1e1 / t143;
t165 = t123 * t133 * t135;
t190 = (t116 * t165 - t115) * t118 + t115;
t94 = 0.1e1 / t97;
t157 = t141 * t159;
t163 = t189 * t143;
t127 = t163 - t157;
t138 = sin(pkin(11));
t140 = cos(pkin(11));
t175 = t139 * t142;
t114 = t127 * t140 + t138 * t175;
t108 = 0.1e1 / t114;
t109 = 0.1e1 / t114 ^ 2;
t103 = -qJD(1) * t155 - qJD(2) * t163 + (qJD(2) * t181 + qJD(1)) * t173;
t122 = t151 ^ 2;
t150 = -t141 * t156 - t142 * t143;
t105 = -t151 * qJD(1) - t150 * qJD(2);
t172 = qJD(2) * t141;
t161 = t136 * t172;
t152 = t105 * t135 + t123 * t161;
t178 = t118 * t133;
t88 = t152 * t178;
t84 = (-t123 * t88 + t139 * t172) * t116 + (t88 * t174 - t105) * t115;
t96 = t94 * t95;
t187 = t84 * t96;
t92 = t122 * t95 + 0.1e1;
t188 = (t103 * t182 - t122 * t187) / t92 ^ 2;
t113 = t127 * t138 - t140 * t175;
t107 = t113 ^ 2;
t100 = t107 * t109 + 0.1e1;
t104 = t150 * qJD(1) + t151 * qJD(2);
t160 = t189 * qJD(1);
t158 = t139 * t160;
t101 = t104 * t138 - t140 * t158;
t179 = t109 * t113;
t102 = t104 * t140 + t138 * t158;
t180 = t102 * t108 * t109;
t186 = (t101 * t179 - t107 * t180) / t100 ^ 2;
t185 = t103 * t95;
t137 = t135 * t136;
t184 = 0.1e1 / t120 ^ 2 * (t105 * t123 * t136 + t121 * t137 * t172) * t134;
t177 = t136 * t141;
t153 = t123 * t177 - t135 * t150;
t89 = t153 * t178;
t183 = t123 * t89;
t176 = t138 * t108;
t171 = 0.2e1 * t188;
t170 = 0.2e1 * t186;
t169 = -0.2e1 * t184;
t168 = t115 * t182;
t167 = t116 * t182;
t166 = t113 * t180;
t164 = t139 * t189;
t162 = qJD(1) * t175;
t112 = t138 * t164 + t140 * t150;
t111 = t138 * t150 - t140 * t164;
t106 = -qJD(1) * t157 - t142 * t172 + (qJD(2) * t156 + t160) * t143;
t98 = 0.1e1 / t100;
t90 = 0.1e1 / t92;
t87 = t190 * t151;
t85 = (t139 * t141 - t183) * t116 + (t89 * t174 + t150) * t115;
t83 = (t153 * t169 + (t105 * t177 + t106 * t135 + (-t150 * t177 + (0.2e1 * t137 * t141 ^ 2 + t135) * t123) * qJD(2)) * t118) * t133;
t1 = [(-t151 * t135 * t169 + (-t103 * t135 - t151 * t161) * t118) * t133, t83, 0, 0, 0, 0; -0.2e1 * t188 * t87 * t182 + t123 * t94 * t171 + (-t105 * t94 + (t103 * t87 + t123 * t84) * t95 - (0.2e1 * t187 * t87 + (t118 * t88 * t165 + t169) * t168 + (0.2e1 * t165 * t184 - t88 + (-t152 * t133 + t88) * t118) * t167 - t190 * t185) * t151) * t90 (-t127 * t94 - t85 * t182) * t171 + (t85 * t185 + t104 * t94 + (-0.2e1 * t151 * t85 * t96 - t127 * t95) * t84 + (-t105 * t89 - t123 * t83 + t150 * t88 + (t88 * t89 + qJD(2)) * t174) * t167 + (t88 * t183 - t106 + (t143 * t83 + (-qJD(2) * t89 - t88) * t141) * t139) * t168) * t90, 0, 0, 0, 0; (-t108 * t111 + t112 * t179) * t170 + ((-t106 * t138 + t140 * t162) * t108 + 0.2e1 * t112 * t166 + (-t111 * t102 - (-t106 * t140 - t138 * t162) * t113 - t112 * t101) * t109) * t98 (-t140 * t179 + t176) * t98 * t103 - (-0.2e1 * t140 * t98 * t166 + t170 * t176 + (t102 * t138 * t98 + (t101 * t98 - 0.2e1 * t113 * t186) * t140) * t109) * t151, 0, 0, 0, 0;];
JaD_rot  = t1;
