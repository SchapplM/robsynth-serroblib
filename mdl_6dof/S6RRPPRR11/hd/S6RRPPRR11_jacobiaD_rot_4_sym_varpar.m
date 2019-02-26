% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR11_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:19
% EndTime: 2019-02-26 21:34:20
% DurationCPUTime: 0.58s
% Computational Cost: add. (1150->84), mult. (3755->189), div. (644->14), fcn. (4884->11), ass. (0->87)
t143 = sin(qJ(2));
t144 = cos(qJ(2));
t145 = cos(qJ(1));
t181 = cos(pkin(6));
t160 = t145 * t181;
t188 = sin(qJ(1));
t125 = t143 * t160 + t188 * t144;
t156 = t181 * t188;
t127 = t145 * t143 + t144 * t156;
t105 = t125 * qJD(1) + t127 * qJD(2);
t155 = t143 * t156;
t173 = t144 * t145;
t128 = -t155 + t173;
t123 = t128 ^ 2;
t141 = sin(pkin(6));
t176 = t141 * t143;
t120 = atan2(-t125, t176);
t116 = sin(t120);
t117 = cos(t120);
t98 = -t116 * t125 + t117 * t176;
t96 = 0.1e1 / t98 ^ 2;
t182 = t128 * t96;
t164 = t188 * t143;
t107 = -qJD(1) * t155 - qJD(2) * t164 + (qJD(2) * t181 + qJD(1)) * t173;
t172 = qJD(2) * t144;
t137 = 0.1e1 / t143;
t138 = 0.1e1 / t143 ^ 2;
t162 = t138 * t172;
t152 = -t107 * t137 + t125 * t162;
t122 = t125 ^ 2;
t136 = 0.1e1 / t141 ^ 2;
t121 = t122 * t136 * t138 + 0.1e1;
t118 = 0.1e1 / t121;
t135 = 0.1e1 / t141;
t178 = t118 * t135;
t89 = t152 * t178;
t85 = (-t125 * t89 + t141 * t172) * t117 + (-t89 * t176 - t107) * t116;
t95 = 0.1e1 / t98;
t97 = t95 * t96;
t186 = t85 * t97;
t93 = t123 * t96 + 0.1e1;
t171 = 0.2e1 * (-t105 * t182 - t123 * t186) / t93 ^ 2;
t166 = t125 * t135 * t137;
t190 = (t117 * t166 + t116) * t118 - t116;
t140 = sin(pkin(11));
t142 = cos(pkin(11));
t165 = t141 * t188;
t113 = t127 * t140 + t142 * t165;
t109 = 0.1e1 / t113;
t110 = 0.1e1 / t113 ^ 2;
t189 = 0.2e1 * t128;
t185 = t105 * t96;
t139 = t137 * t138;
t184 = 0.1e1 / t121 ^ 2 * (t107 * t125 * t138 - t122 * t139 * t172) * t136;
t157 = t144 * t160;
t124 = t164 - t157;
t177 = t138 * t144;
t153 = t124 * t137 + t125 * t177;
t90 = t153 * t178;
t183 = t125 * t90;
t161 = t188 * qJD(1);
t104 = -qJD(1) * t157 - t145 * t172 + (qJD(2) * t156 + t161) * t143;
t175 = t141 * t145;
t163 = qJD(1) * t175;
t103 = -t104 * t140 + t142 * t163;
t180 = t103 * t109 * t110;
t112 = -t127 * t142 + t140 * t165;
t179 = t110 * t112;
t174 = t142 * t109;
t108 = t112 ^ 2;
t101 = t108 * t110 + 0.1e1;
t102 = t104 * t142 + t140 * t163;
t170 = 0.2e1 / t101 ^ 2 * (t102 * t179 - t108 * t180);
t169 = -0.2e1 * t184;
t168 = t116 * t182;
t167 = t117 * t182;
t159 = 0.2e1 * t112 * t180;
t158 = t141 * t161;
t115 = -t124 * t140 + t142 * t175;
t114 = t124 * t142 + t140 * t175;
t106 = t127 * qJD(1) + qJD(2) * t125;
t99 = 0.1e1 / t101;
t91 = 0.1e1 / t93;
t88 = t190 * t128;
t86 = (t141 * t144 - t183) * t117 + (-t90 * t176 + t124) * t116;
t84 = (t153 * t169 + (t107 * t177 + t106 * t137 + (-t124 * t177 + (-0.2e1 * t139 * t144 ^ 2 - t137) * t125) * qJD(2)) * t118) * t135;
t1 = [(t137 * t184 * t189 + (t105 * t137 + t128 * t162) * t118) * t135, t84, 0, 0, 0, 0; (-t107 * t95 + (t105 * t88 + t125 * t85) * t96 + (0.2e1 * t186 * t88 - (-t118 * t89 * t166 + t169) * t168 - (t166 * t169 - t89 + (-t152 * t135 + t89) * t118) * t167 + t190 * t185) * t128) * t91 + (t125 * t95 + t88 * t182) * t171 (t127 * t95 + t86 * t182) * t171 + (t86 * t185 + t104 * t95 + (t86 * t97 * t189 + t127 * t96) * t85 - (-t107 * t90 + t124 * t89 - t125 * t84 + (-t89 * t90 - qJD(2)) * t176) * t167 - (t89 * t183 + t106 + (-t143 * t84 + (-qJD(2) * t90 - t89) * t144) * t141) * t168) * t91, 0, 0, 0, 0; (-t109 * t114 + t115 * t179) * t170 + ((t106 * t142 - t140 * t158) * t109 + t115 * t159 + (-t114 * t103 - (-t106 * t140 - t142 * t158) * t112 - t115 * t102) * t110) * t99 (t140 * t179 + t174) * t99 * t105 + (t170 * t174 + t140 * t99 * t159 + (t103 * t142 * t99 + (-t102 * t99 + t112 * t170) * t140) * t110) * t128, 0, 0, 0, 0;];
JaD_rot  = t1;
