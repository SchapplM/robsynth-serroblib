% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobiaD_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:35
% EndTime: 2019-02-26 20:28:36
% DurationCPUTime: 0.75s
% Computational Cost: add. (813->89), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->90)
t126 = sin(qJ(6));
t130 = cos(qJ(4));
t131 = cos(qJ(1));
t169 = t130 * t131;
t128 = sin(qJ(1));
t129 = cos(qJ(6));
t172 = t128 * t129;
t109 = -t126 * t169 - t172;
t104 = 0.1e1 / t109;
t105 = 0.1e1 / t109 ^ 2;
t107 = t126 * t128 - t129 * t169;
t178 = t107 * t126;
t140 = -t104 * t129 + t105 * t178;
t103 = t107 ^ 2;
t102 = t103 * t105 + 0.1e1;
t99 = 0.1e1 / t102;
t191 = t140 * t99;
t166 = qJD(4) * t128;
t121 = t128 ^ 2;
t127 = sin(qJ(4));
t120 = t127 ^ 2;
t123 = 0.1e1 / t130 ^ 2;
t175 = t120 * t123;
t117 = t121 * t175 + 0.1e1;
t115 = 0.1e1 / t117;
t122 = 0.1e1 / t130;
t154 = t123 * t166;
t167 = qJD(1) * t131;
t156 = t127 * t167;
t165 = qJD(4) * t130;
t88 = ((-t128 * t165 - t156) * t122 - t120 * t154) * t115;
t147 = -t88 - t166;
t190 = t128 * t88 + qJD(4);
t173 = t128 * t127;
t114 = atan2(-t173, t130);
t113 = cos(t114);
t112 = sin(t114);
t158 = t112 * t173;
t97 = t113 * t130 - t158;
t94 = 0.1e1 / t97;
t95 = 0.1e1 / t97 ^ 2;
t189 = t115 - 0.1e1;
t125 = t131 ^ 2;
t155 = t128 * t167;
t159 = t95 * t165;
t177 = t113 * t127;
t83 = -t190 * t177 + (t147 * t130 - t156) * t112;
t187 = t83 * t94 * t95;
t93 = t120 * t125 * t95 + 0.1e1;
t188 = (t125 * t127 * t159 + (-t125 * t187 - t95 * t155) * t120) / t93 ^ 2;
t91 = 0.1e1 / t93;
t186 = t91 * t95;
t185 = t94 * t91;
t179 = t105 * t107;
t145 = qJD(1) * t130 + qJD(6);
t164 = qJD(4) * t131;
t139 = t127 * t164 + t145 * t128;
t146 = qJD(6) * t130 + qJD(1);
t170 = t129 * t131;
t90 = t139 * t126 - t146 * t170;
t183 = t104 * t105 * t90;
t141 = t146 * t126;
t89 = t139 * t129 + t131 * t141;
t184 = 0.1e1 / t102 ^ 2 * (-t103 * t183 + t89 * t179);
t181 = t131 * t95;
t176 = t120 * t122;
t174 = t122 * t127;
t171 = t128 * t130;
t168 = qJD(1) * t128;
t163 = 0.2e1 * t187;
t162 = -0.2e1 * t184;
t119 = t127 * t120;
t142 = t119 * t122 * t123 + t174;
t161 = 0.2e1 / t117 ^ 2 * (t142 * t121 * qJD(4) + t155 * t175);
t160 = t94 * t188;
t157 = t115 * t176;
t153 = 0.2e1 * t95 * t188;
t152 = 0.1e1 + t175;
t151 = -0.2e1 * t107 * t183;
t150 = t127 * t161;
t149 = t128 * t161;
t144 = t128 * t157;
t143 = t152 * t131;
t111 = t126 * t171 - t170;
t110 = -t126 * t131 - t129 * t171;
t101 = t152 * t128 * t115;
t87 = (t189 * t127 * t112 + t113 * t144) * t131;
t86 = -t112 * t171 - t177 - (-t112 * t130 - t113 * t173) * t101;
t84 = t152 * t149 + (-qJD(1) * t143 - 0.2e1 * t142 * t166) * t115;
t1 = [t122 * t131 * t150 + (-qJD(4) * t143 + t168 * t174) * t115, 0, 0, t84, 0, 0; (-t165 * t185 + (0.2e1 * t160 + (qJD(1) * t87 + t83) * t186) * t127) * t128 + (t87 * t153 * t127 + (-t87 * t159 + (t87 * t163 + ((t88 * t144 - t189 * t165 + t150) * t112 + (t149 * t176 + t127 * t88 + (-t119 * t154 - (t88 + 0.2e1 * t166) * t127) * t115) * t113) * t181) * t127 + (-t94 + (-(-t121 + t125) * t113 * t157 + t189 * t158) * t95) * t127 * qJD(1)) * t91) * t131, 0, 0 (-t168 * t185 + (-0.2e1 * t160 + (-qJD(4) * t86 - t83) * t186) * t131) * t130 + (t86 * t131 * t153 + (-t94 * t164 - ((t101 * t167 - t128 * t84) * t113 + (-t190 * t101 - t147) * t112) * t127 * t181 + (t131 * t163 + t95 * t168) * t86) * t91 - ((-t84 - t167) * t112 + (-t147 * t101 - t190) * t113) * t169 * t186) * t127, 0, 0; 0.2e1 * (-t104 * t110 - t111 * t179) * t184 + (t111 * t151 + (-t110 * t90 + t111 * t89 + t146 * t107 * t172 + (-t127 * t166 + t145 * t131) * t178) * t105 + (-t145 * t170 + (qJD(4) * t127 * t129 + t141) * t128) * t104) * t99, 0, 0, t130 * t164 * t191 + (-t168 * t191 + (t140 * t162 + ((qJD(6) * t104 + t151) * t126 + (t126 * t89 + (qJD(6) * t107 + t90) * t129) * t105) * t99) * t131) * t127, 0, t162 + 0.2e1 * (t105 * t89 * t99 + (-t105 * t184 - t99 * t183) * t107) * t107;];
JaD_rot  = t1;
