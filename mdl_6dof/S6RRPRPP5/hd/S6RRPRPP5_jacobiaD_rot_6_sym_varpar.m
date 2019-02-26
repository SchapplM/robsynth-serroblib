% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobiaD_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:07
% EndTime: 2019-02-26 21:37:08
% DurationCPUTime: 0.67s
% Computational Cost: add. (822->93), mult. (2545->206), div. (484->12), fcn. (2996->9), ass. (0->93)
t131 = sin(qJ(2));
t132 = sin(qJ(1));
t134 = cos(qJ(2));
t178 = t132 * t134;
t119 = atan2(t178, -t131);
t118 = cos(t119);
t117 = sin(t119);
t163 = t117 * t178;
t103 = -t118 * t131 + t163;
t100 = 0.1e1 / t103;
t130 = sin(qJ(4));
t133 = cos(qJ(4));
t179 = t132 * t133;
t135 = cos(qJ(1));
t181 = t131 * t135;
t114 = t130 * t181 + t179;
t110 = 0.1e1 / t114;
t123 = 0.1e1 / t131;
t101 = 0.1e1 / t103 ^ 2;
t111 = 0.1e1 / t114 ^ 2;
t124 = 0.1e1 / t131 ^ 2;
t126 = t132 ^ 2;
t128 = t134 ^ 2;
t184 = t124 * t128;
t122 = t126 * t184 + 0.1e1;
t120 = 0.1e1 / t122;
t195 = t120 - 0.1e1;
t129 = t135 ^ 2;
t183 = t128 * t129;
t99 = t101 * t183 + 0.1e1;
t97 = 0.1e1 / t99;
t194 = t101 * t97;
t171 = qJD(2) * t132;
t160 = t124 * t171;
t173 = qJD(1) * t135;
t161 = t134 * t173;
t94 = (-(-t131 * t171 + t161) * t123 + t128 * t160) * t120;
t153 = t94 - t171;
t154 = t132 * t94 - qJD(2);
t186 = t118 * t134;
t89 = t154 * t186 + (t153 * t131 + t161) * t117;
t193 = t100 * t101 * t89;
t177 = t133 * t135;
t180 = t132 * t130;
t113 = t131 * t177 - t180;
t109 = t113 ^ 2;
t192 = t109 * t111;
t112 = t110 * t111;
t191 = t109 * t112;
t190 = t110 * t133;
t189 = t111 * t113;
t188 = t113 * t130;
t187 = t117 * t131;
t185 = t123 * t128;
t182 = t130 * t135;
t175 = qJD(1) * t132;
t174 = qJD(1) * t134;
t172 = qJD(2) * t131;
t170 = qJD(2) * t134;
t149 = t128 * t132 * t173;
t169 = -0.2e1 * (-t183 * t193 + (-t129 * t131 * t170 - t149) * t101) / t99 ^ 2;
t168 = -0.2e1 * t193;
t108 = 0.1e1 + t192;
t151 = -qJD(1) * t131 - qJD(4);
t143 = t151 * t132 + t135 * t170;
t152 = qJD(4) * t131 + qJD(1);
t95 = t143 * t133 - t152 * t182;
t165 = t95 * t189;
t96 = t143 * t130 + t152 * t177;
t167 = 0.2e1 / t108 ^ 2 * (-t96 * t191 + t165);
t127 = t134 * t128;
t145 = qJD(2) * (-t124 * t127 - t134) * t123;
t166 = 0.2e1 * (t124 * t149 + t126 * t145) / t122 ^ 2;
t164 = t97 * t172;
t162 = t120 * t185;
t159 = 0.1e1 + t184;
t158 = t100 * t169;
t157 = 0.2e1 * t112 * t113 * t96;
t156 = t132 * t166;
t155 = t134 * t166;
t150 = t132 * t162;
t148 = t159 * t135;
t147 = t151 * t135;
t146 = t111 * t188 - t190;
t144 = t146 * t135;
t116 = -t131 * t180 + t177;
t115 = -t131 * t179 - t182;
t107 = t159 * t132 * t120;
t105 = 0.1e1 / t108;
t93 = (-t195 * t134 * t117 - t118 * t150) * t135;
t92 = -t132 * t187 - t186 + (t118 * t178 + t187) * t107;
t90 = -t159 * t156 + (qJD(1) * t148 + 0.2e1 * t132 * t145) * t120;
t1 = [t123 * t135 * t155 + (t123 * t132 * t174 + qJD(2) * t148) * t120, t90, 0, 0, 0, 0; (-t100 * t164 + (t158 + (-qJD(1) * t93 - t89) * t194) * t134) * t132 + (t93 * t134 * t97 * t168 + (-t93 * t164 + (t93 * t169 + ((t94 * t150 + t195 * t172 + t155) * t117 + (t156 * t185 + t134 * t94 + (t127 * t160 + (-t94 + 0.2e1 * t171) * t134) * t120) * t118) * t97 * t135) * t134) * t101 + (t100 + ((t126 - t129) * t118 * t162 + t195 * t163) * t101) * t97 * t174) * t135 (-t100 * t97 * t175 + (t158 + (-qJD(2) * t92 - t89) * t194) * t135) * t131 + (t92 * t135 * t101 * t169 + ((qJD(2) * t100 + t92 * t168) * t135 + (-t92 * t175 + ((t107 * t173 + t132 * t90) * t118 + ((-t107 * t132 + 0.1e1) * t94 + (t107 - t132) * qJD(2)) * t117) * t134 * t135) * t101) * t97 + ((t90 - t173) * t117 + (t153 * t107 - t154) * t118) * t181 * t194) * t134, 0, 0, 0, 0; (-t110 * t115 + t116 * t189) * t167 + (t116 * t157 + t147 * t190 + (t152 * t130 - t133 * t170) * t110 * t132 + (-t115 * t96 - t116 * t95 + t152 * t113 * t179 - (-t132 * t170 + t147) * t188) * t111) * t105, t134 * t144 * t167 + (t144 * t172 + (t146 * t175 + ((-qJD(4) * t110 + t157) * t130 + (-t130 * t95 + (-qJD(4) * t113 - t96) * t133) * t111) * t135) * t134) * t105, 0 (t110 * t114 + t192) * t167 + (-0.2e1 * t165 + (t111 * t114 - t110 + 0.2e1 * t191) * t96) * t105, 0, 0;];
JaD_rot  = t1;
