% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_jacobiaD_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:51
% EndTime: 2019-02-26 20:59:52
% DurationCPUTime: 0.70s
% Computational Cost: add. (822->93), mult. (2545->208), div. (484->12), fcn. (2996->9), ass. (0->94)
t137 = cos(qJ(1));
t199 = 0.2e1 * t137;
t134 = sin(qJ(1));
t136 = cos(qJ(3));
t133 = sin(qJ(3));
t152 = qJD(1) * t133 + qJD(4);
t171 = qJD(3) * t137;
t198 = t152 * t134 - t136 * t171;
t177 = t137 * t136;
t121 = atan2(t177, -t133);
t119 = sin(t121);
t120 = cos(t121);
t105 = t119 * t177 - t120 * t133;
t102 = 0.1e1 / t105;
t132 = sin(qJ(4));
t135 = cos(qJ(4));
t180 = t134 * t135;
t116 = t132 * t137 + t133 * t180;
t112 = 0.1e1 / t116;
t125 = 0.1e1 / t133;
t103 = 0.1e1 / t105 ^ 2;
t113 = 0.1e1 / t116 ^ 2;
t126 = 0.1e1 / t133 ^ 2;
t131 = t137 ^ 2;
t130 = t136 ^ 2;
t184 = t126 * t130;
t124 = t131 * t184 + 0.1e1;
t122 = 0.1e1 / t124;
t197 = t122 - 0.1e1;
t128 = t134 ^ 2;
t183 = t128 * t130;
t101 = t103 * t183 + 0.1e1;
t174 = qJD(1) * t137;
t149 = t130 * t134 * t174;
t172 = qJD(3) * t136;
t160 = t126 * t171;
t175 = qJD(1) * t136;
t162 = t134 * t175;
t96 = (-(-t133 * t171 - t162) * t125 + t130 * t160) * t122;
t154 = t96 - t171;
t155 = t137 * t96 - qJD(3);
t186 = t120 * t136;
t91 = t155 * t186 + (t154 * t133 - t162) * t119;
t194 = t102 * t103 * t91;
t196 = 0.1e1 / t101 ^ 2 * (-t183 * t194 + (-t128 * t133 * t172 + t149) * t103);
t99 = 0.1e1 / t101;
t195 = t103 * t99;
t129 = t136 * t130;
t145 = qJD(3) * (-t126 * t129 - t136) * t125;
t192 = (-t126 * t149 + t131 * t145) / t124 ^ 2;
t178 = t135 * t137;
t181 = t134 * t132;
t115 = -t133 * t181 + t178;
t111 = t115 ^ 2;
t191 = t111 * t113;
t114 = t112 * t113;
t190 = t111 * t114;
t189 = t112 * t132;
t188 = t113 * t115;
t187 = t115 * t135;
t185 = t125 * t130;
t182 = t133 * t137;
t179 = t134 * t136;
t176 = qJD(1) * t134;
t173 = qJD(3) * t133;
t170 = 0.2e1 * t196;
t169 = 0.2e1 * t194;
t110 = 0.1e1 + t191;
t153 = -qJD(4) * t133 - qJD(1);
t97 = t153 * t180 + (-t134 * t172 - t152 * t137) * t132;
t166 = t97 * t188;
t98 = t152 * t178 + (t153 * t132 + t135 * t172) * t134;
t168 = 0.2e1 / t110 ^ 2 * (-t98 * t190 + t166);
t167 = t102 * t196;
t165 = t99 * t173;
t164 = t136 * t192;
t163 = t122 * t185;
t161 = t136 * t174;
t158 = 0.1e1 + t184;
t157 = 0.2e1 * t114 * t115 * t98;
t156 = t192 * t199;
t151 = t137 * t163;
t150 = t197 * t136 * t119;
t148 = t158 * t134;
t147 = t153 * t137;
t146 = t113 * t187 + t189;
t118 = t133 * t178 - t181;
t117 = -t132 * t182 - t180;
t109 = t158 * t137 * t122;
t107 = 0.1e1 / t110;
t95 = (t120 * t151 + t150) * t134;
t94 = -t119 * t182 - t186 + (t119 * t133 + t120 * t177) * t109;
t92 = -t158 * t156 + (-qJD(1) * t148 + t145 * t199) * t122;
t1 = [-0.2e1 * t134 * t125 * t164 + (-qJD(3) * t148 + t125 * t161) * t122, 0, t92, 0, 0, 0; (-t102 * t165 + (-0.2e1 * t167 + (-qJD(1) * t95 - t91) * t195) * t136) * t137 + (t95 * t136 * t99 * t169 + (t95 * t165 + (t95 * t170 + ((t96 * t151 + t197 * t173 + 0.2e1 * t164) * t119 + (t156 * t185 + t136 * t96 + (t129 * t160 - (t96 - 0.2e1 * t171) * t136) * t122) * t120) * t99 * t134) * t136) * t103 + (-t102 + ((t128 - t131) * t120 * t163 - t137 * t150) * t103) * t99 * t175) * t134, 0 (-t102 * t99 * t174 + (0.2e1 * t167 + (qJD(3) * t94 + t91) * t195) * t134) * t133 + (((-qJD(3) * t102 + t94 * t169) * t134 + (-t94 * t174 + (-(-t109 * t176 + t137 * t92) * t120 - ((-t109 * t137 + 0.1e1) * t96 + (t109 - t137) * qJD(3)) * t119) * t179) * t103) * t99 + (t94 * t170 - ((t92 + t176) * t119 + (t154 * t109 - t155) * t120) * t99 * t133) * t103 * t134) * t136, 0, 0, 0; (-t112 * t117 + t118 * t188) * t168 + (t118 * t157 + t112 * t135 * t147 + t198 * t189 + (-t115 * t132 * t147 - t117 * t98 - t118 * t97 + t187 * t198) * t113) * t107, 0, t146 * t168 * t179 + (-t146 * t161 + (t146 * t173 + ((-qJD(4) * t112 + t157) * t135 + (-t135 * t97 + (qJD(4) * t115 + t98) * t132) * t113) * t136) * t134) * t107 (t112 * t116 + t191) * t168 + (-0.2e1 * t166 + (t113 * t116 - t112 + 0.2e1 * t190) * t98) * t107, 0, 0;];
JaD_rot  = t1;
