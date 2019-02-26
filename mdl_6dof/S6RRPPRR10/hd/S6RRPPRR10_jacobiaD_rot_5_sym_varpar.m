% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR10_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:42
% EndTime: 2019-02-26 21:33:43
% DurationCPUTime: 0.64s
% Computational Cost: add. (1161->91), mult. (2519->208), div. (480->12), fcn. (2968->9), ass. (0->92)
t134 = sin(qJ(1));
t135 = cos(qJ(2));
t133 = sin(qJ(2));
t152 = qJD(1) * t133 + qJD(5);
t136 = cos(qJ(1));
t172 = qJD(2) * t136;
t197 = t152 * t134 - t135 * t172;
t173 = qJD(2) * t135;
t196 = t134 * t173 + t152 * t136;
t180 = t134 * t135;
t118 = atan2(-t180, t133);
t117 = cos(t118);
t116 = sin(t118);
t166 = t116 * t180;
t105 = t117 * t133 - t166;
t102 = 0.1e1 / t105;
t125 = pkin(10) + qJ(5);
t123 = sin(t125);
t124 = cos(t125);
t181 = t134 * t124;
t182 = t133 * t136;
t113 = t123 * t182 + t181;
t109 = 0.1e1 / t113;
t126 = 0.1e1 / t133;
t103 = 0.1e1 / t105 ^ 2;
t110 = 0.1e1 / t113 ^ 2;
t127 = 0.1e1 / t133 ^ 2;
t129 = t134 ^ 2;
t131 = t135 ^ 2;
t185 = t127 * t131;
t122 = t129 * t185 + 0.1e1;
t120 = 0.1e1 / t122;
t195 = t120 - 0.1e1;
t112 = t123 * t134 - t124 * t182;
t108 = t112 ^ 2;
t101 = t108 * t110 + 0.1e1;
t189 = t110 * t112;
t153 = qJD(5) * t133 + qJD(1);
t147 = t153 * t136;
t94 = -t197 * t123 + t124 * t147;
t191 = t109 * t110 * t94;
t93 = t123 * t147 + t197 * t124;
t194 = (-t108 * t191 + t93 * t189) / t101 ^ 2;
t132 = t136 ^ 2;
t184 = t131 * t132;
t100 = t103 * t184 + 0.1e1;
t96 = 0.1e1 / t100;
t193 = t103 * t96;
t174 = qJD(2) * t134;
t163 = t127 * t174;
t176 = qJD(1) * t136;
t164 = t135 * t176;
t95 = ((t133 * t174 - t164) * t126 + t131 * t163) * t120;
t154 = -t95 + t174;
t155 = -t134 * t95 + qJD(2);
t187 = t117 * t135;
t89 = t155 * t187 + (t154 * t133 - t164) * t116;
t192 = t102 * t103 * t89;
t190 = t109 * t124;
t188 = t112 * t123;
t186 = t126 * t131;
t183 = t133 * t134;
t178 = qJD(1) * t134;
t177 = qJD(1) * t135;
t175 = qJD(2) * t133;
t150 = t131 * t134 * t176;
t171 = 0.2e1 * (-t184 * t192 + (-t132 * t133 * t173 - t150) * t103) / t100 ^ 2;
t170 = 0.2e1 * t194;
t169 = 0.2e1 * t192;
t130 = t135 * t131;
t145 = qJD(2) * (-t127 * t130 - t135) * t126;
t168 = 0.2e1 * (t127 * t150 + t129 * t145) / t122 ^ 2;
t167 = t96 * t175;
t165 = t120 * t186;
t160 = 0.1e1 + t185;
t159 = t102 * t171;
t158 = 0.2e1 * t112 * t191;
t157 = t134 * t168;
t156 = t135 * t168;
t151 = t134 * t165;
t149 = t160 * t136;
t148 = t153 * t134;
t146 = t110 * t188 + t190;
t98 = 0.1e1 / t101;
t144 = t146 * t98;
t115 = -t123 * t183 + t124 * t136;
t114 = t123 * t136 + t133 * t181;
t107 = t160 * t134 * t120;
t92 = (t195 * t135 * t116 + t117 * t151) * t136;
t91 = t116 * t183 + t187 + (-t116 * t133 - t117 * t180) * t107;
t90 = -t160 * t157 + (qJD(1) * t149 + 0.2e1 * t134 * t145) * t120;
t1 = [t126 * t136 * t156 + (t126 * t134 * t177 + qJD(2) * t149) * t120, t90, 0, 0, 0, 0; (t102 * t167 + (t159 + (qJD(1) * t92 + t89) * t193) * t135) * t134 + (t92 * t135 * t96 * t169 + (t92 * t167 + (t92 * t171 + ((t95 * t151 + t195 * t175 + t156) * t116 + (t157 * t186 + t135 * t95 + (t130 * t163 - (t95 - 0.2e1 * t174) * t135) * t120) * t117) * t96 * t136) * t135) * t103 + (-t102 + ((t129 - t132) * t117 * t165 + t195 * t166) * t103) * t96 * t177) * t136 (t102 * t96 * t178 + (t159 + (qJD(2) * t91 + t89) * t193) * t136) * t133 + (t91 * t136 * t103 * t171 + ((-qJD(2) * t102 + t91 * t169) * t136 + (t91 * t178 + (-(-t107 * t176 - t134 * t90) * t117 - ((t107 * t134 - 0.1e1) * t95 + (-t107 + t134) * qJD(2)) * t116) * t135 * t136) * t103) * t96 - ((-t90 + t176) * t116 + (t154 * t107 - t155) * t117) * t182 * t193) * t135, 0, 0, 0, 0; (-t109 * t114 + t115 * t189) * t170 + (t115 * t158 - t109 * t123 * t148 + t196 * t190 + (t112 * t124 * t148 - t114 * t94 - t115 * t93 + t196 * t188) * t110) * t98, t133 * t144 * t172 + (t144 * t178 + (t146 * t170 + ((qJD(5) * t109 + t158) * t123 + (-t123 * t93 + (-qJD(5) * t112 + t94) * t124) * t110) * t98) * t136) * t135, 0, 0, -0.2e1 * t194 + 0.2e1 * (t110 * t93 * t98 + (-t110 * t194 - t98 * t191) * t112) * t112, 0;];
JaD_rot  = t1;
