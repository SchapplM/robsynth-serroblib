% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:56
% EndTime: 2019-02-26 22:03:57
% DurationCPUTime: 0.59s
% Computational Cost: add. (4795->70), mult. (2858->159), div. (686->14), fcn. (3330->7), ass. (0->76)
t128 = sin(qJ(1));
t122 = t128 ^ 2;
t120 = qJ(2) + qJ(3) + pkin(10);
t118 = sin(t120);
t113 = t118 ^ 2;
t119 = cos(t120);
t116 = 0.1e1 / t119 ^ 2;
t166 = t113 * t116;
t108 = t122 * t166 + 0.1e1;
t112 = t118 * t113;
t114 = t119 ^ 2;
t115 = 0.1e1 / t119;
t121 = qJD(2) + qJD(3);
t164 = t115 * t118;
t136 = t121 * (t112 * t115 / t114 + t164);
t129 = cos(qJ(1));
t156 = qJD(1) * t129;
t148 = t128 * t156;
t170 = 0.1e1 / t108 ^ 2 * (t122 * t136 + t148 * t166);
t177 = -0.2e1 * t170;
t161 = t121 * t128;
t106 = 0.1e1 / t108;
t149 = t118 * t156;
t150 = t116 * t161;
t92 = (-(-t119 * t161 - t149) * t115 + t113 * t150) * t106;
t142 = t92 - t161;
t160 = 0.1e1 / t128 * t129;
t146 = 0.1e1 + t166;
t176 = t128 * t146;
t127 = t129 ^ 2;
t175 = qJD(1) * (0.1e1 / t122 * t127 + 0.1e1) * t160;
t158 = t128 * t118;
t105 = atan2(-t158, -t119);
t104 = cos(t105);
t103 = sin(t105);
t151 = t103 * t158;
t100 = -t104 * t119 - t151;
t97 = 0.1e1 / t100;
t98 = 0.1e1 / t100 ^ 2;
t174 = t106 - 0.1e1;
t162 = t119 * t121;
t140 = t118 * t127 * t162;
t143 = -t128 * t92 + t121;
t167 = t104 * t118;
t87 = t143 * t167 + (t119 * t142 - t149) * t103;
t172 = t87 * t97 * t98;
t95 = t113 * t127 * t98 + 0.1e1;
t173 = (t98 * t140 + (-t127 * t172 - t148 * t98) * t113) / t95 ^ 2;
t93 = 0.1e1 / t95;
t171 = t93 * t98;
t169 = t129 * t98;
t168 = t103 * t128;
t165 = t113 * t128;
t163 = t118 * t129;
t124 = 0.1e1 / t128 ^ 2;
t159 = t124 * t127;
t157 = qJD(1) * t128;
t155 = 0.2e1 * t172;
t154 = t97 * t173;
t111 = t114 * t159 + 0.1e1;
t153 = 0.2e1 * (-t114 * t175 - t124 * t140) / t111 ^ 2;
t152 = t93 * t162;
t147 = 0.2e1 * t98 * t173;
t145 = 0.1e1 + t159;
t144 = t115 * t177;
t141 = t104 * t106 * t113 * t115;
t139 = t146 * t129;
t138 = t145 * t118;
t109 = 0.1e1 / t111;
t101 = t106 * t176;
t91 = (t103 * t118 * t174 - t128 * t141) * t129;
t90 = t118 * t153 * t160 + (qJD(1) * t138 - t160 * t162) * t109;
t89 = -t119 * t168 + t167 + (t103 * t119 - t104 * t158) * t101;
t88 = t176 * t177 + (qJD(1) * t139 + 0.2e1 * t128 * t136) * t106;
t85 = (-t97 * t93 * t157 + (-0.2e1 * t154 + (-t121 * t89 - t87) * t171) * t129) * t119 + (t89 * t129 * t147 + (-t129 * t121 * t97 - (-t104 * t128 * t88 - t142 * t103 + (-t103 * t121 - t104 * t156 + t168 * t92) * t101) * t98 * t163 + (t129 * t155 + t98 * t157) * t89 - ((t88 - t156) * t103 + (t101 * t142 + t143) * t104) * t119 * t169) * t93) * t118;
t1 = [t144 * t163 + (t121 * t139 - t157 * t164) * t106, t88, t88, 0, 0, 0; (-t97 * t152 + (0.2e1 * t154 + (qJD(1) * t91 + t87) * t171) * t118) * t128 + (-t91 * t98 * t152 + (t91 * t147 + (t91 * t155 + ((0.2e1 * t118 * t170 + t162 + (-t115 * t165 * t92 - t162) * t106) * t103 + (t144 * t165 + t118 * t92 + (t112 * t150 - (t92 - 0.2e1 * t161) * t118) * t106) * t104) * t169) * t93) * t118 + (-t97 + (-(t122 - t127) * t141 + t174 * t151) * t98) * t118 * t93 * qJD(1)) * t129, t85, t85, 0, 0, 0; t145 * t119 * t153 + (0.2e1 * t119 * t175 + t121 * t138) * t109, t90, t90, 0, 0, 0;];
JaD_rot  = t1;
