% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR5_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:22
% EndTime: 2019-02-26 22:18:22
% DurationCPUTime: 0.58s
% Computational Cost: add. (3078->72), mult. (2858->158), div. (686->14), fcn. (3330->7), ass. (0->75)
t122 = sin(qJ(1));
t115 = t122 ^ 2;
t121 = qJ(2) + qJ(3);
t112 = sin(t121);
t107 = t112 ^ 2;
t113 = cos(t121);
t110 = 0.1e1 / t113 ^ 2;
t156 = t107 * t110;
t102 = t115 * t156 + 0.1e1;
t106 = t112 * t107;
t108 = t113 ^ 2;
t109 = 0.1e1 / t113;
t114 = qJD(2) + qJD(3);
t155 = t109 * t112;
t131 = t114 * (t106 * t109 / t108 + t155);
t123 = cos(qJ(1));
t147 = qJD(1) * t123;
t139 = t122 * t147;
t164 = 0.1e1 / t102 ^ 2 * (t115 * t131 + t139 * t156);
t172 = -0.2e1 * t164;
t100 = 0.1e1 / t102;
t137 = 0.1e1 + t156;
t170 = t122 * t137;
t95 = t100 * t170;
t171 = t122 * t95 - 0.1e1;
t151 = 0.1e1 / t122 * t123;
t120 = t123 ^ 2;
t169 = qJD(1) * (0.1e1 / t115 * t120 + 0.1e1) * t151;
t149 = t122 * t112;
t99 = atan2(-t149, -t113);
t97 = sin(t99);
t143 = t97 * t149;
t98 = cos(t99);
t94 = -t113 * t98 - t143;
t91 = 0.1e1 / t94;
t92 = 0.1e1 / t94 ^ 2;
t153 = t113 * t114;
t135 = t112 * t120 * t153;
t152 = t114 * t122;
t162 = t113 * t97;
t140 = t110 * t152;
t86 = (-(-t112 * t147 - t113 * t152) * t109 + t107 * t140) * t100;
t81 = (t86 - t152) * t162 + (-t97 * t147 + (-t122 * t86 + t114) * t98) * t112;
t167 = t81 * t91 * t92;
t89 = t107 * t120 * t92 + 0.1e1;
t168 = (t92 * t135 + (-t120 * t167 - t139 * t92) * t107) / t89 ^ 2;
t87 = 0.1e1 / t89;
t165 = t87 * t92;
t163 = t112 * t97;
t161 = t114 * t95;
t159 = t123 * t92;
t158 = t98 * t112;
t157 = t107 * t109;
t154 = t112 * t123;
t117 = 0.1e1 / t122 ^ 2;
t150 = t117 * t120;
t148 = qJD(1) * t122;
t146 = 0.2e1 * t167;
t105 = t108 * t150 + 0.1e1;
t145 = 0.2e1 / t105 ^ 2 * (-t108 * t169 - t117 * t135);
t144 = t91 * t168;
t142 = t87 * t153;
t141 = t122 * t157;
t138 = 0.2e1 * t92 * t168;
t136 = 0.1e1 + t150;
t134 = t137 * t123;
t133 = t136 * t112;
t130 = -t141 * t98 + t163;
t103 = 0.1e1 / t105;
t85 = (t100 * t130 - t163) * t123;
t84 = t112 * t145 * t151 + (qJD(1) * t133 - t151 * t153) * t103;
t83 = (-t122 + t95) * t162 - t171 * t158;
t82 = t170 * t172 + (qJD(1) * t134 + 0.2e1 * t122 * t131) * t100;
t79 = (-t91 * t87 * t148 + (-0.2e1 * t144 + (-t114 * t83 - t81) * t165) * t123) * t113 + (t83 * t123 * t138 + (-t123 * t114 * t91 - ((-t122 * t82 - t147 * t95) * t98 + (t171 * t86 + t152 - t161) * t97) * t92 * t154 + (t123 * t146 + t148 * t92) * t83 - ((t82 - t147) * t97 + (t86 * t95 + t114 + (-t86 - t161) * t122) * t98) * t113 * t159) * t87) * t112;
t1 = [t109 * t154 * t172 + (t114 * t134 - t148 * t155) * t100, t82, t82, 0, 0, 0; (-t91 * t142 + (0.2e1 * t144 + (qJD(1) * t85 + t81) * t165) * t112) * t122 + (-t85 * t92 * t142 + (t85 * t138 + (t85 * t146 + (t86 * t158 + t97 * t153 + 0.2e1 * t130 * t164 + ((-t141 * t86 - t153) * t97 + (t106 * t140 - (t86 - 0.2e1 * t152) * t112) * t98) * t100) * t159) * t87) * t112 + (-t91 + (-t143 + (t143 - (t115 - t120) * t98 * t157) * t100) * t92) * t112 * t87 * qJD(1)) * t123, t79, t79, 0, 0, 0; t136 * t113 * t145 + (0.2e1 * t113 * t169 + t114 * t133) * t103, t84, t84, 0, 0, 0;];
JaD_rot  = t1;
