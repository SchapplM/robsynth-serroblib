% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:16
% EndTime: 2019-02-26 21:38:17
% DurationCPUTime: 0.59s
% Computational Cost: add. (4795->72), mult. (2858->158), div. (686->14), fcn. (3330->7), ass. (0->75)
t125 = sin(qJ(1));
t119 = t125 ^ 2;
t117 = qJ(2) + pkin(10) + qJ(4);
t115 = sin(t117);
t110 = t115 ^ 2;
t116 = cos(t117);
t113 = 0.1e1 / t116 ^ 2;
t161 = t110 * t113;
t105 = t119 * t161 + 0.1e1;
t109 = t115 * t110;
t111 = t116 ^ 2;
t112 = 0.1e1 / t116;
t118 = qJD(2) + qJD(4);
t159 = t112 * t115;
t133 = t118 * (t109 * t112 / t111 + t159);
t126 = cos(qJ(1));
t151 = qJD(1) * t126;
t143 = t125 * t151;
t166 = 0.1e1 / t105 ^ 2 * (t119 * t133 + t143 * t161);
t175 = -0.2e1 * t166;
t103 = 0.1e1 / t105;
t141 = 0.1e1 + t161;
t173 = t125 * t141;
t98 = t103 * t173;
t174 = t125 * t98 - 0.1e1;
t155 = 0.1e1 / t125 * t126;
t124 = t126 ^ 2;
t172 = qJD(1) * (0.1e1 / t119 * t124 + 0.1e1) * t155;
t153 = t125 * t115;
t102 = atan2(-t153, -t116);
t101 = cos(t102);
t100 = sin(t102);
t146 = t100 * t153;
t97 = -t101 * t116 - t146;
t94 = 0.1e1 / t97;
t95 = 0.1e1 / t97 ^ 2;
t171 = t103 - 0.1e1;
t157 = t116 * t118;
t137 = t115 * t124 * t157;
t144 = t115 * t151;
t156 = t118 * t125;
t162 = t101 * t115;
t145 = t113 * t156;
t89 = (-(-t116 * t156 - t144) * t112 + t110 * t145) * t103;
t84 = (-t125 * t89 + t118) * t162 + (-t144 + (t89 - t156) * t116) * t100;
t169 = t84 * t94 * t95;
t92 = t110 * t124 * t95 + 0.1e1;
t170 = (t95 * t137 + (-t124 * t169 - t95 * t143) * t110) / t92 ^ 2;
t90 = 0.1e1 / t92;
t168 = t90 * t95;
t165 = t118 * t98;
t163 = t126 * t95;
t160 = t110 * t125;
t158 = t115 * t126;
t121 = 0.1e1 / t125 ^ 2;
t154 = t121 * t124;
t152 = qJD(1) * t125;
t150 = 0.2e1 * t169;
t108 = t111 * t154 + 0.1e1;
t149 = 0.2e1 / t108 ^ 2 * (-t111 * t172 - t121 * t137);
t148 = t94 * t170;
t147 = t90 * t157;
t142 = 0.2e1 * t95 * t170;
t140 = 0.1e1 + t154;
t139 = t112 * t175;
t138 = t101 * t103 * t110 * t112;
t136 = t141 * t126;
t135 = t140 * t115;
t106 = 0.1e1 / t108;
t88 = (t171 * t115 * t100 - t125 * t138) * t126;
t87 = t115 * t149 * t155 + (qJD(1) * t135 - t155 * t157) * t106;
t86 = -t174 * t162 + (-t125 + t98) * t116 * t100;
t85 = t173 * t175 + (qJD(1) * t136 + 0.2e1 * t125 * t133) * t103;
t82 = (-t94 * t90 * t152 + (-0.2e1 * t148 + (-t118 * t86 - t84) * t168) * t126) * t116 + (t86 * t126 * t142 + (-t126 * t118 * t94 - ((-t125 * t85 - t151 * t98) * t101 + (t174 * t89 + t156 - t165) * t100) * t95 * t158 + (t126 * t150 + t95 * t152) * t86 - ((t85 - t151) * t100 + (t89 * t98 + t118 + (-t89 - t165) * t125) * t101) * t116 * t163) * t90) * t115;
t1 = [t139 * t158 + (t118 * t136 - t152 * t159) * t103, t85, 0, t85, 0, 0; (-t94 * t147 + (0.2e1 * t148 + (qJD(1) * t88 + t84) * t168) * t115) * t125 + (-t88 * t95 * t147 + (t88 * t142 + (t88 * t150 + ((0.2e1 * t115 * t166 + t157 + (-t112 * t89 * t160 - t157) * t103) * t100 + (t139 * t160 + t115 * t89 + (t109 * t145 - (t89 - 0.2e1 * t156) * t115) * t103) * t101) * t163) * t90) * t115 + (-t94 + (-(t119 - t124) * t138 + t171 * t146) * t95) * t115 * t90 * qJD(1)) * t126, t82, 0, t82, 0, 0; t140 * t116 * t149 + (0.2e1 * t116 * t172 + t118 * t135) * t106, t87, 0, t87, 0, 0;];
JaD_rot  = t1;
