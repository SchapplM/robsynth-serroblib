% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR12_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiaD_rot_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:12
% EndTime: 2019-02-26 21:44:13
% DurationCPUTime: 0.52s
% Computational Cost: add. (1114->72), mult. (3196->173), div. (656->14), fcn. (4222->9), ass. (0->74)
t125 = cos(qJ(2));
t126 = cos(qJ(1));
t156 = cos(pkin(6));
t141 = t126 * t156;
t139 = t125 * t141;
t123 = sin(qJ(2));
t124 = sin(qJ(1));
t152 = t123 * t124;
t104 = -t139 + t152;
t122 = sin(pkin(6));
t114 = 0.1e1 / t122;
t119 = 0.1e1 / t125;
t144 = t104 * t114 * t119;
t153 = t122 * t125;
t94 = atan2(-t104, -t153);
t92 = sin(t94);
t93 = cos(t94);
t101 = t104 ^ 2;
t115 = 0.1e1 / t122 ^ 2;
t120 = 0.1e1 / t125 ^ 2;
t99 = t101 * t115 * t120 + 0.1e1;
t95 = 0.1e1 / t99;
t166 = (t93 * t144 - t92) * t95 + t92;
t87 = -t104 * t92 - t93 * t153;
t84 = 0.1e1 / t87;
t116 = 0.1e1 / t124;
t117 = 0.1e1 / t124 ^ 2;
t85 = 0.1e1 / t87 ^ 2;
t149 = qJD(2) * t123;
t158 = t125 * t92;
t163 = t104 * t93;
t143 = t120 * t149;
t160 = t114 * t95;
t134 = -t123 * t141 - t124 * t125;
t142 = t124 * t156;
t135 = -t126 * t123 - t125 * t142;
t90 = -t135 * qJD(1) - t134 * qJD(2);
t77 = (t104 * t143 + t119 * t90) * t160;
t74 = -t77 * t163 - t92 * t90 + (t93 * t149 + t77 * t158) * t122;
t165 = t74 * t84 * t85;
t154 = t120 * t123;
t136 = t104 * t154 - t119 * t134;
t78 = t136 * t160;
t164 = t77 * t78;
t162 = t135 * t85;
t161 = t135 * t93;
t159 = t119 * t95;
t157 = t92 * t135;
t155 = t117 * t126;
t151 = t126 * t125;
t150 = qJD(1) * t126;
t102 = t135 ^ 2;
t81 = t102 * t85 + 0.1e1;
t138 = qJD(2) * t156 + qJD(1);
t88 = -qJD(1) * t139 - qJD(2) * t151 + t138 * t152;
t148 = 0.2e1 * (-t102 * t165 + t88 * t162) / t81 ^ 2;
t147 = 0.2e1 * t165;
t121 = t119 * t120;
t146 = -0.2e1 * (t101 * t121 * t149 + t104 * t120 * t90) * t115 / t99 ^ 2;
t140 = t123 * t142;
t108 = -t140 + t151;
t103 = t108 ^ 2;
t100 = t103 * t115 * t117 + 0.1e1;
t118 = t116 * t117;
t89 = t134 * qJD(1) + t135 * qJD(2);
t145 = 0.2e1 * (-t103 * t118 * t150 + t108 * t117 * t89) * t115 / t100 ^ 2;
t133 = t119 * t146 + t95 * t143;
t97 = 0.1e1 / t100;
t91 = -qJD(1) * t140 - t124 * t149 + t138 * t151;
t79 = 0.1e1 / t81;
t76 = t166 * t135;
t75 = -t78 * t163 + t92 * t134 + (t123 * t93 + t78 * t158) * t122;
t73 = (t136 * t146 + (t90 * t154 + t119 * t91 + (-t134 * t154 + (0.2e1 * t121 * t123 ^ 2 + t119) * t104) * qJD(2)) * t95) * t114;
t1 = [(-t133 * t135 - t88 * t159) * t114, t73, 0, 0, 0, 0; t104 * t84 * t148 + (-t90 * t84 + (t104 * t74 + t76 * t88) * t85) * t79 - (t76 * t147 * t79 + (t76 * t148 + ((t77 * t95 * t144 + t146) * t157 + ((t95 - 0.1e1) * t77 + (-t133 * t104 - t90 * t159) * t114) * t161 - t166 * t88) * t79) * t85) * t135 (-t108 * t84 - t75 * t162) * t148 + (-t75 * t135 * t147 + t89 * t84 + (-t108 * t74 + t75 * t88 + (t104 * t164 - t91) * t157 + (-t104 * t73 + t134 * t77 - t78 * t90) * t161) * t85 + ((-qJD(2) * t78 - t77) * t92 * t123 + (t73 * t92 + (qJD(2) + t164) * t93) * t125) * t122 * t162) * t79, 0, 0, 0, 0; ((t108 * t155 - t116 * t134) * t145 + (-t89 * t155 - t116 * t91 + (-t134 * t155 + (0.2e1 * t118 * t126 ^ 2 + t116) * t108) * qJD(1)) * t97) * t114 (t116 * t88 * t97 - (t117 * t97 * t150 + t116 * t145) * t135) * t114, 0, 0, 0, 0;];
JaD_rot  = t1;
