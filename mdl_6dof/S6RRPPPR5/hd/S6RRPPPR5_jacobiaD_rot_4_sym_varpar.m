% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR5_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiaD_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:13
% EndTime: 2019-02-26 21:24:13
% DurationCPUTime: 0.55s
% Computational Cost: add. (906->74), mult. (3124->186), div. (647->14), fcn. (4112->9), ass. (0->89)
t119 = cos(qJ(2));
t118 = sin(qJ(1));
t159 = cos(pkin(9));
t136 = t118 * t159;
t116 = sin(pkin(9));
t120 = cos(qJ(1));
t154 = t120 * t116;
t102 = t119 * t154 - t136;
t135 = t120 * t159;
t155 = t118 * t116;
t103 = t119 * t135 + t155;
t177 = qJD(1) * t103;
t117 = sin(qJ(2));
t109 = 0.1e1 / t117;
t110 = 0.1e1 / t117 ^ 2;
t111 = t109 * t110;
t176 = qJD(2) * (0.2e1 * t111 * t119 ^ 2 + t109);
t107 = 0.1e1 / t116;
t99 = t119 * t155 + t135;
t164 = t107 * t99;
t145 = t109 * t164;
t156 = t117 * t116;
t91 = atan2(-t99, t156);
t87 = sin(t91);
t88 = cos(t91);
t108 = 0.1e1 / t116 ^ 2;
t96 = t99 ^ 2;
t94 = t108 * t110 * t96 + 0.1e1;
t89 = 0.1e1 / t94;
t175 = (t88 * t145 + t87) * t89 - t87;
t168 = t87 * t99;
t83 = t88 * t156 - t168;
t80 = 0.1e1 / t83;
t72 = t175 * t102;
t174 = 0.2e1 * t72;
t113 = 0.1e1 / t120;
t114 = 0.1e1 / t120 ^ 2;
t81 = 0.1e1 / t83 ^ 2;
t151 = qJD(2) * t117;
t140 = t116 * t151;
t84 = t99 * qJD(1) + t120 * t140;
t169 = t81 * t84;
t150 = qJD(2) * t119;
t161 = t117 * t87;
t167 = t88 * t99;
t142 = t110 * t150;
t130 = t99 * t142;
t86 = t102 * qJD(1) - t118 * t140;
t73 = (-t109 * t86 + t130) * t89 * t107;
t70 = -t73 * t167 - t86 * t87 + (t88 * t150 - t73 * t161) * t116;
t172 = t70 * t80 * t81;
t97 = t102 ^ 2;
t77 = t81 * t97 + 0.1e1;
t173 = (-t102 * t169 - t97 * t172) / t77 ^ 2;
t115 = t113 * t114;
t141 = t111 * t150;
t153 = qJD(1) * t118;
t98 = t103 ^ 2;
t162 = t114 * t98;
t101 = -t119 * t136 + t154;
t129 = t159 * t151;
t85 = t101 * qJD(1) - t120 * t129;
t95 = t110 * t162 + 0.1e1;
t171 = (-t141 * t162 + (t103 * t114 * t85 + t115 * t98 * t153) * t110) / t95 ^ 2;
t170 = (t110 * t86 * t99 - t96 * t141) * t108 / t94 ^ 2;
t166 = t102 * t87;
t165 = t102 * t88;
t92 = 0.1e1 / t95;
t163 = t113 * t92;
t158 = t110 * t119;
t79 = (t158 * t164 + t118) * t89;
t160 = t118 - t79;
t157 = t114 * t118;
t152 = qJD(1) * t120;
t149 = 0.2e1 * t171;
t148 = -0.2e1 * t170;
t147 = t81 * t173;
t146 = t81 * t166;
t144 = t92 * t157;
t139 = 0.2e1 * t80 * t173;
t138 = t99 * t148;
t137 = -0.2e1 * t109 * t171;
t133 = t109 * t144;
t132 = t158 * t163;
t128 = t144 * t158;
t75 = 0.1e1 / t77;
t71 = -t79 * t167 + (t119 * t88 + t160 * t161) * t116;
t69 = t118 * t148 + t89 * t152 + (t138 * t158 + (t86 * t158 - t99 * t176) * t89) * t107;
t1 = [(t109 * t84 * t89 + (0.2e1 * t109 * t170 + t89 * t142) * t102) * t107, t69, 0, 0, 0, 0; t99 * t139 + (-t86 * t80 + (t70 * t99 + t72 * t84) * t81) * t75 + (t147 * t174 + (t172 * t174 - (-t73 * t89 * t145 + t148) * t146 - ((t89 - 0.1e1) * t73 + (-t89 * t130 + (t86 * t89 + t138) * t109) * t107) * t81 * t165 + t175 * t169) * t75) * t102, t71 * t75 * t169 + (-(-t69 * t167 + (t73 * t168 - t86 * t88) * t79) * t81 * t75 + 0.2e1 * (t75 * t172 + t147) * t71) * t102 + (t120 * t139 * t117 + ((-t120 * qJD(2) * t80 - (t160 * qJD(2) - t73) * t146) * t119 + (t80 * t153 + (t120 * t70 - (-t69 + t152) * t166 - (t160 * t73 - qJD(2)) * t165) * t81) * t117) * t75) * t116, 0, 0, 0, 0; t85 * t133 + (qJD(1) * t133 - qJD(2) * t132 + t113 * t137) * t101 + ((t118 * t129 - t177) * t163 + (0.2e1 * t118 ^ 2 * t115 + t113) * t92 * t177) * t109 + (-qJD(2) * t128 + t137 * t157) * t103, -t85 * t132 + t159 * t149 + (-qJD(1) * t128 + (t149 * t158 + t92 * t176) * t113) * t103, 0, 0, 0, 0;];
JaD_rot  = t1;
