% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:22:36
% EndTime: 2019-02-26 21:22:36
% DurationCPUTime: 0.65s
% Computational Cost: add. (1897->83), mult. (2191->183), div. (456->12), fcn. (2616->9), ass. (0->85)
t107 = qJ(2) + pkin(9);
t106 = cos(t107);
t104 = t106 ^ 2;
t105 = sin(t107);
t153 = 0.1e1 / t105 ^ 2 * t104;
t112 = sin(qJ(1));
t127 = 0.1e1 + t153;
t108 = t112 ^ 2;
t99 = t108 * t153 + 0.1e1;
t97 = 0.1e1 / t99;
t124 = t127 * t97;
t80 = t112 * t124;
t169 = t112 * t80 - 0.1e1;
t100 = 0.1e1 / t105;
t166 = t106 * t153;
t122 = qJD(2) * (-t106 - t166) * t100;
t113 = cos(qJ(1));
t145 = qJD(1) * t113;
t133 = t112 * t145;
t168 = (t108 * t122 + t133 * t153) / t99 ^ 2;
t110 = sin(pkin(10));
t142 = qJD(2) * t113;
t131 = t106 * t142;
t149 = t112 * t110;
t111 = cos(pkin(10));
t151 = t111 * t113;
t93 = -t105 * t149 + t151;
t85 = t93 * qJD(1) + t110 * t131;
t148 = t112 * t111;
t152 = t110 * t113;
t91 = t105 * t152 + t148;
t88 = 0.1e1 / t91 ^ 2;
t167 = t85 * t88;
t150 = t112 * t106;
t96 = atan2(-t150, t105);
t94 = sin(t96);
t136 = t94 * t150;
t95 = cos(t96);
t78 = t105 * t95 - t136;
t75 = 0.1e1 / t78;
t87 = 0.1e1 / t91;
t76 = 0.1e1 / t78 ^ 2;
t165 = t97 - 0.1e1;
t143 = qJD(2) * t112;
t158 = t105 * t94;
t71 = ((t105 * t143 - t106 * t145) * t100 + t143 * t153) * t97;
t66 = (-t71 + t143) * t158 + (-t94 * t145 + (-t112 * t71 + qJD(2)) * t95) * t106;
t164 = t66 * t75 * t76;
t109 = t113 ^ 2;
t74 = t104 * t109 * t76 + 0.1e1;
t72 = 0.1e1 / t74;
t162 = t72 * t76;
t161 = t75 * t72;
t160 = t87 * t167;
t123 = t105 * t151 - t149;
t159 = t88 * t123;
t156 = t113 * t76;
t155 = qJD(2) * t80;
t154 = t100 * t104;
t147 = qJD(1) * t106;
t146 = qJD(1) * t112;
t144 = qJD(2) * t105;
t134 = t76 * t144;
t141 = 0.2e1 * (-t109 * t106 * t134 + (-t109 * t164 - t76 * t133) * t104) / t74 ^ 2;
t140 = 0.2e1 * t164;
t86 = t123 ^ 2;
t83 = t86 * t88 + 0.1e1;
t92 = t105 * t148 + t152;
t84 = t92 * qJD(1) - t111 * t131;
t139 = 0.2e1 * (-t84 * t159 - t86 * t160) / t83 ^ 2;
t138 = 0.2e1 * t168;
t137 = t100 * t112 * t97;
t135 = t165 * t106;
t132 = t106 * t143;
t130 = t75 * t141;
t129 = t76 * t141;
t128 = -0.2e1 * t123 * t160;
t126 = t106 * t138;
t125 = t104 * t137;
t81 = 0.1e1 / t83;
t121 = (-t110 * t159 + t111 * t87) * t81;
t70 = (t95 * t125 + t94 * t135) * t113;
t68 = -t169 * t95 * t106 + (t112 - t80) * t158;
t67 = t124 * t145 + 0.2e1 * (t122 * t97 - t127 * t168) * t112;
t1 = [t137 * t147 + (qJD(2) * t124 + t100 * t126) * t113, t67, 0, 0, 0, 0; (t144 * t161 + (t130 + (qJD(1) * t70 + t66) * t162) * t106) * t112 + (t70 * t129 * t106 + (t70 * t134 + (t70 * t140 + ((t71 * t125 + t165 * t144 + t126) * t94 + (-t71 * t135 + (t138 * t154 + (0.2e1 * t106 + t166) * t97 * qJD(2)) * t112) * t95) * t156) * t106 + (-t75 + (t165 * t136 - (-t108 + t109) * t97 * t95 * t154) * t76) * t147) * t72) * t113 (t146 * t161 + (t130 + (qJD(2) * t68 + t66) * t162) * t113) * t105 + (t68 * t113 * t129 + (-t75 * t142 + (t113 * t140 + t76 * t146) * t68 + (-((-t112 * t67 - t145 * t80) * t95 + (t169 * t71 + t143 - t155) * t94) * t106 - ((-t67 + t145) * t94 + (-t71 * t80 - qJD(2) + (t71 + t155) * t112) * t95) * t105) * t156) * t72) * t106, 0, 0, 0, 0; (-t93 * t159 - t87 * t92) * t139 + ((t123 * qJD(1) + t111 * t132) * t87 + t93 * t128 + (-t92 * t85 + (-t91 * qJD(1) - t110 * t132) * t123 - t93 * t84) * t88) * t81, t105 * t121 * t142 + (t121 * t146 + ((t87 * t139 + t81 * t167) * t111 + (-t139 * t159 + (-t84 * t88 + t128) * t81) * t110) * t113) * t106, 0, 0, 0, 0;];
JaD_rot  = t1;
