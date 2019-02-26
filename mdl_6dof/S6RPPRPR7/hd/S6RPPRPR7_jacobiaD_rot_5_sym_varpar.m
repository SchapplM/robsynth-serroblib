% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:09
% EndTime: 2019-02-26 20:29:10
% DurationCPUTime: 0.66s
% Computational Cost: add. (1897->83), mult. (2191->188), div. (456->12), fcn. (2616->9), ass. (0->86)
t104 = pkin(9) + qJ(4);
t103 = cos(t104);
t101 = t103 ^ 2;
t102 = sin(t104);
t152 = t101 / t102 ^ 2;
t110 = cos(qJ(1));
t128 = 0.1e1 + t152;
t106 = t110 ^ 2;
t96 = t106 * t152 + 0.1e1;
t94 = 0.1e1 / t96;
t148 = t110 * t94;
t77 = t128 * t148;
t167 = t110 * t77 - 0.1e1;
t108 = cos(pkin(10));
t109 = sin(qJ(1));
t137 = qJD(4) * t109;
t126 = t103 * t137;
t107 = sin(pkin(10));
t144 = t109 * t107;
t145 = t108 * t110;
t90 = t102 * t145 - t144;
t82 = t90 * qJD(1) + t108 * t126;
t143 = t109 * t108;
t146 = t107 * t110;
t88 = t102 * t143 + t146;
t85 = 0.1e1 / t88 ^ 2;
t166 = t82 * t85;
t165 = t103 * t152;
t87 = t102 * t144 - t145;
t154 = t85 * t87;
t83 = t87 ^ 2;
t80 = t83 * t85 + 0.1e1;
t78 = 0.1e1 / t80;
t84 = 0.1e1 / t88;
t164 = (-t107 * t84 + t108 * t154) * t78;
t142 = t110 * t103;
t93 = atan2(-t142, t102);
t91 = sin(t93);
t92 = cos(t93);
t75 = t92 * t102 - t91 * t142;
t72 = 0.1e1 / t75;
t97 = 0.1e1 / t102;
t73 = 0.1e1 / t75 ^ 2;
t163 = t94 - 0.1e1;
t105 = t109 ^ 2;
t139 = qJD(1) * t110;
t127 = t109 * t139;
t138 = qJD(4) * t102;
t129 = t73 * t138;
t136 = qJD(4) * t110;
t140 = qJD(1) * t109;
t151 = t102 * t91;
t68 = ((t102 * t136 + t103 * t140) * t97 + t136 * t152) * t94;
t63 = (-t68 + t136) * t151 + (t91 * t140 + (-t110 * t68 + qJD(4)) * t92) * t103;
t161 = t63 * t72 * t73;
t71 = t101 * t105 * t73 + 0.1e1;
t162 = (-t105 * t103 * t129 + (-t105 * t161 + t73 * t127) * t101) / t71 ^ 2;
t155 = t84 * t166;
t89 = t102 * t146 + t143;
t81 = t89 * qJD(1) + t107 * t126;
t160 = (t81 * t154 - t83 * t155) / t80 ^ 2;
t69 = 0.1e1 / t71;
t158 = t69 * t73;
t157 = t72 * t69;
t118 = qJD(4) * (-t103 - t165) * t97;
t156 = (t106 * t118 - t127 * t152) / t96 ^ 2;
t153 = t101 * t97;
t150 = t109 * t73;
t147 = qJD(4) * t77;
t141 = qJD(1) * t103;
t135 = -0.2e1 * t161;
t134 = 0.2e1 * t160;
t133 = t72 * t162;
t132 = t103 * t156;
t131 = t97 * t148;
t130 = t103 * t163;
t125 = t103 * t136;
t124 = -0.2e1 * t73 * t162;
t123 = 0.2e1 * t87 * t155;
t122 = t101 * t131;
t121 = t91 * t130;
t120 = t128 * t94;
t67 = (-t92 * t122 - t121) * t109;
t65 = -t167 * t92 * t103 + (t110 - t77) * t151;
t64 = -t120 * t140 + 0.2e1 * (t118 * t94 - t128 * t156) * t110;
t1 = [t131 * t141 + (-qJD(4) * t120 - 0.2e1 * t97 * t132) * t109, 0, 0, t64, 0, 0; (t138 * t157 + (0.2e1 * t133 + (qJD(1) * t67 + t63) * t158) * t103) * t110 + (t67 * t124 * t103 + (-t67 * t129 + (t67 * t135 + ((t68 * t122 + t163 * t138 + 0.2e1 * t132) * t91 + (-t68 * t130 + (0.2e1 * t153 * t156 + (0.2e1 * t103 + t165) * t94 * qJD(4)) * t110) * t92) * t150) * t103 + (t72 + (-t110 * t121 + (t105 - t106) * t94 * t92 * t153) * t73) * t141) * t69) * t109, 0, 0 (t139 * t157 + (-0.2e1 * t133 + (-qJD(4) * t65 - t63) * t158) * t109) * t102 + (t65 * t109 * t124 + (t72 * t137 + (t109 * t135 + t73 * t139) * t65 + (((-t110 * t64 + t140 * t77) * t92 + (t167 * t68 + t136 - t147) * t91) * t103 + ((-t64 - t140) * t91 + (-t68 * t77 - qJD(4) + (t68 + t147) * t110) * t92) * t102) * t150) * t69) * t103, 0, 0; (t90 * t154 - t84 * t89) * t134 + ((-t87 * qJD(1) + t107 * t125) * t84 + t90 * t123 + (-t89 * t82 - (-t88 * qJD(1) + t108 * t125) * t87 - t90 * t81) * t85) * t78, 0, 0, t102 * t137 * t164 + (-t139 * t164 + ((-0.2e1 * t84 * t160 - t78 * t166) * t107 + (t134 * t154 + (-t81 * t85 + t123) * t78) * t108) * t109) * t103, 0, 0;];
JaD_rot  = t1;
