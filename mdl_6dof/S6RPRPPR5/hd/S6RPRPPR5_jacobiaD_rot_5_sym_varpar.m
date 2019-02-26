% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:24
% EndTime: 2019-02-26 20:41:25
% DurationCPUTime: 0.68s
% Computational Cost: add. (1897->83), mult. (2191->183), div. (456->12), fcn. (2616->9), ass. (0->85)
t104 = pkin(9) + qJ(3);
t103 = cos(t104);
t101 = t103 ^ 2;
t102 = sin(t104);
t154 = t101 / t102 ^ 2;
t109 = sin(qJ(1));
t130 = 0.1e1 + t154;
t105 = t109 ^ 2;
t96 = t105 * t154 + 0.1e1;
t94 = 0.1e1 / t96;
t121 = t130 * t94;
t77 = t109 * t121;
t166 = t109 * t77 - 0.1e1;
t163 = t103 * t154;
t97 = 0.1e1 / t102;
t119 = qJD(3) * (-t103 - t163) * t97;
t110 = cos(qJ(1));
t142 = qJD(1) * t110;
t129 = t109 * t142;
t165 = (t105 * t119 + t129 * t154) / t96 ^ 2;
t107 = sin(pkin(10));
t139 = qJD(3) * t110;
t127 = t103 * t139;
t146 = t109 * t107;
t108 = cos(pkin(10));
t148 = t108 * t110;
t90 = -t102 * t146 + t148;
t82 = t90 * qJD(1) + t107 * t127;
t145 = t109 * t108;
t149 = t107 * t110;
t88 = t102 * t149 + t145;
t85 = 0.1e1 / t88 ^ 2;
t164 = t82 * t85;
t147 = t109 * t103;
t93 = atan2(-t147, t102);
t91 = sin(t93);
t133 = t91 * t147;
t92 = cos(t93);
t75 = t102 * t92 - t133;
t72 = 0.1e1 / t75;
t84 = 0.1e1 / t88;
t73 = 0.1e1 / t75 ^ 2;
t162 = t94 - 0.1e1;
t140 = qJD(3) * t109;
t153 = t102 * t91;
t68 = ((t102 * t140 - t103 * t142) * t97 + t140 * t154) * t94;
t63 = (-t68 + t140) * t153 + (-t91 * t142 + (-t109 * t68 + qJD(3)) * t92) * t103;
t161 = t63 * t72 * t73;
t106 = t110 ^ 2;
t71 = t101 * t106 * t73 + 0.1e1;
t69 = 0.1e1 / t71;
t159 = t69 * t73;
t158 = t72 * t69;
t157 = t84 * t164;
t120 = t102 * t148 - t146;
t156 = t85 * t120;
t155 = t101 * t97;
t151 = t110 * t73;
t150 = qJD(3) * t77;
t144 = qJD(1) * t103;
t143 = qJD(1) * t109;
t141 = qJD(3) * t102;
t131 = t73 * t141;
t138 = 0.2e1 * (-t106 * t103 * t131 + (-t106 * t161 - t73 * t129) * t101) / t71 ^ 2;
t137 = 0.2e1 * t161;
t83 = t120 ^ 2;
t80 = t83 * t85 + 0.1e1;
t89 = t102 * t145 + t149;
t81 = t89 * qJD(1) - t108 * t127;
t136 = 0.2e1 * (-t81 * t156 - t83 * t157) / t80 ^ 2;
t135 = 0.2e1 * t165;
t134 = t109 * t94 * t97;
t132 = t162 * t103;
t128 = t103 * t140;
t126 = t72 * t138;
t125 = t73 * t138;
t124 = -0.2e1 * t120 * t157;
t123 = t103 * t135;
t122 = t101 * t134;
t78 = 0.1e1 / t80;
t118 = (-t107 * t156 + t108 * t84) * t78;
t67 = (t92 * t122 + t91 * t132) * t110;
t65 = -t166 * t92 * t103 + (t109 - t77) * t153;
t64 = t121 * t142 + 0.2e1 * (t119 * t94 - t130 * t165) * t109;
t1 = [t134 * t144 + (qJD(3) * t121 + t97 * t123) * t110, 0, t64, 0, 0, 0; (t141 * t158 + (t126 + (qJD(1) * t67 + t63) * t159) * t103) * t109 + (t67 * t125 * t103 + (t67 * t131 + (t67 * t137 + ((t68 * t122 + t162 * t141 + t123) * t91 + (-t68 * t132 + (t135 * t155 + (0.2e1 * t103 + t163) * t94 * qJD(3)) * t109) * t92) * t151) * t103 + (-t72 + (t162 * t133 - (-t105 + t106) * t94 * t92 * t155) * t73) * t144) * t69) * t110, 0 (t143 * t158 + (t126 + (qJD(3) * t65 + t63) * t159) * t110) * t102 + (t65 * t110 * t125 + (-t72 * t139 + (t110 * t137 + t73 * t143) * t65 + (-((-t109 * t64 - t142 * t77) * t92 + (t166 * t68 + t140 - t150) * t91) * t103 - ((-t64 + t142) * t91 + (-t68 * t77 - qJD(3) + (t68 + t150) * t109) * t92) * t102) * t151) * t69) * t103, 0, 0, 0; (-t90 * t156 - t84 * t89) * t136 + ((t120 * qJD(1) + t108 * t128) * t84 + t90 * t124 + (-t89 * t82 + (-t88 * qJD(1) - t107 * t128) * t120 - t90 * t81) * t85) * t78, 0, t102 * t118 * t139 + (t118 * t143 + ((t84 * t136 + t78 * t164) * t108 + (-t136 * t156 + (-t81 * t85 + t124) * t78) * t107) * t110) * t103, 0, 0, 0;];
JaD_rot  = t1;
