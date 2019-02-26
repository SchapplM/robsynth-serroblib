% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:01
% EndTime: 2019-02-26 20:28:01
% DurationCPUTime: 0.67s
% Computational Cost: add. (514->81), mult. (2191->187), div. (456->12), fcn. (2616->9), ass. (0->88)
t102 = sin(qJ(4));
t94 = 0.1e1 / t102 ^ 2;
t104 = cos(qJ(4));
t98 = t104 ^ 2;
t151 = t94 * t98;
t103 = sin(qJ(1));
t126 = 0.1e1 + t151;
t96 = t103 ^ 2;
t92 = t96 * t151 + 0.1e1;
t89 = 0.1e1 / t92;
t115 = t126 * t89;
t76 = t103 * t115;
t163 = t103 * t76 - 0.1e1;
t101 = cos(pkin(9));
t105 = cos(qJ(1));
t133 = qJD(4) * t105;
t122 = t104 * t133;
t100 = sin(pkin(9));
t141 = t105 * t100;
t143 = t103 * t101;
t84 = -t102 * t143 - t141;
t78 = qJD(1) * t84 + t101 * t122;
t140 = t105 * t101;
t144 = t103 * t100;
t86 = t102 * t140 - t144;
t81 = 0.1e1 / t86 ^ 2;
t162 = t78 * t81;
t161 = t104 * t151;
t142 = t103 * t104;
t91 = atan2(t142, t102);
t87 = sin(t91);
t127 = t87 * t142;
t88 = cos(t91);
t74 = t88 * t102 + t127;
t70 = 0.1e1 / t74;
t80 = 0.1e1 / t86;
t93 = 0.1e1 / t102;
t71 = 0.1e1 / t74 ^ 2;
t160 = t89 - 0.1e1;
t135 = qJD(4) * t103;
t137 = qJD(1) * t105;
t149 = t102 * t87;
t64 = ((-t102 * t135 + t104 * t137) * t93 - t135 * t151) * t89;
t59 = (-t64 - t135) * t149 + (t87 * t137 + (t103 * t64 + qJD(4)) * t88) * t104;
t159 = t59 * t70 * t71;
t85 = t102 * t141 + t143;
t153 = t81 * t85;
t154 = t80 * t162;
t79 = t85 ^ 2;
t73 = t79 * t81 + 0.1e1;
t83 = -t102 * t144 + t140;
t77 = qJD(1) * t83 + t100 * t122;
t158 = (t77 * t153 - t79 * t154) / t73 ^ 2;
t99 = t105 ^ 2;
t150 = t98 * t99;
t67 = t71 * t150 + 0.1e1;
t65 = 0.1e1 / t67;
t156 = t65 * t71;
t114 = (t104 + t161) * t93;
t116 = t103 * t98 * t137;
t155 = (-qJD(4) * t114 * t96 + t116 * t94) / t92 ^ 2;
t152 = t89 * t93;
t147 = t105 * t71;
t146 = qJD(4) * t76;
t145 = qJD(4) * t89;
t139 = qJD(1) * t103;
t138 = qJD(1) * t104;
t136 = qJD(4) * t102;
t134 = qJD(4) * t104;
t132 = -0.2e1 * (-t150 * t159 + (-t102 * t99 * t134 - t116) * t71) / t67 ^ 2;
t131 = -0.2e1 * t159;
t130 = 0.2e1 * t158;
t129 = 0.2e1 * t155;
t128 = t103 * t152;
t125 = t160 * t104;
t124 = t65 * t136;
t123 = t103 * t134;
t121 = t70 * t132;
t120 = t71 * t132;
t119 = -0.2e1 * t93 * t155;
t118 = 0.2e1 * t85 * t154;
t117 = t98 * t128;
t68 = 0.1e1 / t73;
t113 = (-t100 * t80 + t101 * t153) * t68;
t63 = (t88 * t117 - t87 * t125) * t105;
t61 = -t163 * t88 * t104 + (-t103 + t76) * t149;
t60 = -t115 * t137 + (0.2e1 * t114 * t145 + t126 * t129) * t103;
t1 = [-t128 * t138 + (-qJD(4) * t115 + t104 * t119) * t105, 0, 0, t60, 0, 0; (-t70 * t124 + (t121 + (-qJD(1) * t63 - t59) * t156) * t104) * t103 + (-t63 * t71 * t124 + (t63 * t120 + (t63 * t131 + ((t104 * t129 - t64 * t117 + t160 * t136) * t87 + (-t64 * t125 + (t98 * t119 + (-0.2e1 * t104 - t161) * t145) * t103) * t88) * t147) * t65) * t104 + (t70 + ((-t96 + t99) * t98 * t88 * t152 + t160 * t127) * t71) * t65 * t138) * t105, 0, 0 (-t70 * t65 * t139 + (t121 + (-qJD(4) * t61 - t59) * t156) * t105) * t102 + (t61 * t105 * t120 + (t70 * t133 + (t105 * t131 - t71 * t139) * t61 + (((t103 * t60 - t137 * t76) * t88 + (t163 * t64 - t135 + t146) * t87) * t104 + ((-t60 - t137) * t87 + (t64 * t76 - qJD(4) + (-t64 + t146) * t103) * t88) * t102) * t147) * t65) * t104, 0, 0; (t84 * t153 - t80 * t83) * t130 + ((-qJD(1) * t85 - t100 * t123) * t80 + t84 * t118 + (-t83 * t78 - (-qJD(1) * t86 - t101 * t123) * t85 - t84 * t77) * t81) * t68, 0, 0, t102 * t113 * t133 + (t113 * t139 + ((-0.2e1 * t80 * t158 - t68 * t162) * t100 + (t130 * t153 + (-t77 * t81 + t118) * t68) * t101) * t105) * t104, 0, 0;];
JaD_rot  = t1;
