% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR10_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:59
% EndTime: 2019-02-26 20:54:00
% DurationCPUTime: 0.66s
% Computational Cost: add. (703->82), mult. (2191->190), div. (456->12), fcn. (2616->9), ass. (0->85)
t102 = sin(qJ(3));
t94 = 0.1e1 / t102 ^ 2;
t104 = cos(qJ(3));
t98 = t104 ^ 2;
t147 = t94 * t98;
t105 = cos(qJ(1));
t126 = 0.1e1 + t147;
t99 = t105 ^ 2;
t92 = t147 * t99 + 0.1e1;
t90 = 0.1e1 / t92;
t115 = t126 * t90;
t76 = t105 * t115;
t161 = t105 * t76 - 0.1e1;
t101 = cos(pkin(10));
t103 = sin(qJ(1));
t133 = qJD(3) * t104;
t122 = t103 * t133;
t100 = sin(pkin(10));
t140 = t103 * t100;
t141 = t102 * t105;
t86 = t101 * t141 - t140;
t78 = qJD(1) * t86 + t101 * t122;
t139 = t103 * t101;
t84 = t100 * t105 + t102 * t139;
t81 = 0.1e1 / t84 ^ 2;
t160 = t78 * t81;
t159 = t104 * t147;
t83 = -t101 * t105 + t102 * t140;
t149 = t81 * t83;
t79 = t83 ^ 2;
t74 = t79 * t81 + 0.1e1;
t72 = 0.1e1 / t74;
t80 = 0.1e1 / t84;
t158 = (-t100 * t80 + t101 * t149) * t72;
t138 = t105 * t104;
t89 = atan2(-t138, t102);
t87 = sin(t89);
t88 = cos(t89);
t71 = t88 * t102 - t138 * t87;
t68 = 0.1e1 / t71;
t93 = 0.1e1 / t102;
t69 = 0.1e1 / t71 ^ 2;
t157 = t90 - 0.1e1;
t135 = qJD(1) * t105;
t116 = t103 * t98 * t135;
t96 = t103 ^ 2;
t146 = t96 * t98;
t132 = qJD(3) * t105;
t137 = qJD(1) * t103;
t145 = t102 * t87;
t136 = qJD(1) * t104;
t64 = ((t102 * t132 + t103 * t136) * t93 + t132 * t147) * t90;
t59 = (-t64 + t132) * t145 + (t87 * t137 + (-t105 * t64 + qJD(3)) * t88) * t104;
t155 = t59 * t68 * t69;
t67 = t146 * t69 + 0.1e1;
t156 = (-t146 * t155 + (-t102 * t133 * t96 + t116) * t69) / t67 ^ 2;
t150 = t80 * t160;
t85 = t100 * t141 + t139;
t77 = qJD(1) * t85 + t100 * t122;
t154 = (t149 * t77 - t150 * t79) / t74 ^ 2;
t65 = 0.1e1 / t67;
t152 = t65 * t69;
t113 = qJD(3) * (-t104 - t159) * t93;
t151 = (t113 * t99 - t116 * t94) / t92 ^ 2;
t148 = t93 * t98;
t144 = t103 * t69;
t142 = qJD(3) * t76;
t134 = qJD(3) * t102;
t131 = -0.2e1 * t155;
t130 = 0.2e1 * t154;
t129 = t68 * t156;
t128 = t90 * t148;
t127 = t104 * t151;
t125 = t104 * t157;
t124 = t65 * t134;
t123 = t104 * t135;
t121 = t104 * t132;
t120 = -0.2e1 * t69 * t156;
t119 = 0.2e1 * t83 * t150;
t118 = t105 * t128;
t117 = t87 * t125;
t63 = (-t118 * t88 - t117) * t103;
t61 = -t161 * t88 * t104 + (t105 - t76) * t145;
t60 = -t115 * t137 + 0.2e1 * (t113 * t90 - t126 * t151) * t105;
t1 = [t93 * t90 * t123 + (-qJD(3) * t115 - 0.2e1 * t127 * t93) * t103, 0, t60, 0, 0, 0; (t68 * t124 + (0.2e1 * t129 + (qJD(1) * t63 + t59) * t152) * t104) * t105 + (-t63 * t69 * t124 + (t63 * t120 + (t63 * t131 + ((t118 * t64 + t134 * t157 + 0.2e1 * t127) * t87 + (-t64 * t125 + (0.2e1 * t148 * t151 + (0.2e1 * t104 + t159) * t90 * qJD(3)) * t105) * t88) * t144) * t65) * t104 + (t68 + ((t96 - t99) * t88 * t128 - t105 * t117) * t69) * t65 * t136) * t103, 0 (t68 * t65 * t135 + (-0.2e1 * t129 + (-qJD(3) * t61 - t59) * t152) * t103) * t102 + (t61 * t103 * t120 + (t103 * qJD(3) * t68 + (t103 * t131 + t135 * t69) * t61 + (((-t105 * t60 + t137 * t76) * t88 + (t161 * t64 + t132 - t142) * t87) * t104 + ((-t60 - t137) * t87 + (-t64 * t76 - qJD(3) + (t64 + t142) * t105) * t88) * t102) * t144) * t65) * t104, 0, 0, 0; (t149 * t86 - t80 * t85) * t130 + ((-qJD(1) * t83 + t100 * t121) * t80 + t86 * t119 + (-t85 * t78 - (-qJD(1) * t84 + t101 * t121) * t83 - t86 * t77) * t81) * t72, 0, -t123 * t158 + (t134 * t158 + ((-0.2e1 * t154 * t80 - t72 * t160) * t100 + (t130 * t149 + (-t77 * t81 + t119) * t72) * t101) * t104) * t103, 0, 0, 0;];
JaD_rot  = t1;
