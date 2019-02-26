% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR8_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobiaD_rot_4_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:07
% EndTime: 2019-02-26 20:43:07
% DurationCPUTime: 0.50s
% Computational Cost: add. (587->69), mult. (1835->162), div. (470->13), fcn. (2177->7), ass. (0->72)
t84 = cos(qJ(1));
t111 = qJD(1) * t84;
t81 = sin(qJ(3));
t70 = 0.1e1 / t81;
t133 = 0.2e1 * t70;
t132 = 0.2e1 * t84;
t108 = qJD(3) * t84;
t83 = cos(qJ(3));
t112 = qJD(1) * t83;
t82 = sin(qJ(1));
t102 = t82 * t112;
t71 = 0.1e1 / t81 ^ 2;
t78 = t83 ^ 2;
t117 = t71 * t78;
t80 = t84 ^ 2;
t69 = t80 * t117 + 0.1e1;
t66 = 0.1e1 / t69;
t50 = ((t81 * t108 + t102) * t70 + t108 * t117) * t66;
t130 = -t50 + t108;
t114 = t84 * t83;
t63 = atan2(-t114, t81);
t61 = sin(t63);
t62 = cos(t63);
t57 = -t61 * t114 + t62 * t81;
t54 = 0.1e1 / t57;
t74 = 0.1e1 / t82;
t55 = 0.1e1 / t57 ^ 2;
t129 = t66 - 0.1e1;
t109 = qJD(3) * t83;
t73 = t82 ^ 2;
t116 = t73 * t78;
t120 = t62 * t83;
t46 = (-t50 * t84 + qJD(3)) * t120 + (t130 * t81 + t102) * t61;
t127 = t46 * t54 * t55;
t53 = t55 * t116 + 0.1e1;
t96 = t78 * t82 * t111;
t128 = (-t116 * t127 + (-t73 * t81 * t109 + t96) * t55) / t53 ^ 2;
t126 = t50 * t83;
t125 = t55 * t82;
t124 = t55 * t83;
t118 = t70 * t83;
t72 = t70 * t71;
t77 = t83 * t78;
t92 = qJD(3) * (-t72 * t77 - t118);
t123 = (-t71 * t96 + t80 * t92) / t69 ^ 2;
t115 = 0.1e1 / t82 ^ 2 * t80;
t68 = t71 * t115 + 0.1e1;
t93 = (-0.1e1 - 0.1e1 / t73 * t80) * t74 * t111;
t122 = (-t72 * t109 * t115 + t71 * t93) / t68 ^ 2;
t121 = t61 * t84;
t119 = t70 * t78;
t113 = qJD(1) * t82;
t110 = qJD(3) * t81;
t107 = -0.2e1 * t127;
t106 = t83 * t128;
t105 = t82 * t124;
t104 = t83 * t123;
t103 = t66 * t119;
t101 = 0.1e1 + t117;
t100 = -0.1e1 - t115;
t99 = t123 * t132;
t98 = t84 * t103;
t97 = t129 * t83 * t61;
t95 = t101 * t82;
t94 = t100 * t83 * t71;
t64 = 0.1e1 / t68;
t60 = t101 * t84 * t66;
t51 = 0.1e1 / t53;
t49 = (-t62 * t98 - t97) * t82;
t48 = t81 * t121 + t120 + (-t62 * t114 - t61 * t81) * t60;
t47 = -t101 * t99 + (-qJD(1) * t95 + t92 * t132) * t66;
t1 = [-0.2e1 * t82 * t70 * t104 + (-qJD(3) * t95 + t111 * t118) * t66, 0, t47, 0, 0, 0; (0.2e1 * t54 * t106 + (t54 * t110 + (qJD(1) * t49 + t46) * t124) * t51) * t84 + (-0.2e1 * t55 * t106 * t49 + (((t129 * t110 + t50 * t98 + 0.2e1 * t104) * t61 + (t99 * t119 + t126 + (-t126 + (t71 * t77 + 0.2e1 * t83) * t108) * t66) * t62) * t105 + (t83 * t107 - t55 * t110) * t49 + (t54 + ((t73 - t80) * t62 * t103 - t84 * t97) * t55) * t112) * t51) * t82, 0, 0.2e1 * (-t48 * t124 - t54 * t81) * t82 * t128 + ((t54 * t111 + (-qJD(3) * t48 - t46) * t125) * t81 + (t82 * qJD(3) * t54 + (-t47 * t62 * t84 + t130 * t61 + (-qJD(3) * t61 + t113 * t62 + t121 * t50) * t60) * t105 + (t82 * t107 + t55 * t111) * t48 + ((-t47 - t113) * t61 + ((t60 * t84 - 0.1e1) * qJD(3) + (-t60 + t84) * t50) * t62) * t81 * t125) * t83) * t51, 0, 0, 0; t100 * t122 * t133 + (qJD(3) * t94 + t93 * t133) * t64, 0, -0.2e1 * t74 * t71 * t114 * t122 + ((-0.2e1 * t72 * t78 - t70) * t74 * t108 + qJD(1) * t94) * t64, 0, 0, 0;];
JaD_rot  = t1;
