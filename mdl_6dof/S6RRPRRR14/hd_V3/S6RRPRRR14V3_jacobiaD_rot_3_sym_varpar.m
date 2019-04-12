% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14V3_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiaD_rot_3_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:04
% EndTime: 2019-04-12 15:12:05
% DurationCPUTime: 0.52s
% Computational Cost: add. (776->71), mult. (1835->159), div. (470->13), fcn. (2177->7), ass. (0->71)
t80 = sin(qJ(1));
t111 = qJD(1) * t80;
t131 = 0.2e1 * t80;
t71 = t80 ^ 2;
t82 = cos(qJ(1));
t75 = t82 ^ 2;
t76 = 0.1e1 / t82;
t129 = (t71 / t75 + 0.1e1) * t76 * t111;
t79 = sin(qJ(2));
t112 = t80 * t79;
t81 = cos(qJ(2));
t61 = atan2(-t112, -t81);
t59 = sin(t61);
t103 = t59 * t112;
t60 = cos(t61);
t55 = -t60 * t81 - t103;
t52 = 0.1e1 / t55;
t72 = 0.1e1 / t81;
t53 = 0.1e1 / t55 ^ 2;
t73 = 0.1e1 / t81 ^ 2;
t128 = -0.2e1 * t79;
t69 = t79 ^ 2;
t116 = t69 * t73;
t66 = t71 * t116 + 0.1e1;
t62 = 0.1e1 / t66;
t127 = t62 - 0.1e1;
t110 = qJD(1) * t82;
t101 = t79 * t110;
t109 = qJD(2) * t80;
t118 = t60 * t79;
t108 = qJD(2) * t81;
t48 = (-(-t80 * t108 - t101) * t72 + t109 * t116) * t62;
t44 = (-t48 * t80 + qJD(2)) * t118 + (-t101 + (t48 - t109) * t81) * t59;
t126 = t44 * t52 * t53;
t125 = t48 * t59;
t124 = t48 * t79;
t123 = t53 * t79;
t122 = t53 * t82;
t113 = t72 * t79;
t68 = t79 * t69;
t74 = t72 * t73;
t90 = qJD(2) * (t68 * t74 + t113);
t94 = t69 * t80 * t110;
t121 = (t71 * t90 + t73 * t94) / t66 ^ 2;
t100 = 0.1e1 + t116;
t58 = t100 * t80 * t62;
t120 = t58 * t80;
t119 = t59 * t81;
t117 = t69 * t72;
t115 = t69 * t75;
t114 = t71 / t82 ^ 2;
t51 = t53 * t115 + 0.1e1;
t107 = 0.2e1 * (-t115 * t126 + (t75 * t79 * t108 - t94) * t53) / t51 ^ 2;
t106 = 0.2e1 * t126;
t67 = t73 * t114 + 0.1e1;
t105 = 0.2e1 * (t74 * qJD(2) * t79 * t114 + t129 * t73) / t67 ^ 2;
t104 = t79 * t122;
t102 = t62 * t117;
t99 = 0.1e1 + t114;
t98 = t79 * t107;
t97 = t121 * t128;
t96 = t121 * t131;
t95 = t80 * t102;
t93 = t100 * t82;
t91 = t99 * t79 * t73;
t64 = 0.1e1 / t67;
t49 = 0.1e1 / t51;
t47 = (t127 * t79 * t59 - t60 * t95) * t82;
t46 = -t80 * t119 + t118 + (-t60 * t112 + t119) * t58;
t45 = -t100 * t96 + (qJD(1) * t93 + t131 * t90) * t62;
t1 = [t82 * t72 * t97 + (qJD(2) * t93 - t111 * t113) * t62, t45, 0, 0, 0, 0; (t52 * t98 + (-t52 * t108 + (qJD(1) * t47 + t44) * t123) * t49) * t80 + (t53 * t98 * t47 + (-((t127 * t108 + t48 * t95 + t97) * t59 + (t96 * t117 - t124 + (t124 + (-t68 * t73 + t128) * t109) * t62) * t60) * t104 + (t79 * t106 - t53 * t108) * t47 + (-t52 + ((-t71 + t75) * t60 * t102 + t127 * t103) * t53) * t79 * qJD(1)) * t49) * t82 (t46 * t123 - t52 * t81) * t82 * t107 + ((-t52 * t111 + (-qJD(2) * t46 - t44) * t122) * t81 + (-t82 * qJD(2) * t52 - (-t45 * t60 * t80 + t59 * t109 + t120 * t125 - t125 + (-qJD(2) * t59 - t110 * t60) * t58) * t104 + (t82 * t106 + t53 * t111) * t46 - ((t45 - t110) * t59 + ((0.1e1 - t120) * qJD(2) + (t58 - t80) * t48) * t60) * t81 * t122) * t79) * t49, 0, 0, 0, 0; t99 * t72 * t105 + (-qJD(2) * t91 - 0.2e1 * t129 * t72) * t64, t76 * t73 * t105 * t112 + ((-0.2e1 * t69 * t74 - t72) * t76 * t109 - qJD(1) * t91) * t64, 0, 0, 0, 0;];
JaD_rot  = t1;
