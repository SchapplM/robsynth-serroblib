% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobiaD_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:35
% EndTime: 2019-02-26 20:28:36
% DurationCPUTime: 0.52s
% Computational Cost: add. (396->67), mult. (1839->157), div. (436->14), fcn. (2165->7), ass. (0->74)
t88 = cos(qJ(1));
t115 = qJD(1) * t88;
t86 = sin(qJ(1));
t79 = 0.1e1 / t86 ^ 2;
t84 = t88 ^ 2;
t122 = t79 * t84;
t85 = sin(qJ(4));
t73 = t85 ^ 2;
t72 = t73 * t122 + 0.1e1;
t77 = t86 ^ 2;
t78 = 0.1e1 / t86;
t95 = (-0.1e1 - 0.1e1 / t77 * t84) * t78 * t115;
t114 = qJD(4) * t85;
t87 = cos(qJ(4));
t99 = t84 * t87 * t114;
t137 = (t73 * t95 + t79 * t99) / t72 ^ 2;
t75 = 0.1e1 / t85 ^ 2;
t82 = t87 ^ 2;
t123 = t75 * t82;
t105 = 0.1e1 + t123;
t135 = t105 * t86;
t113 = qJD(4) * t86;
t106 = t87 * t115;
t71 = t77 * t123 + 0.1e1;
t66 = 0.1e1 / t71;
t74 = 0.1e1 / t85;
t53 = ((-t85 * t113 + t106) * t74 - t113 * t123) * t66;
t134 = -t53 - t113;
t119 = t86 * t87;
t70 = atan2(t119, t85);
t64 = sin(t70);
t108 = t64 * t119;
t65 = cos(t70);
t60 = t65 * t85 + t108;
t57 = 0.1e1 / t60;
t58 = 0.1e1 / t60 ^ 2;
t133 = -0.2e1 * t87;
t132 = t66 - 0.1e1;
t120 = t82 * t86;
t100 = t115 * t120;
t121 = t82 * t84;
t124 = t65 * t87;
t49 = (t53 * t86 + qJD(4)) * t124 + (t134 * t85 + t106) * t64;
t130 = t49 * t57 * t58;
t56 = t58 * t121 + 0.1e1;
t131 = (-t121 * t130 + (-t99 - t100) * t58) / t56 ^ 2;
t129 = t53 * t87;
t128 = t58 * t87;
t127 = t58 * t88;
t81 = t87 * t82;
t96 = (t87 + 0.1e1 / t73 * t81) * t74;
t126 = (-qJD(4) * t77 * t96 + t75 * t100) / t71 ^ 2;
t125 = t64 * t86;
t118 = t87 * t88;
t117 = qJD(1) * t86;
t116 = qJD(1) * t87;
t112 = qJD(4) * t88;
t111 = -0.2e1 * t130;
t110 = 0.2e1 * t126;
t109 = t58 * t118;
t107 = t66 * t74 * t82;
t104 = 0.1e1 + t122;
t103 = t131 * t133;
t102 = -0.2e1 * t74 * t126;
t101 = t86 * t107;
t98 = t105 * t88;
t97 = t104 * t87;
t68 = 0.1e1 / t72;
t63 = t66 * t135;
t54 = 0.1e1 / t56;
t52 = (-t132 * t87 * t64 + t65 * t101) * t88;
t51 = -t85 * t125 + t124 - (t65 * t119 - t64 * t85) * t63;
t50 = t110 * t135 + (-qJD(1) * t98 + 0.2e1 * t113 * t96) * t66;
t1 = [t102 * t118 + (-t74 * t86 * t116 - qJD(4) * t98) * t66, 0, 0, t50, 0, 0; (t57 * t103 + (-t57 * t114 + (-qJD(1) * t52 - t49) * t128) * t54) * t86 + (t58 * t103 * t52 + (((-t53 * t101 + t87 * t110 + t132 * t114) * t64 + (t102 * t120 + t129 + (-t129 + (-t75 * t81 + t133) * t113) * t66) * t65) * t109 + (t87 * t111 - t58 * t114) * t52 + (t57 + ((-t77 + t84) * t65 * t107 + t132 * t108) * t58) * t116) * t54) * t88, 0, 0, 0.2e1 * (-t51 * t128 - t57 * t85) * t88 * t131 + ((-t57 * t117 + (-qJD(4) * t51 - t49) * t127) * t85 + (t57 * t112 + (t50 * t65 * t86 + t134 * t64 - (-qJD(4) * t64 + t115 * t65 - t125 * t53) * t63) * t109 + (t88 * t111 - t58 * t117) * t51 + ((-t50 - t115) * t64 + ((t63 * t86 - 0.1e1) * qJD(4) + (t63 - t86) * t53) * t65) * t85 * t127) * t87) * t54, 0, 0; -0.2e1 * t104 * t85 * t137 + (qJD(4) * t97 + 0.2e1 * t85 * t95) * t68, 0, 0, 0.2e1 * t78 * t118 * t137 + (t78 * t85 * t112 + qJD(1) * t97) * t68, 0, 0;];
JaD_rot  = t1;
