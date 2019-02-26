% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR12_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_jacobiaD_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:02
% EndTime: 2019-02-26 20:55:03
% DurationCPUTime: 0.49s
% Computational Cost: add. (585->68), mult. (1839->158), div. (436->14), fcn. (2165->7), ass. (0->74)
t85 = sin(qJ(1));
t116 = qJD(1) * t85;
t87 = cos(qJ(1));
t137 = 0.2e1 * t87;
t77 = t85 ^ 2;
t82 = 0.1e1 / t87 ^ 2;
t119 = t77 * t82;
t84 = sin(qJ(3));
t72 = t84 ^ 2;
t70 = t72 * t119 + 0.1e1;
t80 = t87 ^ 2;
t81 = 0.1e1 / t87;
t95 = (t77 / t80 + 0.1e1) * t81 * t116;
t113 = qJD(3) * t84;
t86 = cos(qJ(3));
t98 = t77 * t86 * t113;
t136 = -0.2e1 * (t72 * t95 + t82 * t98) / t70 ^ 2;
t111 = qJD(3) * t87;
t115 = qJD(1) * t86;
t105 = t85 * t115;
t74 = 0.1e1 / t84 ^ 2;
t79 = t86 ^ 2;
t121 = t74 * t79;
t71 = t80 * t121 + 0.1e1;
t68 = 0.1e1 / t71;
t73 = 0.1e1 / t84;
t52 = ((t84 * t111 + t105) * t73 + t111 * t121) * t68;
t134 = -t52 + t111;
t117 = t87 * t86;
t65 = atan2(-t117, t84);
t63 = sin(t65);
t64 = cos(t65);
t59 = -t63 * t117 + t64 * t84;
t56 = 0.1e1 / t59;
t57 = 0.1e1 / t59 ^ 2;
t133 = t68 - 0.1e1;
t120 = t77 * t79;
t124 = t64 * t86;
t48 = (-t52 * t87 + qJD(3)) * t124 + (t134 * t84 + t105) * t63;
t131 = t48 * t56 * t57;
t55 = t57 * t120 + 0.1e1;
t114 = qJD(1) * t87;
t99 = t79 * t85 * t114;
t132 = (-t120 * t131 + (-t98 + t99) * t57) / t55 ^ 2;
t130 = t52 * t86;
t129 = t57 * t85;
t128 = t57 * t86;
t122 = t73 * t86;
t78 = t86 * t79;
t94 = qJD(3) * (-t73 / t72 * t78 - t122);
t127 = (-t74 * t99 + t80 * t94) / t71 ^ 2;
t125 = t63 * t87;
t123 = t73 * t79;
t118 = t85 * t86;
t112 = qJD(3) * t85;
t110 = -0.2e1 * t131;
t109 = t86 * t132;
t108 = t57 * t118;
t107 = t86 * t127;
t106 = t68 * t123;
t104 = 0.1e1 + t121;
t103 = 0.1e1 + t119;
t102 = t127 * t137;
t101 = t87 * t106;
t100 = t133 * t86 * t63;
t97 = t104 * t85;
t96 = t103 * t86;
t66 = 0.1e1 / t70;
t62 = t104 * t87 * t68;
t53 = 0.1e1 / t55;
t51 = (-t64 * t101 - t100) * t85;
t50 = t84 * t125 + t124 + (-t64 * t117 - t63 * t84) * t62;
t49 = -t104 * t102 + (-qJD(1) * t97 + t94 * t137) * t68;
t1 = [-0.2e1 * t85 * t73 * t107 + (-qJD(3) * t97 + t114 * t122) * t68, 0, t49, 0, 0, 0; (0.2e1 * t56 * t109 + (t56 * t113 + (qJD(1) * t51 + t48) * t128) * t53) * t87 + (-0.2e1 * t57 * t109 * t51 + (((t52 * t101 + t133 * t113 + 0.2e1 * t107) * t63 + (t102 * t123 + t130 + (-t130 + (t74 * t78 + 0.2e1 * t86) * t111) * t68) * t64) * t108 + (t86 * t110 - t57 * t113) * t51 + (t56 + ((t77 - t80) * t64 * t106 - t87 * t100) * t57) * t115) * t53) * t85, 0, 0.2e1 * (-t50 * t128 - t56 * t84) * t85 * t132 + ((t56 * t114 + (-qJD(3) * t50 - t48) * t129) * t84 + (t56 * t112 + (-t49 * t64 * t87 + t134 * t63 + (-qJD(3) * t63 + t116 * t64 + t125 * t52) * t62) * t108 + (t85 * t110 + t57 * t114) * t50 + ((-t49 - t116) * t63 + ((t62 * t87 - 0.1e1) * qJD(3) + (-t62 + t87) * t52) * t64) * t84 * t129) * t86) * t53, 0, 0, 0; t103 * t84 * t136 + (qJD(3) * t96 + 0.2e1 * t84 * t95) * t66, 0, t81 * t118 * t136 + (-t81 * t84 * t112 + qJD(1) * t96) * t66, 0, 0, 0;];
JaD_rot  = t1;
