% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPPP1_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_3_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_3_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_rot_3_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:29:48
% EndTime: 2019-02-26 19:29:48
% DurationCPUTime: 0.26s
% Computational Cost: add. (361->28), mult. (1032->88), div. (172->14), fcn. (1422->9), ass. (0->50)
t100 = cos(qJ(1));
t120 = cos(pkin(4));
t99 = sin(qJ(1));
t116 = t99 * t120;
t96 = sin(pkin(6));
t98 = cos(pkin(6));
t110 = -t100 * t98 + t96 * t116;
t79 = t110 ^ 2;
t97 = sin(pkin(4));
t88 = 0.1e1 / t97 ^ 2;
t94 = 0.1e1 / t99 ^ 2;
t74 = t79 * t94 * t88 + 0.1e1;
t135 = -0.2e1 / t74;
t93 = 0.1e1 / t99;
t133 = qJD(1) * t93 * t94;
t115 = t100 * t120;
t82 = -t96 * t115 - t99 * t98;
t76 = t82 * qJD(1);
t134 = -0.2e1 * (-t100 * t79 * t133 - t110 * t76 * t94) * t88 / t74 ^ 2;
t87 = 0.1e1 / t97;
t125 = t87 / t98;
t124 = t88 / t98 ^ 2;
t80 = -t98 * t115 + t99 * t96;
t122 = t97 * t98;
t68 = atan2(-t80, -t122);
t66 = sin(t68);
t67 = cos(t68);
t130 = t80 ^ 2;
t73 = t130 * t124 + 0.1e1;
t69 = 0.1e1 / t73;
t131 = (t67 * t80 * t125 - t66) * t69 + t66;
t64 = -t67 * t122 - t66 * t80;
t61 = 0.1e1 / t64;
t62 = 0.1e1 / t64 ^ 2;
t83 = t100 * t96 + t98 * t116;
t77 = t83 * qJD(1);
t56 = t131 * t77;
t129 = t56 * t61 * t62;
t70 = 0.1e1 / t73 ^ 2;
t128 = t70 * t80;
t75 = t80 * qJD(1);
t127 = t75 * t62;
t126 = t77 * t83;
t123 = t124 * t125;
t119 = t69 * t125;
t118 = t100 * t87 * t94;
t78 = t83 ^ 2;
t60 = t78 * t62 + 0.1e1;
t57 = t131 * t83;
t1 = [-0.2e1 * t123 * t126 * t128 - t75 * t119, 0, 0, 0; 0.2e1 * (-t57 * t62 * t83 + t61 * t80) / t60 ^ 2 * (-t83 * t127 + t78 * t129) + (-t77 * t61 + (-t80 * t56 - t57 * t75) * t62 + (0.2e1 * t57 * t129 - t131 * t127 - (-t66 * t124 * t128 + (-0.2e1 * t119 + (0.2e1 * t130 * t123 + t125) * t70) * t67) * t62 * t126) * t83) / t60, 0, 0, 0; t82 * t93 * t87 * t134 + t76 * t118 * t135 + (t100 ^ 2 * t87 * t133 * t135 + t118 * t134) * t110, 0, 0, 0;];
JaD_rot  = t1;
