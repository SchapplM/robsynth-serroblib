% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function JaD_rot = S4RPPP1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_rot_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:29:58
% EndTime: 2019-02-26 19:29:58
% DurationCPUTime: 0.24s
% Computational Cost: add. (295->28), mult. (1032->91), div. (172->14), fcn. (1422->9), ass. (0->52)
t100 = cos(qJ(1));
t121 = cos(pkin(4));
t99 = sin(qJ(1));
t116 = t99 * t121;
t96 = sin(pkin(6));
t98 = cos(pkin(6));
t110 = t100 * t96 + t98 * t116;
t78 = t110 ^ 2;
t97 = sin(pkin(4));
t91 = 0.1e1 / t97 ^ 2;
t94 = 0.1e1 / t99 ^ 2;
t74 = t78 * t94 * t91 + 0.1e1;
t133 = -0.2e1 / t74;
t120 = qJD(1) * t100;
t115 = t100 * t121;
t80 = -t98 * t115 + t99 * t96;
t75 = t80 * qJD(1);
t93 = 0.1e1 / t99;
t95 = t93 * t94;
t132 = -0.2e1 * (-t110 * t75 * t94 - t78 * t95 * t120) * t91 / t74 ^ 2;
t90 = 0.1e1 / t97;
t125 = 0.1e1 / t96 * t90;
t124 = 0.1e1 / t96 ^ 2 * t91;
t122 = t97 * t96;
t81 = t96 * t115 + t99 * t98;
t70 = atan2(-t81, t122);
t66 = sin(t70);
t67 = cos(t70);
t130 = t81 ^ 2;
t73 = t130 * t124 + 0.1e1;
t68 = 0.1e1 / t73;
t109 = (t67 * t81 * t125 + t66) * t68 - t66;
t64 = t67 * t122 - t66 * t81;
t61 = 0.1e1 / t64;
t62 = 0.1e1 / t64 ^ 2;
t112 = t96 * t116;
t77 = -qJD(1) * t112 + t98 * t120;
t56 = t109 * t77;
t129 = t56 * t61 * t62;
t69 = 0.1e1 / t73 ^ 2;
t128 = t69 * t81;
t76 = t81 * qJD(1);
t127 = t76 * t62;
t84 = t100 * t98 - t112;
t126 = t77 * t84;
t123 = t124 * t125;
t119 = t68 * t125;
t118 = t100 * t90 * t94;
t79 = t84 ^ 2;
t60 = t79 * t62 + 0.1e1;
t57 = t109 * t84;
t1 = [0.2e1 * t123 * t126 * t128 + t76 * t119, 0, 0, 0; 0.2e1 * (t57 * t62 * t84 + t61 * t81) / t60 ^ 2 * (-t84 * t127 - t79 * t129) + (-t77 * t61 + (t81 * t56 + t57 * t76) * t62 + (0.2e1 * t57 * t129 + t109 * t127 - (-t66 * t124 * t128 + (0.2e1 * t119 + (-0.2e1 * t130 * t123 - t125) * t69) * t67) * t62 * t126) * t84) / t60, 0, 0, 0; t80 * t93 * t90 * t132 + t75 * t118 * t133 + (t100 ^ 2 * t90 * t95 * qJD(1) * t133 + t118 * t132) * t110, 0, 0, 0;];
JaD_rot  = t1;
