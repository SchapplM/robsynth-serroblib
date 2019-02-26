% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR3_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiaD_rot_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:00
% EndTime: 2019-02-26 19:46:00
% DurationCPUTime: 0.26s
% Computational Cost: add. (543->41), mult. (1611->102), div. (383->14), fcn. (2156->9), ass. (0->52)
t119 = cos(pkin(6));
t94 = sin(pkin(10));
t113 = t94 * t119;
t118 = cos(pkin(10));
t96 = sin(qJ(2));
t97 = cos(qJ(2));
t87 = -t96 * t113 + t118 * t97;
t82 = 0.1e1 / t87 ^ 2;
t95 = sin(pkin(6));
t105 = t95 ^ 2;
t126 = t82 * t94 ^ 2 * t105;
t120 = t95 * t97;
t112 = t119 * t118;
t83 = -t97 * t112 + t94 * t96;
t70 = atan2(-t83, -t120);
t68 = sin(t70);
t69 = cos(t70);
t66 = -t69 * t120 - t68 * t83;
t63 = 0.1e1 / t66;
t91 = 0.1e1 / t97;
t64 = 0.1e1 / t66 ^ 2;
t92 = 0.1e1 / t97 ^ 2;
t111 = t68 * t120 - t69 * t83;
t116 = t69 * t95 * t96;
t122 = t92 * t96;
t115 = t83 * t122;
t80 = t83 ^ 2;
t90 = 0.1e1 / t105;
t75 = t80 * t90 * t92 + 0.1e1;
t72 = 0.1e1 / t75;
t89 = 0.1e1 / t95;
t123 = t72 * t89;
t85 = t96 * t112 + t94 * t97;
t77 = t85 * qJD(2);
t58 = (qJD(2) * t115 + t77 * t91) * t123;
t56 = qJD(2) * t116 + t111 * t58 - t68 * t77;
t125 = t56 * t63 * t64;
t109 = -t97 * t113 - t118 * t96;
t124 = t64 * t109;
t78 = t109 * qJD(2);
t114 = t109 / t87 * t82 * t78;
t110 = t85 * t91 + t115;
t93 = t91 * t92;
t81 = t109 ^ 2;
t79 = t87 * qJD(2);
t76 = qJD(2) * t83;
t74 = 0.1e1 + t126;
t62 = t81 * t64 + 0.1e1;
t59 = t110 * t123;
t57 = t111 * t59 - t68 * t85 + t116;
t55 = (-0.2e1 * t110 / t75 ^ 2 * (qJD(2) * t80 * t93 * t96 + t77 * t83 * t92) * t90 + (t77 * t122 - t76 * t91 + (t85 * t122 + (0.2e1 * t93 * t96 ^ 2 + t91) * t83) * qJD(2)) * t72) * t89;
t1 = [0, t55, 0, 0, 0, 0; 0, 0.2e1 * (-t57 * t124 - t63 * t87) / t62 ^ 2 * (-t79 * t124 - t81 * t125) + (t78 * t63 + (-t87 * t56 - t57 * t79) * t64 - (0.2e1 * t57 * t125 + (-(qJD(2) * t120 - t55 * t83 - t59 * t77 + (t59 * t120 - t85) * t58) * t69 - (t58 * t59 * t83 + t76 + (t55 * t97 + (-qJD(2) * t59 - t58) * t96) * t95) * t68) * t64) * t109) / t62, 0, 0, 0, 0; 0 (0.2e1 / t74 ^ 2 * t114 * t126 + (-t79 * t82 - 0.2e1 * t114) / t74) * t94 * t95, 0, 0, 0, 0;];
JaD_rot  = t1;
