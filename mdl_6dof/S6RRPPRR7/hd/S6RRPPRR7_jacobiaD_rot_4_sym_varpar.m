% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR7_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:57
% EndTime: 2019-02-26 21:31:57
% DurationCPUTime: 0.25s
% Computational Cost: add. (351->43), mult. (876->111), div. (130->12), fcn. (1072->9), ass. (0->58)
t106 = cos(pkin(6));
t101 = 0.1e1 / t106 ^ 2;
t110 = cos(qJ(1));
t104 = t110 ^ 2;
t105 = sin(pkin(6));
t99 = t105 ^ 2;
t95 = t104 * t99 * t101 + 0.1e1;
t94 = 0.1e1 / t95 ^ 2;
t139 = t101 * t94;
t128 = t110 * t105;
t92 = atan2(-t128, -t106);
t90 = sin(t92);
t91 = cos(t92);
t76 = -t91 * t106 - t90 * t128;
t73 = 0.1e1 / t76;
t107 = sin(qJ(2));
t127 = t110 * t107;
t108 = sin(qJ(1));
t109 = cos(qJ(2));
t129 = t108 * t109;
t88 = t106 * t129 + t127;
t82 = 0.1e1 / t88;
t100 = 0.1e1 / t106;
t74 = 0.1e1 / t76 ^ 2;
t83 = 0.1e1 / t88 ^ 2;
t84 = t82 * t83;
t126 = t110 * t109;
t130 = t108 * t107;
t119 = t106 * t130 - t126;
t85 = t119 ^ 2;
t138 = t84 * t85;
t137 = t85 * t83;
t122 = t106 * t126;
t86 = t122 - t130;
t136 = t86 * t119;
t135 = t100 * t99;
t134 = t100 * t139;
t103 = t108 ^ 2;
t133 = t103 * t74;
t132 = t108 * t74;
t131 = t110 * t74;
t125 = qJD(1) * t110;
t87 = -t106 * t127 - t129;
t78 = t87 * qJD(1) - t88 * qJD(2);
t123 = t119 * t83 * t78;
t77 = -qJD(1) * t122 - qJD(2) * t126 + (qJD(2) * t106 + qJD(1)) * t130;
t81 = 0.1e1 + t137;
t124 = 0.2e1 * (t77 * t138 - t123) / t81 ^ 2;
t93 = 0.1e1 / t95;
t121 = t91 * t93 * t135;
t120 = (-t93 + 0.1e1) * t90 * t105;
t69 = (t110 * t121 + t120) * t108;
t98 = t105 * t99;
t79 = 0.1e1 / t81;
t75 = t73 * t74;
t72 = t99 * t133 + 0.1e1;
t68 = qJD(1) * t69;
t1 = [(-t100 * t105 * t93 - 0.2e1 * t103 * t98 * t134) * t125, 0, 0, 0, 0, 0; (0.2e1 * (t110 * t73 - t69 * t132) / t72 ^ 2 * (-t103 * t68 * t75 + t125 * t132) * t99 + ((-0.2e1 * t108 * t69 * t75 + t131) * t68 + (t69 * t131 + (t104 * t74 * t121 + t73 + t120 * t131 + (-t110 * t90 * t98 * t139 + (0.2e1 * t104 * t99 ^ 2 * t134 + (-0.2e1 * t93 + t94) * t135) * t91) * t133) * t108) * qJD(1)) / t72) * t105, 0, 0, 0, 0, 0; (-t83 * t136 - t82 * t87) * t124 + ((t119 * qJD(1) - t86 * qJD(2)) * t82 + 0.2e1 * t84 * t77 * t136 + (t87 * t77 + (-t88 * qJD(1) + t87 * qJD(2)) * t119 - t86 * t78) * t83) * t79 (t82 * t88 + t137) * t124 + (0.2e1 * t123 + (-t83 * t88 - 0.2e1 * t138 + t82) * t77) * t79, 0, 0, 0, 0;];
JaD_rot  = t1;
