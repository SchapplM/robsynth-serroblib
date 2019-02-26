% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR5
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
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR5_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:48
% EndTime: 2019-02-26 21:30:48
% DurationCPUTime: 0.26s
% Computational Cost: add. (351->41), mult. (876->109), div. (130->12), fcn. (1072->9), ass. (0->57)
t104 = cos(pkin(6));
t103 = sin(pkin(6));
t108 = cos(qJ(1));
t128 = t108 * t103;
t91 = atan2(-t128, -t104);
t89 = sin(t91);
t90 = cos(t91);
t75 = -t90 * t104 - t89 * t128;
t72 = 0.1e1 / t75;
t107 = cos(qJ(2));
t126 = t108 * t107;
t105 = sin(qJ(2));
t106 = sin(qJ(1));
t130 = t106 * t105;
t117 = t104 * t130 - t126;
t82 = 0.1e1 / t117;
t98 = 0.1e1 / t104;
t99 = 0.1e1 / t104 ^ 2;
t73 = 0.1e1 / t75 ^ 2;
t83 = 0.1e1 / t117 ^ 2;
t127 = t108 * t105;
t129 = t106 * t107;
t119 = t104 * t129 + t127;
t81 = t119 ^ 2;
t138 = t81 * t83;
t84 = t82 * t83;
t137 = t81 * t84;
t118 = t104 * t127 + t129;
t136 = t118 * t119;
t102 = t108 ^ 2;
t97 = t103 ^ 2;
t94 = t102 * t97 * t99 + 0.1e1;
t93 = 0.1e1 / t94 ^ 2;
t135 = t93 * t103 * t97;
t134 = t97 * t98;
t101 = t106 ^ 2;
t133 = t101 * t73;
t132 = t106 * t73;
t131 = t108 * t73;
t125 = qJD(1) * t108;
t85 = -t104 * t126 + t130;
t76 = t85 * qJD(1) + t117 * qJD(2);
t123 = t119 * t83 * t76;
t77 = t118 * qJD(1) + t119 * qJD(2);
t80 = 0.1e1 + t138;
t124 = 0.2e1 * (-t77 * t137 - t123) / t80 ^ 2;
t92 = 0.1e1 / t94;
t122 = t92 * t134;
t121 = t90 * t122;
t120 = (-t92 + 0.1e1) * t89 * t103;
t68 = (t108 * t121 + t120) * t106;
t100 = t98 * t99;
t78 = 0.1e1 / t80;
t74 = t72 * t73;
t71 = t97 * t133 + 0.1e1;
t67 = qJD(1) * t68;
t1 = [(-0.2e1 * t100 * t101 * t135 - t103 * t92 * t98) * t125, 0, 0, 0, 0, 0; (0.2e1 * (t108 * t72 - t68 * t132) / t71 ^ 2 * (-t101 * t67 * t74 + t125 * t132) * t97 + ((-0.2e1 * t106 * t68 * t74 + t131) * t67 + (t68 * t131 + (t102 * t73 * t121 + t72 + t120 * t131 + (-t108 * t89 * t99 * t135 + (-0.2e1 * t122 + (0.2e1 * t100 * t102 * t97 ^ 2 + t134) * t93) * t90) * t133) * t106) * qJD(1)) / t71) * t103, 0, 0, 0, 0, 0; (t83 * t136 + t82 * t85) * t124 + (-(t119 * qJD(1) + t118 * qJD(2)) * t82 + 0.2e1 * t84 * t77 * t136 + (t85 * t77 + (t117 * qJD(1) + t85 * qJD(2)) * t119 + t118 * t76) * t83) * t78 (t117 * t82 + t138) * t124 + (0.2e1 * t123 + (t117 * t83 + 0.2e1 * t137 - t82) * t77) * t78, 0, 0, 0, 0;];
JaD_rot  = t1;
