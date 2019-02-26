% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR10_jacobiaD_rot_2_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_rot_2_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_rot_2_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobiaD_rot_2_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:48
% EndTime: 2019-02-26 22:35:48
% DurationCPUTime: 0.28s
% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
t99 = sin(pkin(6));
t93 = t99 ^ 2;
t100 = cos(pkin(6));
t95 = 0.1e1 / t100 ^ 2;
t104 = cos(qJ(1));
t98 = t104 ^ 2;
t89 = t98 * t93 * t95 + 0.1e1;
t102 = sin(qJ(1));
t97 = t102 ^ 2;
t126 = 0.1e1 / t89 ^ 2 * t97;
t131 = t126 * t95;
t122 = t104 * t99;
t88 = atan2(t122, t100);
t84 = sin(t88);
t85 = cos(t88);
t72 = t85 * t100 + t84 * t122;
t67 = 0.1e1 / t72;
t103 = cos(qJ(2));
t118 = t104 * t103;
t101 = sin(qJ(2));
t121 = t102 * t101;
t113 = t100 * t121 - t118;
t77 = 0.1e1 / t113;
t94 = 0.1e1 / t100;
t68 = 0.1e1 / t72 ^ 2;
t78 = 0.1e1 / t113 ^ 2;
t119 = t104 * t101;
t120 = t102 * t103;
t81 = -t100 * t119 - t120;
t82 = t100 * t120 + t119;
t71 = t81 * qJD(1) - t82 * qJD(2);
t128 = t71 * t77 * t78;
t115 = t100 * t118;
t70 = -qJD(1) * t115 - qJD(2) * t118 + (qJD(2) * t100 + qJD(1)) * t121;
t129 = t70 * t78;
t76 = t82 ^ 2;
t75 = t76 * t78 + 0.1e1;
t130 = (t76 * t128 - t82 * t129) / t75 ^ 2;
t127 = t81 * t82;
t125 = t93 * t94;
t124 = t102 * t68;
t123 = t104 * t68;
t117 = qJD(1) * t104;
t86 = 0.1e1 / t89;
t116 = (t86 - 0.1e1) * t99;
t114 = -0.2e1 * t94 * t131;
t80 = t115 - t121;
t63 = (-t104 * t85 * t86 * t125 + t84 * t116) * t102;
t92 = t99 * t93;
t73 = 0.1e1 / t75;
t69 = t67 * t68;
t66 = t97 * t93 * t68 + 0.1e1;
t62 = qJD(1) * t63;
t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0, 0;];
JaD_rot  = t1;
