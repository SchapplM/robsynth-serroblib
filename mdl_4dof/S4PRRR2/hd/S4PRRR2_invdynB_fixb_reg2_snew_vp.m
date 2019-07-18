% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRRR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynB_fixb_reg2_snew_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:25
% EndTime: 2019-07-18 13:27:27
% DurationCPUTime: 0.68s
% Computational Cost: add. (1578->100), mult. (2193->117), div. (0->0), fcn. (1476->6), ass. (0->61)
t113 = sin(qJ(2));
t116 = cos(qJ(2));
t110 = g(2) + qJDD(1);
t111 = sin(qJ(4));
t112 = sin(qJ(3));
t114 = cos(qJ(4));
t115 = cos(qJ(3));
t119 = (t111 * t115 + t112 * t114) * t110;
t131 = (t111 * t112 - t114 * t115) * t110;
t109 = qJD(2) + qJD(3);
t105 = qJD(4) + t109;
t103 = t105 ^ 2;
t108 = qJDD(2) + qJDD(3);
t104 = qJDD(4) + t108;
t85 = t114 * t103 + t111 * t104;
t88 = t111 * t103 - t114 * t104;
t61 = t112 * t88 - t115 * t85;
t64 = t112 * t85 + t115 * t88;
t52 = t113 * t64 + t116 * t61;
t153 = -qJ(1) * t52 - t113 * t119 - t116 * t131;
t149 = -t113 * t61 + t116 * t64;
t152 = qJ(1) * t149 - t113 * t131 + t116 * t119;
t107 = t109 ^ 2;
t92 = t115 * t107 + t112 * t108;
t95 = t112 * t107 - t115 * t108;
t74 = t113 * t95 - t116 * t92;
t151 = -qJ(1) * t74 + (-t112 * t113 + t115 * t116) * t110;
t101 = t116 * g(1) + t113 * g(3);
t97 = qJDD(2) * pkin(1) + t101;
t102 = t113 * g(1) - t116 * g(3);
t117 = qJD(2) ^ 2;
t98 = -t117 * pkin(1) + t102;
t76 = t112 * t98 - t115 * t97;
t69 = t108 * pkin(2) - t76;
t77 = t112 * t97 + t115 * t98;
t70 = -t107 * pkin(2) + t77;
t54 = t111 * t70 - t114 * t69;
t55 = t111 * t69 + t114 * t70;
t127 = t111 * t54 + t114 * t55;
t44 = t111 * t55 - t114 * t54;
t145 = -t112 * t44 + t115 * t127;
t40 = t112 * t127 + t115 * t44;
t150 = t113 * t40 - t116 * t145;
t37 = -t113 * t145 - t116 * t40;
t141 = t113 * t92 + t116 * t95;
t146 = qJ(1) * t141 + (t112 * t116 + t113 * t115) * t110;
t126 = t112 * t76 + t115 * t77;
t57 = t112 * t77 - t115 * t76;
t144 = t113 * t57 - t116 * t126;
t47 = -t113 * t126 - t116 * t57;
t130 = t113 * t110;
t129 = t116 * t110;
t128 = pkin(2) * t110 * t112;
t99 = t113 * qJDD(2) + t116 * t117;
t124 = qJ(1) * t99 + t129;
t100 = -t116 * qJDD(2) + t113 * t117;
t123 = qJ(1) * t100 + t130;
t78 = -t116 * t101 - t113 * t102;
t122 = t113 * t101 - t116 * t102;
t96 = (-pkin(2) * t115 - pkin(1)) * t110;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t100, t99, 0, t78, 0, 0, 0, 0, 0, 0, t141, -t74, 0, t47, 0, 0, 0, 0, 0, 0, t149, -t52, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t99, t100, 0, -t122, 0, 0, 0, 0, 0, 0, t74, t141, 0, -t144, 0, 0, 0, 0, 0, 0, t52, t149, 0, -t150; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, t110, 0, g(3), qJ(1) * g(3), 0, 0, -t99, 0, t100, 0, t124, -t123, t122, qJ(1) * t122, 0, 0, t74, 0, t141, 0, t151, -t146, t144, pkin(1) * t129 + qJ(1) * t144, 0, 0, t52, 0, t149, 0, t153, -t152, t150, qJ(1) * t150 - t113 * t128 - t116 * t96; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, -g(1), -g(3), 0, 0, 0, 0, 0, 0, 0, -qJDD(2), -t101, t102, 0, 0, 0, 0, 0, 0, 0, -t108, pkin(1) * t95 + t76, pkin(1) * t92 + t77, 0, -pkin(1) * t57, 0, 0, 0, 0, 0, -t104, pkin(1) * t64 + pkin(2) * t88 + t54, -pkin(1) * t61 + pkin(2) * t85 + t55, 0, -pkin(1) * t40 - pkin(2) * t44; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, -g(1), -qJ(1) * g(1), 0, 0, -t100, 0, -t99, 0, t123, t124, t78, qJ(1) * t78, 0, 0, -t141, 0, t74, 0, t146, t151, t47, pkin(1) * t130 + qJ(1) * t47, 0, 0, -t149, 0, t52, 0, t152, t153, t37, qJ(1) * t37 - t113 * t96 + t116 * t128;];
tauB_reg  = t1;
