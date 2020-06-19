% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S2RR3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
%
% Output:
% m_new_reg [(3*3)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S2RR3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynm_fixb_reg2_snew_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynm_fixb_reg2_snew_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynm_fixb_reg2_snew_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynm_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:28
% EndTime: 2020-06-19 09:14:29
% DurationCPUTime: 0.55s
% Computational Cost: add. (435->61), mult. (695->64), div. (0->0), fcn. (476->4), ass. (0->39)
t110 = sin(qJ(1));
t112 = cos(qJ(1));
t108 = qJD(1) + qJD(2);
t106 = t108 ^ 2;
t107 = qJDD(1) + qJDD(2);
t109 = sin(qJ(2));
t111 = cos(qJ(2));
t95 = t111 * t106 + t109 * t107;
t98 = t109 * t106 - t111 * t107;
t126 = t110 * t95 + t112 * t98;
t128 = pkin(3) * t98 - t109 * g(3);
t91 = pkin(3) * t95 - t111 * g(3);
t134 = pkin(2) * t126 + t110 * t91 + t112 * t128;
t127 = t110 * t98 - t112 * t95;
t133 = pkin(2) * t127 + t110 * t128 - t112 * t91;
t104 = t112 * g(1) + t110 * g(2);
t113 = qJD(1) ^ 2;
t100 = -t113 * pkin(1) - t104;
t103 = t110 * g(1) - t112 * g(2);
t114 = qJDD(1) * pkin(1) + t103;
t86 = t109 * t100 - t111 * t114;
t87 = t111 * t100 + t109 * t114;
t119 = t109 * t86 + t111 * t87;
t78 = t109 * t87 - t111 * t86;
t120 = t112 * t78;
t130 = t110 * t119 + t120;
t121 = t110 * t78;
t129 = t112 * t119 - t121;
t118 = -t110 * t103 - t112 * t104;
t102 = t112 * qJDD(1) - t110 * t113;
t117 = -pkin(2) * t102 - t110 * g(3);
t115 = t112 * t103 - t110 * t104;
t101 = t110 * qJDD(1) + t112 * t113;
t93 = -pkin(2) * t101 + t112 * g(3);
t82 = -pkin(1) * t95 - t87;
t81 = -pkin(1) * t98 - t86;
t76 = pkin(1) * t78;
t75 = pkin(1) * g(3) + pkin(3) * t119;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t102, 0, -t101, 0, t117, -t93, -t115, -pkin(2) * t115, 0, 0, -t126, 0, t127, 0, t134, -t133, -t130, -pkin(2) * t130 - pkin(3) * t120 - t110 * t75; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t101, 0, t102, 0, t93, t117, t118, pkin(2) * t118, 0, 0, -t127, 0, -t126, 0, t133, t134, t129, pkin(2) * t129 - pkin(3) * t121 + t112 * t75; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t103, t104, 0, 0, 0, 0, 0, 0, 0, t107, t81, t82, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t113, 0, 0, -g(3), -t103, 0, 0, 0, -t98, 0, -t95, 0, t128, t91, -t78, -pkin(3) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, qJDD(1), 0, g(3), 0, -t104, 0, 0, 0, t95, 0, -t98, 0, -t91, t128, t119, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t103, t104, 0, 0, 0, 0, 0, 0, 0, t107, t81, t82, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, -t106, 0, 0, -g(3), t86, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, t107, 0, g(3), 0, t87, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, -t86, -t87, 0, 0;];
m_new_reg = t1;
