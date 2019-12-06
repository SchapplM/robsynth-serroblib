% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:58
% EndTime: 2019-12-05 16:42:02
% DurationCPUTime: 0.61s
% Computational Cost: add. (1469->116), mult. (2171->137), div. (0->0), fcn. (1318->8), ass. (0->86)
t79 = qJD(2) + qJD(3);
t77 = t79 ^ 2;
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t68 = t88 * t77 * t91;
t61 = qJDD(4) - t68;
t120 = t91 * t61;
t82 = t88 ^ 2;
t127 = t82 * t77;
t94 = qJD(4) ^ 2;
t64 = t94 + t127;
t41 = -t88 * t64 + t120;
t113 = qJD(4) * t79;
t78 = qJDD(2) + qJDD(3);
t122 = t88 * t78;
t51 = 0.2e1 * t91 * t113 + t122;
t89 = sin(qJ(3));
t92 = cos(qJ(3));
t137 = pkin(2) * (t89 * t41 + t92 * t51);
t136 = pkin(3) * t51 + pkin(7) * t41;
t134 = qJ(5) * t51;
t86 = sin(pkin(8));
t87 = cos(pkin(8));
t62 = t86 * g(1) - t87 * g(2);
t63 = -t87 * g(1) - t86 * g(2);
t90 = sin(qJ(2));
t93 = cos(qJ(2));
t100 = -t90 * t62 - t93 * t63;
t35 = -qJD(2) ^ 2 * pkin(2) - t100;
t107 = t93 * t62 - t90 * t63;
t98 = qJDD(2) * pkin(2) + t107;
t24 = t92 * t35 + t89 * t98;
t22 = -t77 * pkin(3) + t78 * pkin(7) + t24;
t84 = -g(3) + qJDD(1);
t73 = t91 * t84;
t15 = t88 * t22 - t73;
t16 = t91 * t22 + t88 * t84;
t4 = t88 * t15 + t91 * t16;
t83 = t91 ^ 2;
t126 = t83 * t77;
t36 = t120 + t88 * (-t94 + t126);
t114 = t88 * qJ(5);
t102 = -t91 * pkin(4) - t114;
t129 = t102 * t77;
t12 = -qJDD(4) * pkin(4) - t94 * qJ(5) + (t22 + t129) * t88 + qJDD(5) - t73;
t133 = 2 * qJD(5);
t103 = qJDD(4) * qJ(5) + (qJD(4) * t133) + t91 * t129 + t16;
t11 = -t94 * pkin(4) + t103;
t2 = t91 * t11 + t88 * t12;
t128 = t79 * t88;
t110 = t88 * t113;
t119 = t91 * t78;
t23 = -t89 * t35 + t92 * t98;
t21 = -t78 * pkin(3) - t77 * pkin(7) - t23;
t96 = t21 - (-t110 + t119) * pkin(4) - t134;
t9 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t128 + t96;
t132 = -pkin(3) * t9 + pkin(7) * t2;
t130 = -pkin(3) * t21 + pkin(7) * t4;
t60 = qJDD(4) + t68;
t125 = t88 * t60;
t121 = t91 * t51;
t66 = -t94 - t126;
t40 = t91 * t66 - t125;
t52 = -0.2e1 * t110 + t119;
t117 = pkin(3) * t52 + pkin(7) * t40;
t115 = t82 + t83;
t54 = t115 * t78;
t57 = t115 * t77;
t116 = pkin(3) * t57 + pkin(7) * t54;
t112 = t88 * t21 - t136;
t111 = -t91 * t21 + t117;
t106 = t88 * (qJ(5) * t57 + t12) + t91 * ((t57 - t94) * pkin(4) + t103) + t116;
t105 = t116 + t4;
t26 = t88 * t52 + t121;
t101 = t88 * t61 + t91 * t64;
t95 = t128 * t133 - t96;
t99 = pkin(4) * t121 + t88 * (-pkin(4) * t110 + t134 + t95) + t136;
t97 = t52 * t114 + t117 + t91 * ((t52 - t110) * pkin(4) + t95);
t58 = (t82 - t83) * t77;
t39 = t125 + t91 * (t94 - t127);
t37 = t91 * t60 + t88 * t66;
t31 = t51 * t88;
t30 = t52 * t91;
t28 = pkin(2) * (t89 * t54 + t92 * t57);
t25 = pkin(2) * (t89 * t40 + t92 * t52);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, t37, -t101, 0, -t91 * t15 + t88 * t16, 0, 0, 0, 0, 0, 0, t37, 0, t101, t88 * t11 - t91 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t107, t100, 0, 0, 0, 0, 0, 0, 0, t78, pkin(2) * (-t89 * t77 + t92 * t78) + t23, pkin(2) * (-t92 * t77 - t89 * t78) - t24, 0, pkin(2) * (t92 * t23 + t89 * t24), t31, t26, t39, t30, t36, 0, t25 + t111, t112 - t137, t28 + t105, pkin(2) * (-t92 * t21 + t89 * t4) + t130, t31, t39, -t26, 0, -t36, t30, t25 + t97, t28 + t106, t99 + t137, pkin(2) * t89 * t2 + (-pkin(2) * t92 + t102) * t9 + t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t23, -t24, 0, 0, t31, t26, t39, t30, t36, 0, t111, t112, t105, t130, t31, t39, -t26, 0, -t36, t30, t97, t106, t99, t102 * t9 + t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t58, t122, t68, t119, qJDD(4), -t15, -t16, 0, 0, -t68, t122, -t58, qJDD(4), -t119, t68, pkin(4) * t60 + qJ(5) * t66 - t12, (-pkin(4) * t88 + qJ(5) * t91) * t78, qJ(5) * t61 + (t64 - t94) * pkin(4) + t103, -pkin(4) * t12 + qJ(5) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t122, -t64, t12;];
tauJ_reg = t1;
