% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:02
% EndTime: 2019-12-31 17:51:04
% DurationCPUTime: 0.54s
% Computational Cost: add. (849->114), mult. (1639->132), div. (0->0), fcn. (794->6), ass. (0->86)
t101 = (qJD(1) * qJD(5));
t114 = 2 * t101;
t71 = sin(qJ(4));
t62 = t71 ^ 2;
t73 = cos(qJ(4));
t63 = t73 ^ 2;
t106 = t62 + t63;
t76 = qJD(1) ^ 2;
t47 = t106 * t76;
t67 = sin(pkin(7));
t94 = pkin(1) * t67 + qJ(3);
t68 = cos(pkin(7));
t113 = pkin(1) * t68 + pkin(2) + pkin(6);
t104 = t71 * qJDD(1);
t102 = qJD(1) * qJD(4);
t95 = t73 * t102;
t40 = -t95 - t104;
t105 = qJD(1) * t73;
t48 = qJD(4) * pkin(4) - qJ(5) * t105;
t66 = qJDD(1) * pkin(2);
t72 = sin(qJ(1));
t74 = cos(qJ(1));
t98 = t72 * g(1) - t74 * g(2);
t35 = qJDD(1) * pkin(1) + t98;
t86 = t74 * g(1) + t72 * g(2);
t36 = -t76 * pkin(1) - t86;
t93 = t68 * t35 - t67 * t36;
t15 = -t76 * qJ(3) + qJDD(3) - t66 - t93;
t13 = -qJDD(1) * pkin(6) + t15;
t64 = -g(3) + qJDD(2);
t7 = t71 * t13 + t73 * t64;
t79 = -t40 * qJ(5) + qJD(4) * t48 + t71 * t114 - t7;
t112 = t62 * t76;
t111 = t63 * t76;
t108 = t73 * t76;
t99 = t71 * t108;
t49 = qJDD(4) + t99;
t110 = t71 * t49;
t50 = qJDD(4) - t99;
t109 = t73 * t50;
t107 = t67 * t35 + t68 * t36;
t56 = t73 * qJDD(1);
t103 = qJ(5) * qJDD(1);
t100 = qJD(3) * qJD(1);
t96 = t71 * t102;
t75 = qJD(4) ^ 2;
t52 = -t75 - t111;
t92 = t52 + t112;
t60 = qJDD(1) * qJ(3);
t91 = -t76 * pkin(2) + t107 + t60;
t89 = -pkin(1) * (t67 * qJDD(1) + t68 * t76) - t107;
t41 = t56 - t96;
t9 = t73 * t13;
t88 = -qJDD(4) * pkin(4) + t41 * qJ(5) - t9;
t87 = -t76 * pkin(6) + t91;
t85 = qJ(5) * t102 + t64;
t3 = -0.2e1 * t73 * t101 + (-pkin(4) * t108 - t85) * t71 - t88;
t4 = -pkin(4) * t112 - t79;
t1 = t73 * t3 + t71 * t4;
t6 = t71 * t64 - t9;
t2 = -t73 * t6 + t71 * t7;
t51 = -t75 - t112;
t21 = t71 * t51 + t109;
t39 = 0.2e1 * t95 + t104;
t84 = -t113 * t21 + t39 * t94;
t22 = t73 * t52 - t110;
t42 = t56 - 0.2e1 * t96;
t83 = -t113 * t22 + t42 * t94;
t45 = t106 * qJDD(1);
t82 = t113 * t45 - t47 * t94;
t80 = -pkin(1) * (-t68 * qJDD(1) + t67 * t76) + t93;
t78 = t40 * pkin(4) - t48 * t105 - qJDD(5) - t87;
t57 = 0.2e1 * t100;
t77 = t57 - t78;
t46 = (-t62 + t63) * t76;
t27 = -t73 * t49 - t71 * t52;
t26 = t109 - t71 * (t75 - t111);
t25 = (t41 - t96) * t73;
t24 = -t71 * t50 + t73 * t51;
t23 = t73 * (-t75 + t112) - t110;
t20 = (-t40 + t95) * t71;
t16 = -t73 * t39 - t71 * t42;
t14 = t57 + t91;
t12 = t57 + t87;
t5 = -qJ(5) * t112 + t77;
t8 = [0, 0, 0, 0, 0, qJDD(1), t98, t86, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t80, t89, 0, pkin(1) * (t67 * t107 + t68 * t93), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t66 - t80, t57 + 0.2e1 * t60 - t89, pkin(1) * (t67 * t14 - t68 * t15) - pkin(2) * t15 + qJ(3) * t14, t25, t16, t26, t20, t23, 0, t71 * t12 + t84, t73 * t12 + t83, -t2 + t82, -t113 * t2 + t94 * t12, t25, t16, t26, t20, t23, 0, -t71 * (-pkin(4) * t39 - 0.2e1 * t100 + t78) + (-t109 - t71 * (t51 + t112)) * qJ(5) + t84, t73 * (-t92 * qJ(5) + t77) - t71 * (-pkin(4) * t42 - qJ(5) * t49) + t83, ((t114 + t103) * t73 + t88) * t73 + (t71 * t103 + t73 * t85 + t79) * t71 + t82, (t71 * pkin(4) + t94) * t5 + (-qJ(5) - t113) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, t24, t27, 0, t71 * t6 + t73 * t7, 0, 0, 0, 0, 0, 0, t24, t27, 0, -t71 * t3 + t73 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t76, t15, 0, 0, 0, 0, 0, 0, t21, t22, -t45, t2, 0, 0, 0, 0, 0, 0, t21, t22, -t45, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t46, t56, -t99, -t104, qJDD(4), -t6, -t7, 0, 0, t99, t46, t56, -t99, -t104, qJDD(4), pkin(4) * t50 + t3, t92 * pkin(4) + t79, -pkin(4) * t56, pkin(4) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t42, -t47, t5;];
tauJ_reg = t8;
