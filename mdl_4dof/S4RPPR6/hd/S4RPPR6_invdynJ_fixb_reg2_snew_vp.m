% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPPR6
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPPR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:44
% EndTime: 2019-12-31 16:40:46
% DurationCPUTime: 0.67s
% Computational Cost: add. (840->135), mult. (2067->175), div. (0->0), fcn. (1244->6), ass. (0->93)
t62 = sin(pkin(6));
t63 = cos(pkin(6));
t64 = sin(qJ(4));
t66 = cos(qJ(4));
t82 = t62 * t64 + t63 * t66;
t40 = t82 * qJD(1);
t96 = qJD(1) * t63;
t97 = qJD(1) * t62;
t42 = -t64 * t96 + t66 * t97;
t111 = t42 * t40;
t118 = qJDD(4) - t111;
t120 = t118 * t64;
t119 = t118 * t66;
t115 = 2 * qJD(2);
t68 = qJD(1) ^ 2;
t112 = cos(qJ(1));
t65 = sin(qJ(1));
t78 = t112 * g(1) + t65 * g(2);
t44 = -t68 * pkin(1) + qJDD(1) * qJ(2) - t78;
t86 = qJD(1) * t115 + t44;
t69 = t62 ^ 2;
t71 = t63 ^ 2;
t117 = t69 + t71;
t110 = t63 * t68;
t114 = t63 * pkin(2);
t99 = t62 * qJ(3);
t83 = -t99 - t114;
t47 = t83 * qJD(1);
t80 = t44 + (t115 + t47) * qJD(1);
t113 = t63 * g(3);
t87 = qJDD(3) + t113;
t116 = (-pkin(3) * t110 - pkin(5) * qJDD(1) + t80) * t62 + t87;
t14 = t82 * qJDD(1);
t49 = t117 * t68;
t37 = t40 ^ 2;
t38 = t42 ^ 2;
t104 = t71 * t68;
t105 = t69 * t68;
t89 = t65 * g(1) - t112 * g(2);
t98 = t68 * qJ(2);
t76 = -qJDD(2) + t89 + t98;
t74 = 0.2e1 * qJD(3) * t97 + t76;
t77 = t99 + t63 * (pkin(2) + pkin(3)) + pkin(1);
t13 = t77 * qJDD(1) + t74 + (-t104 - t105) * pkin(5);
t109 = t64 * t13;
t20 = qJDD(4) + t111;
t108 = t64 * t20;
t107 = t66 * t13;
t106 = t66 * t20;
t103 = t86 * t63;
t58 = t69 * qJDD(1);
t60 = t71 * qJDD(1);
t102 = qJ(2) * (t60 + t58) + pkin(1) * t49;
t92 = t63 * qJDD(1);
t101 = -t117 * t63 * t98 + pkin(1) * t92;
t100 = qJ(2) * t62 * t49;
t95 = qJDD(1) * pkin(1);
t94 = t40 * qJD(4);
t93 = t42 * qJD(4);
t59 = t62 * qJDD(1);
t90 = t62 * t92;
t85 = -t62 * g(3) + t103;
t88 = t62 * (t86 * t62 + t113) + t63 * t85;
t35 = t47 * t96;
t81 = t35 + t85;
t12 = -pkin(3) * t104 - pkin(5) * t92 + t81;
t3 = -t66 * t116 + t64 * t12;
t4 = t116 * t64 + t66 * t12;
t1 = -t66 * t3 + t64 * t4;
t2 = t64 * t3 + t66 * t4;
t39 = t66 * t59 - t64 * t92;
t79 = pkin(1) - t83;
t75 = t79 * qJDD(1);
t67 = qJD(4) ^ 2;
t36 = t76 + t95;
t31 = -t38 - t67;
t30 = -t38 + t67;
t29 = t37 - t67;
t25 = t39 - t94;
t24 = t39 - 0.2e1 * t94;
t23 = -t14 - t93;
t22 = 0.2e1 * t93 + t14;
t18 = -t67 - t37;
t17 = t75 + t74;
t16 = -t37 - t38;
t15 = t80 * t62 + t87;
t11 = -t64 * t31 - t106;
t10 = t66 * t31 - t108;
t8 = -t66 * t14 + t64 * t39;
t7 = -t64 * t14 - t66 * t39;
t6 = t66 * t18 - t120;
t5 = t64 * t18 + t119;
t9 = [0, 0, 0, 0, 0, qJDD(1), t89, t78, 0, 0, t58, 0.2e1 * t90, 0, t60, 0, 0, t63 * t36 + t101, t100 + (-t36 - t95) * t62, t88 + t102, pkin(1) * t36 + qJ(2) * t88, t58, 0, -0.2e1 * t90, 0, 0, t60, ((pkin(1) + 0.2e1 * t99 + 0.2e1 * t114) * qJDD(1) + t74) * t63 + t101, t63 * (pkin(2) * t49 + t103 + t35) + (qJ(3) * t49 + qJDD(3) + (qJD(1) * t47 + t86) * t62) * t62 + t102, -t100 + (t74 + 0.2e1 * t75) * t62, qJ(2) * (t62 * t15 + t63 * t81) + t79 * t17, t62 * (t66 * t25 - t64 * t93) + t63 * (-t64 * t25 - t66 * t93), t62 * (-t66 * t22 - t64 * t24) + t63 * (t64 * t22 - t66 * t24), t62 * (-t64 * t30 + t119) + t63 * (-t66 * t30 - t120), t62 * (-t64 * t23 + t66 * t94) + t63 * (-t66 * t23 - t64 * t94), t62 * (t66 * t29 - t108) + t63 * (-t64 * t29 - t106), (t62 * (-t40 * t66 + t42 * t64) + t63 * (t40 * t64 + t42 * t66)) * qJD(4), t62 * (-pkin(5) * t5 + t109) + t63 * (-pkin(5) * t6 + t107) + qJ(2) * (t62 * t5 + t63 * t6) + t77 * t22, t62 * (-pkin(5) * t10 + t107) + t63 * (-pkin(5) * t11 - t109) + qJ(2) * (t62 * t10 + t63 * t11) + t77 * t24, t62 * (-pkin(5) * t7 - t1) + t63 * (-pkin(5) * t8 - t2) + qJ(2) * (t62 * t7 + t63 * t8) + t77 * t16, t77 * t13 + (-pkin(5) + qJ(2)) * (t62 * t1 + t63 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t59, -t49, -t36, 0, 0, 0, 0, 0, 0, -t92, -t49, -t59, -t17, 0, 0, 0, 0, 0, 0, -t22, -t24, -t16, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t110, t59, -t105, t15, 0, 0, 0, 0, 0, 0, t5, t10, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t38 - t37, t39, -t111, -t14, qJDD(4), -t3, -t4, 0, 0;];
tauJ_reg = t9;
