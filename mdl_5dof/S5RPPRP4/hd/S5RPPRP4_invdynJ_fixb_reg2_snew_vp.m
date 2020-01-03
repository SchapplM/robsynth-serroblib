% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:20
% EndTime: 2019-12-31 17:52:22
% DurationCPUTime: 0.62s
% Computational Cost: add. (1311->137), mult. (2453->162), div. (0->0), fcn. (1016->6), ass. (0->86)
t76 = sin(qJ(4));
t78 = cos(qJ(4));
t81 = qJD(1) ^ 2;
t57 = t78 * t81 * t76;
t52 = qJDD(4) + t57;
t113 = pkin(4) * t52;
t112 = pkin(1) + pkin(2);
t68 = qJDD(1) * qJ(2);
t77 = sin(qJ(1));
t79 = cos(qJ(1));
t93 = t79 * g(1) + t77 * g(2);
t90 = (2 * qJD(2) * qJD(1)) - t93;
t88 = t68 + t90;
t37 = -t112 * t81 + t88;
t73 = sin(pkin(7));
t74 = cos(pkin(7));
t98 = t77 * g(1) - t79 * g(2);
t92 = qJDD(2) - t98;
t86 = -t81 * qJ(2) + t92;
t82 = -t112 * qJDD(1) + t86;
t15 = t74 * t37 + t73 * t82;
t13 = -t81 * pkin(3) - qJDD(1) * pkin(6) + t15;
t71 = g(3) + qJDD(3);
t10 = t78 * t13 + t76 * t71;
t100 = qJD(1) * qJD(5);
t102 = t78 * qJDD(1);
t101 = qJD(1) * qJD(4);
t96 = t76 * t101;
t45 = t96 - t102;
t105 = qJD(1) * t76;
t51 = qJD(4) * pkin(4) + qJ(5) * t105;
t87 = t45 * qJ(5) - qJD(4) * t51 - 0.2e1 * t78 * t100 + t10;
t69 = t76 ^ 2;
t111 = t69 * t81;
t70 = t78 ^ 2;
t110 = t70 * t81;
t109 = t76 * t52;
t53 = qJDD(4) - t57;
t108 = t78 * t53;
t107 = t69 + t70;
t106 = qJ(5) * t76;
t104 = qJDD(1) * pkin(1);
t103 = t76 * qJDD(1);
t99 = 2 * t100;
t95 = t78 * t101;
t63 = t78 * t71;
t9 = t76 * t13 - t63;
t94 = t73 * t37 - t74 * t82;
t4 = t78 * t10 + t76 * t9;
t43 = -t95 - t103;
t91 = -t43 - t95;
t12 = qJDD(1) * pkin(3) - t81 * pkin(6) + t94;
t80 = qJD(4) ^ 2;
t56 = -t80 - t110;
t33 = t78 * t56 - t109;
t44 = -0.2e1 * t96 + t102;
t20 = t73 * t33 - t74 * t44;
t85 = pkin(3) * t44 - pkin(6) * t33 + qJ(2) * (t74 * t33 + t73 * t44) - t112 * t20;
t55 = -t80 - t111;
t34 = -t76 * t55 - t108;
t42 = 0.2e1 * t95 + t103;
t21 = t73 * t34 + t74 * t42;
t84 = -pkin(3) * t42 - pkin(6) * t34 + qJ(2) * (t74 * t34 - t73 * t42) - t112 * t21;
t48 = t107 * qJDD(1);
t49 = t107 * t81;
t25 = -t73 * t48 + t74 * t49;
t83 = -pkin(3) * t49 + pkin(6) * t48 + qJ(2) * (-t74 * t48 - t73 * t49) - t112 * t25;
t8 = -t45 * pkin(4) - qJ(5) * t110 - t51 * t105 + qJDD(5) + t12;
t6 = t91 * qJ(5) + t76 * t99 + t113 - t9;
t50 = (t69 - t70) * t81;
t47 = t74 * qJDD(1) + t73 * t81;
t46 = -t73 * qJDD(1) + t74 * t81;
t38 = -t86 + t104;
t32 = -t76 * t53 + t78 * t55;
t31 = -t109 - t78 * (t80 - t111);
t30 = t78 * t52 + t76 * t56;
t29 = -t76 * (-t80 + t110) - t108;
t28 = (-t43 + t95) * t76;
t27 = (-t45 - t96) * t78;
t22 = t78 * t42 + t76 * t44;
t7 = -pkin(4) * t110 + t87;
t5 = t73 * t15 - t74 * t94;
t3 = -t74 * t12 + t73 * t4;
t2 = -t76 * t6 + t78 * t7;
t1 = t73 * t2 - t74 * t8;
t11 = [0, 0, 0, 0, 0, qJDD(1), t98, t93, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -t92 + 0.2e1 * t104, 0, 0.2e1 * t68 + t90, qJ(2) * (-t81 * pkin(1) + t88) + pkin(1) * t38, 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t46 + t112 * t47 + t94, qJ(2) * t47 + t112 * t46 + t15, 0, qJ(2) * (t74 * t15 + t73 * t94) - t112 * t5, t28, t22, t31, t27, t29, 0, t78 * t12 + t85, -t76 * t12 + t84, -t4 + t83, qJ(2) * (t73 * t12 + t74 * t4) + pkin(3) * t12 - pkin(6) * t4 - t112 * t3, t28, t22, t31, t27, t29, 0, t52 * t106 - t78 * (-pkin(4) * t44 + qJ(5) * t56 - t8) + t85, -t76 * (-qJ(5) * t55 + t8) - t78 * (pkin(4) * t42 - qJ(5) * t53) + t84, -t78 * (-qJ(5) * t102 + (t49 - t110) * pkin(4) + t87) + (t63 + (-t13 + t99) * t76 + (t91 + t103) * qJ(5) + t113) * t76 + t83, qJ(2) * (t74 * t2 + t73 * t8) + t6 * t106 - t78 * (-pkin(4) * t8 + qJ(5) * t7) + pkin(3) * t8 - pkin(6) * t2 - t112 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t81, -t38, 0, 0, 0, 0, 0, 0, -t47, -t46, 0, t5, 0, 0, 0, 0, 0, 0, t20, t21, t25, t3, 0, 0, 0, 0, 0, 0, t20, t21, t25, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, 0, 0, 0, 0, t30, t32, 0, t76 * t10 - t78 * t9, 0, 0, 0, 0, 0, 0, t30, t32, 0, t78 * t6 + t76 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t50, -t103, t57, -t102, qJDD(4), -t9, -t10, 0, 0, -t57, t50, -t103, t57, -t102, qJDD(4), t6 + t113, (t55 + t110) * pkin(4) - t87, pkin(4) * t103, pkin(4) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t42, -t49, t8;];
tauJ_reg = t11;
