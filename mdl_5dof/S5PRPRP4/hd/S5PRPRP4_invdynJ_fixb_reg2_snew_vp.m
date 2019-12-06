% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRP4
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:06
% EndTime: 2019-12-05 15:36:09
% DurationCPUTime: 0.67s
% Computational Cost: add. (1048->126), mult. (1936->168), div. (0->0), fcn. (1161->8), ass. (0->83)
t80 = sin(qJ(4));
t82 = cos(qJ(4));
t85 = qJD(2) ^ 2;
t63 = t80 * t85 * t82;
t58 = qJDD(4) - t63;
t110 = t82 * t58;
t70 = t80 ^ 2;
t117 = t70 * t85;
t84 = qJD(4) ^ 2;
t59 = t84 + t117;
t34 = -t80 * t59 + t110;
t103 = qJD(2) * qJD(4);
t105 = t80 * qJDD(2);
t49 = 0.2e1 * t82 * t103 + t105;
t75 = sin(pkin(8));
t77 = cos(pkin(8));
t21 = t75 * t34 + t77 * t49;
t124 = pkin(2) * t21 + pkin(6) * t34;
t81 = sin(qJ(2));
t83 = cos(qJ(2));
t123 = t83 * t21 + t81 * (t77 * t34 - t75 * t49);
t71 = t82 ^ 2;
t116 = t71 * t85;
t119 = t110 + t80 * (-t84 + t116);
t108 = t80 * qJ(5);
t93 = -t82 * pkin(4) - t108;
t107 = t85 * t93;
t72 = -g(3) + qJDD(1);
t76 = sin(pkin(7));
t78 = cos(pkin(7));
t94 = -t78 * g(1) - t76 * g(2);
t38 = t81 * t72 + t83 * t94;
t36 = -t85 * pkin(2) + t38;
t37 = t83 * t72 - t81 * t94;
t86 = qJDD(2) * pkin(2) + t37;
t18 = t77 * t36 + t75 * t86;
t16 = -t85 * pkin(3) + qJDD(2) * pkin(6) + t18;
t89 = -t76 * g(1) + t78 * g(2) + qJDD(3);
t41 = t82 * t89;
t8 = -qJDD(4) * pkin(4) - t84 * qJ(5) + (t16 + t107) * t80 + qJDD(5) - t41;
t118 = 2 * qJD(5);
t57 = qJDD(4) + t63;
t113 = t80 * t57;
t12 = t82 * t16 + t80 * t89;
t109 = t70 + t71;
t106 = qJD(2) * t80;
t104 = t82 * qJDD(2);
t61 = -t84 - t116;
t33 = t82 * t61 - t113;
t100 = t80 * t103;
t50 = -0.2e1 * t100 + t104;
t20 = t75 * t33 + t77 * t50;
t102 = pkin(2) * t20 + pkin(3) * t50 + pkin(6) * t33;
t54 = t109 * qJDD(2);
t55 = t109 * t85;
t24 = t75 * t54 + t77 * t55;
t101 = pkin(2) * t24 + pkin(3) * t55 + pkin(6) * t54;
t11 = t80 * t16 - t41;
t4 = t80 * t11 + t82 * t12;
t98 = t75 * t36 - t77 * t86;
t95 = qJDD(4) * qJ(5) + (qJD(4) * t118) + t82 * t107 + t12;
t92 = t82 * t49 + t80 * t50;
t91 = t80 * t58 + t82 * t59;
t15 = -qJDD(2) * pkin(3) - t85 * pkin(6) + t98;
t90 = pkin(3) - t93;
t88 = t15 - (-t100 + t104) * pkin(4) - qJ(5) * t49;
t87 = t106 * t118 - t88;
t56 = (t70 - t71) * t85;
t53 = -t75 * qJDD(2) - t77 * t85;
t52 = t77 * qJDD(2) - t75 * t85;
t32 = t113 + t82 * (t84 - t117);
t31 = t82 * t57 + t80 * t61;
t30 = t49 * t80;
t29 = t50 * t82;
t13 = t81 * (t77 * t54 - t75 * t55) + t83 * t24;
t9 = t81 * (t77 * t33 - t75 * t50) + t83 * t20;
t7 = -t84 * pkin(4) + t95;
t6 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t106 + t88;
t5 = t75 * t18 - t77 * t98;
t3 = t82 * t7 + t80 * t8;
t2 = -t77 * t15 + t75 * t4;
t1 = t75 * t3 - t77 * t6;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, 0, 0, 0, 0, t83 * qJDD(2) - t81 * t85, -t81 * qJDD(2) - t83 * t85, 0, t83 * t37 + t81 * t38, 0, 0, 0, 0, 0, 0, t83 * t52 + t81 * t53, -t81 * t52 + t83 * t53, 0, t81 * (t77 * t18 + t75 * t98) + t83 * t5, 0, 0, 0, 0, 0, 0, t9, -t123, t13, t81 * (t75 * t15 + t77 * t4) + t83 * t2, 0, 0, 0, 0, 0, 0, t9, t13, t123, t81 * (t77 * t3 + t75 * t6) + t83 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t37, -t38, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t52 - t98, pkin(2) * t53 - t18, 0, pkin(2) * t5, t30, t92, t32, t29, t119, 0, -t82 * t15 + t102, -pkin(3) * t49 + t80 * t15 - t124, t4 + t101, pkin(2) * t2 - pkin(3) * t15 + pkin(6) * t4, t30, t32, -t92, 0, -t119, t29, t50 * t108 + t82 * ((t50 - t100) * pkin(4) + t87) + t102, t82 * ((t55 - t84) * pkin(4) + t95) + (qJ(5) * t55 + t8) * t80 + t101, t80 * (-pkin(4) * t100 + t87) + t90 * t49 + t124, pkin(2) * t1 + pkin(6) * t3 - t90 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, 0, 0, 0, 0, t31, -t91, 0, -t82 * t11 + t80 * t12, 0, 0, 0, 0, 0, 0, t31, 0, t91, t80 * t7 - t82 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t56, t105, t63, t104, qJDD(4), -t11, -t12, 0, 0, -t63, t105, -t56, qJDD(4), -t104, t63, pkin(4) * t57 + qJ(5) * t61 - t8, (-pkin(4) * t80 + qJ(5) * t82) * qJDD(2), qJ(5) * t58 + (t59 - t84) * pkin(4) + t95, -pkin(4) * t8 + qJ(5) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t105, -t59, t8;];
tauJ_reg = t10;
