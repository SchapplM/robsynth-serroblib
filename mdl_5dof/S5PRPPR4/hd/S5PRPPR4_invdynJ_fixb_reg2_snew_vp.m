% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:55
% EndTime: 2019-12-31 17:36:58
% DurationCPUTime: 0.76s
% Computational Cost: add. (1196->141), mult. (2731->192), div. (0->0), fcn. (1774->8), ass. (0->96)
t73 = sin(pkin(8));
t74 = cos(pkin(8));
t75 = sin(qJ(5));
t77 = cos(qJ(5));
t88 = t73 * t75 + t74 * t77;
t43 = t88 * qJD(2);
t102 = qJD(2) * t74;
t103 = qJD(2) * t73;
t45 = -t102 * t75 + t103 * t77;
t115 = t45 * t43;
t124 = qJDD(5) - t115;
t126 = t124 * t75;
t125 = t124 * t77;
t119 = 2 * qJD(3);
t79 = qJD(2) ^ 2;
t116 = cos(qJ(2));
t104 = sin(pkin(7));
t105 = cos(pkin(7));
t53 = -g(1) * t105 - g(2) * t104;
t76 = sin(qJ(2));
t84 = g(1) * t104 - g(2) * t105;
t82 = -t116 * t53 - t76 * t84;
t32 = -t79 * pkin(2) + qJDD(2) * qJ(3) - t82;
t91 = qJD(2) * t119 + t32;
t68 = t73 ^ 2;
t69 = t74 ^ 2;
t123 = t68 + t69;
t114 = t74 * t79;
t106 = t73 * qJ(4);
t117 = t74 * pkin(3);
t89 = -t106 - t117;
t50 = t89 * qJD(2);
t87 = t32 + (t119 + t50) * qJD(2);
t71 = -g(3) + qJDD(1);
t60 = t74 * t71;
t98 = qJDD(4) - t60;
t122 = (-pkin(4) * t114 - pkin(6) * qJDD(2) + t87) * t73 + t98;
t66 = t73 * qJDD(2);
t61 = qJ(4) * t66;
t95 = qJD(4) * t103;
t121 = 0.2e1 * t61 + 0.2e1 * t95;
t18 = t88 * qJDD(2);
t52 = t123 * t79;
t86 = t74 * (pkin(3) + pkin(4)) + pkin(2) + t106;
t40 = t43 ^ 2;
t41 = t45 ^ 2;
t70 = t79 * qJ(3);
t72 = qJDD(2) * pkin(2);
t92 = t116 * t84 - t76 * t53;
t27 = qJDD(3) - t70 - t72 - t92;
t99 = t74 * qJDD(2);
t63 = pkin(3) * t99;
t85 = -t27 + t63;
t17 = -t61 - t85 - 0.2e1 * t95;
t11 = -pkin(4) * t99 + pkin(6) * t52 + t17;
t113 = t75 * t11;
t25 = qJDD(5) + t115;
t112 = t75 * t25;
t111 = t77 * t11;
t110 = t77 * t25;
t65 = t68 * qJDD(2);
t67 = t69 * qJDD(2);
t109 = qJ(3) * (t67 + t65) + pkin(2) * t52;
t108 = -t123 * t74 * t70 + pkin(2) * t99;
t107 = qJ(3) * t73 * t52;
t101 = t43 * qJD(5);
t100 = t45 * qJD(5);
t20 = t73 * t71 + t91 * t74;
t96 = t73 * t99;
t19 = t73 * t91 - t60;
t93 = t73 * t19 + t74 * t20;
t13 = t50 * t102 + t20;
t6 = -t69 * t79 * pkin(4) - pkin(6) * t99 + t13;
t3 = -t77 * t122 + t75 * t6;
t4 = t122 * t75 + t77 * t6;
t1 = -t77 * t3 + t75 * t4;
t2 = t75 * t3 + t77 * t4;
t42 = t66 * t77 - t75 * t99;
t78 = qJD(5) ^ 2;
t36 = -t41 - t78;
t35 = -t41 + t78;
t34 = t40 - t78;
t31 = t42 - t101;
t30 = t42 - 0.2e1 * t101;
t29 = -t18 - t100;
t28 = 0.2e1 * t100 + t18;
t23 = -t78 - t40;
t21 = -t40 - t41;
t15 = -t75 * t36 - t110;
t14 = t77 * t36 - t112;
t12 = t73 * t87 + t98;
t10 = -t77 * t18 + t75 * t42;
t9 = -t75 * t18 - t77 * t42;
t8 = t77 * t23 - t126;
t7 = t75 * t23 + t125;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74 * t19 + t73 * t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74 * t12 + t73 * t13, 0, 0, 0, 0, 0, 0, -t74 * t7 + t73 * t8, -t74 * t14 + t73 * t15, t73 * t10 - t74 * t9, -t74 * t1 + t73 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t92, t82, 0, 0, t65, 0.2e1 * t96, 0, t67, 0, 0, -t74 * t27 + t108, t107 + (t27 - t72) * t73, t93 + t109, -pkin(2) * t27 + qJ(3) * t93, t65, 0, -0.2e1 * t96, 0, 0, t67, (-t27 + 0.2e1 * t63 + t121) * t74 + t108, t74 * (pkin(3) * t52 + t13) + (qJ(4) * t52 + (qJD(2) * t50 + t91) * t73 + t98) * t73 + t109, -t107 + ((pkin(2) + t117) * qJDD(2) + t85 + t121) * t73, qJ(3) * (t73 * t12 + t74 * t13) + (-pkin(2) + t89) * t17, t73 * (-t100 * t75 + t77 * t31) + t74 * (-t100 * t77 - t75 * t31), t73 * (-t77 * t28 - t75 * t30) + t74 * (t75 * t28 - t77 * t30), t73 * (-t75 * t35 + t125) + t74 * (-t77 * t35 - t126), t73 * (t101 * t77 - t75 * t29) + t74 * (-t101 * t75 - t77 * t29), t73 * (t77 * t34 - t112) + t74 * (-t75 * t34 - t110), (t73 * (-t43 * t77 + t45 * t75) + t74 * (t43 * t75 + t45 * t77)) * qJD(5), t73 * (-pkin(6) * t7 - t113) + t74 * (-pkin(6) * t8 - t111) + qJ(3) * (t73 * t7 + t74 * t8) + t86 * t28, t73 * (-pkin(6) * t14 - t111) + t74 * (-pkin(6) * t15 + t113) + qJ(3) * (t73 * t14 + t74 * t15) + t86 * t30, t73 * (-pkin(6) * t9 - t1) + t74 * (-pkin(6) * t10 - t2) + qJ(3) * (t74 * t10 + t73 * t9) + t86 * t21, -t86 * t11 + (qJ(3) - pkin(6)) * (t73 * t1 + t74 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t66, -t52, t27, 0, 0, 0, 0, 0, 0, -t99, -t52, -t66, t17, 0, 0, 0, 0, 0, 0, -t28, -t30, -t21, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 * t114, t66, -t68 * t79, t12, 0, 0, 0, 0, 0, 0, t7, t14, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t41 - t40, t42, -t115, -t18, qJDD(5), -t3, -t4, 0, 0;];
tauJ_reg = t5;
