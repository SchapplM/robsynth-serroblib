% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:18
% EndTime: 2019-12-31 16:49:21
% DurationCPUTime: 0.71s
% Computational Cost: add. (1638->166), mult. (3310->230), div. (0->0), fcn. (2075->8), ass. (0->101)
t88 = sin(qJ(3));
t106 = qJD(1) * t88;
t87 = sin(qJ(4));
t89 = cos(qJ(4));
t90 = cos(qJ(3));
t52 = -t89 * t90 * qJD(1) + t87 * t106;
t54 = (t90 * t87 + t88 * t89) * qJD(1);
t39 = t54 * t52;
t78 = qJDD(3) + qJDD(4);
t121 = -t39 + t78;
t124 = t121 * t87;
t123 = t121 * t89;
t103 = qJD(1) * qJD(3);
t101 = t90 * t103;
t105 = t88 * qJDD(1);
t59 = t101 + t105;
t102 = t88 * t103;
t74 = t90 * qJDD(1);
t60 = t74 - t102;
t30 = -t52 * qJD(4) + t89 * t59 + t87 * t60;
t79 = qJD(3) + qJD(4);
t46 = t79 * t52;
t122 = t30 - t46;
t92 = qJD(1) ^ 2;
t117 = sin(qJ(1));
t118 = cos(qJ(1));
t97 = t118 * g(1) + t117 * g(2);
t57 = -t92 * pkin(1) - t97;
t84 = sin(pkin(7));
t85 = cos(pkin(7));
t96 = t117 * g(1) - t118 * g(2);
t94 = qJDD(1) * pkin(1) + t96;
t107 = t85 * t57 + t84 * t94;
t33 = -t92 * pkin(2) + qJDD(1) * pkin(5) + t107;
t82 = -g(3) + qJDD(2);
t27 = t88 * t33 - t90 * t82;
t120 = -t27 + (-t59 + t101) * pkin(6);
t47 = t52 ^ 2;
t48 = t54 ^ 2;
t77 = t79 ^ 2;
t28 = t90 * t33 + t88 * t82;
t66 = qJD(3) * pkin(3) - pkin(6) * t106;
t81 = t90 ^ 2;
t76 = t81 * t92;
t14 = -pkin(3) * t76 + t60 * pkin(6) - qJD(3) * t66 + t28;
t70 = t90 * t92 * t88;
t104 = qJDD(3) + t70;
t93 = t104 * pkin(3) + t120;
t6 = t87 * t14 - t89 * t93;
t111 = t89 * t14;
t7 = t87 * t93 + t111;
t2 = -t89 * t6 + t87 * t7;
t119 = t88 * t2;
t116 = t79 * t87;
t115 = t79 * t89;
t100 = -t84 * t57 + t85 * t94;
t32 = -qJDD(1) * pkin(2) - t92 * pkin(5) - t100;
t17 = -t60 * pkin(3) - pkin(6) * t76 + t66 * t106 + t32;
t114 = t87 * t17;
t36 = t39 + t78;
t113 = t87 * t36;
t112 = t88 * t104;
t110 = t89 * t17;
t109 = t89 * t36;
t65 = qJDD(3) - t70;
t108 = t90 * t65;
t3 = t87 * t6 + t89 * t7;
t12 = t88 * t27 + t90 * t28;
t99 = t87 * t59 - t89 * t60;
t95 = (-qJD(4) + t79) * t54 - t99;
t91 = qJD(3) ^ 2;
t80 = t88 ^ 2;
t75 = t80 * t92;
t68 = -t76 - t91;
t67 = -t75 - t91;
t63 = t75 + t76;
t62 = (t80 + t81) * qJDD(1);
t61 = t74 - 0.2e1 * t102;
t58 = 0.2e1 * t101 + t105;
t44 = -t48 + t77;
t43 = t47 - t77;
t42 = -t48 - t77;
t41 = -t88 * t67 - t108;
t40 = t90 * t68 - t112;
t38 = t48 - t47;
t34 = -t77 - t47;
t31 = -t47 - t48;
t29 = -t54 * qJD(4) - t99;
t25 = -t87 * t42 - t109;
t24 = t89 * t42 - t113;
t23 = t30 + t46;
t18 = (qJD(4) + t79) * t54 + t99;
t16 = t89 * t34 - t124;
t15 = t87 * t34 + t123;
t11 = -t88 * t24 + t90 * t25;
t10 = t87 * t23 + t89 * t95;
t9 = -t89 * t23 + t87 * t95;
t8 = -t88 * t15 + t90 * t16;
t4 = t90 * t10 - t88 * t9;
t1 = t90 * t3 - t119;
t5 = [0, 0, 0, 0, 0, qJDD(1), t96, t97, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t85 * qJDD(1) - t84 * t92) + t100, pkin(1) * (-t84 * qJDD(1) - t85 * t92) - t107, 0, pkin(1) * (t85 * t100 + t84 * t107), (t59 + t101) * t88, t90 * t58 + t88 * t61, t112 + t90 * (-t75 + t91), (t60 - t102) * t90, t88 * (t76 - t91) + t108, 0, -t90 * t32 + pkin(2) * t61 + pkin(5) * t40 + pkin(1) * (t84 * t40 + t85 * t61), t88 * t32 - pkin(2) * t58 + pkin(5) * t41 + pkin(1) * (t84 * t41 - t85 * t58), pkin(2) * t63 + pkin(5) * t62 + pkin(1) * (t84 * t62 + t85 * t63) + t12, -pkin(2) * t32 + pkin(5) * t12 + pkin(1) * (t84 * t12 - t85 * t32), t88 * (-t116 * t54 + t89 * t30) + t90 * (t54 * t115 + t87 * t30), t88 * (-t122 * t87 - t89 * t18) + t90 * (t122 * t89 - t87 * t18), t88 * (-t87 * t44 + t123) + t90 * (t89 * t44 + t124), t88 * (t52 * t115 - t87 * t29) + t90 * (t52 * t116 + t89 * t29), t88 * (t89 * t43 - t113) + t90 * (t87 * t43 + t109), (t88 * (-t52 * t89 + t54 * t87) + t90 * (-t52 * t87 - t54 * t89)) * t79, t88 * (-pkin(6) * t15 + t114) + t90 * (-pkin(3) * t18 + pkin(6) * t16 - t110) - pkin(2) * t18 + pkin(5) * t8 + pkin(1) * (-t85 * t18 + t84 * t8), t88 * (-pkin(6) * t24 + t110) + t90 * (-pkin(3) * t122 + pkin(6) * t25 + t114) - pkin(2) * t122 + pkin(5) * t11 + pkin(1) * (t84 * t11 - t122 * t85), t88 * (-pkin(6) * t9 - t2) + t90 * (-pkin(3) * t31 + pkin(6) * t10 + t3) - pkin(2) * t31 + pkin(5) * t4 + pkin(1) * (-t85 * t31 + t84 * t4), -pkin(6) * t119 + t90 * (-pkin(3) * t17 + pkin(6) * t3) - pkin(2) * t17 + pkin(5) * t1 + pkin(1) * (t84 * t1 - t85 * t17); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, t104 * t90 + t88 * t68, -t88 * t65 + t90 * t67, 0, -t90 * t27 + t88 * t28, 0, 0, 0, 0, 0, 0, t90 * t15 + t88 * t16, t90 * t24 + t88 * t25, t88 * t10 + t90 * t9, t90 * t2 + t88 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t75 - t76, t105, t70, t74, qJDD(3), -t27, -t28, 0, 0, t39, t38, t23, -t39, t95, t78, pkin(3) * t15 - t6, -t111 - t87 * t120 + (-t104 * t87 + t24) * pkin(3), pkin(3) * t9, pkin(3) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t38, t23, -t39, t95, t78, -t6, -t7, 0, 0;];
tauJ_reg = t5;
