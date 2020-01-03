% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:12
% EndTime: 2019-12-31 17:41:14
% DurationCPUTime: 0.60s
% Computational Cost: add. (855->158), mult. (1852->157), div. (0->0), fcn. (928->6), ass. (0->99)
t76 = sin(qJ(3));
t78 = cos(qJ(3));
t80 = qJD(2) ^ 2;
t107 = t76 * t80 * t78;
t47 = qJDD(3) + t107;
t141 = pkin(4) * t47;
t48 = qJDD(3) - t107;
t124 = t78 * t48;
t110 = qJD(2) * qJD(3);
t103 = t78 * t110;
t113 = t76 * qJDD(2);
t35 = 0.2e1 * t103 + t113;
t68 = t76 ^ 2;
t131 = t68 * t80;
t79 = qJD(3) ^ 2;
t52 = t79 + t131;
t122 = pkin(2) * t35 + pkin(6) * (-t76 * t52 + t124);
t112 = t78 * qJDD(2);
t64 = t76 * t110;
t36 = -t64 + t112;
t115 = qJD(2) * t76;
t46 = -qJD(3) * pkin(4) - qJ(5) * t115;
t140 = t36 * qJ(5) - qJD(3) * t46;
t69 = t78 ^ 2;
t130 = t69 * t80;
t55 = -t79 - t130;
t139 = pkin(3) * t47 + qJ(4) * t55;
t138 = qJDD(2) * pkin(6);
t18 = t124 + t76 * (-t79 + t130);
t136 = t46 * t115 + qJDD(5);
t135 = t36 * pkin(4) + t136;
t114 = -g(3) + qJDD(1);
t63 = t78 * t114;
t98 = qJDD(3) * pkin(3) + t79 * qJ(4) - qJDD(4) + t63;
t117 = t76 * qJ(4);
t94 = -t78 * pkin(3) - t117;
t33 = t94 * qJD(2);
t116 = qJD(2) * t33;
t132 = cos(qJ(2));
t73 = sin(pkin(7));
t74 = cos(pkin(7));
t45 = -t74 * g(1) - t73 * g(2);
t77 = sin(qJ(2));
t96 = t73 * g(1) - t74 * g(2);
t85 = -t132 * t45 - t77 * t96;
t12 = -t80 * pkin(2) + t138 - t85;
t99 = t12 + t116;
t6 = t99 * t76 - t98;
t134 = pkin(3) + pkin(4);
t120 = t68 + t69;
t34 = t120 * t138;
t129 = t76 * t12;
t128 = t76 * t47;
t125 = t78 * t35;
t9 = t76 * t114 + t78 * t12;
t37 = -0.2e1 * t64 + t112;
t123 = pkin(6) * (t78 * t55 - t128) + pkin(2) * t37;
t43 = t120 * t80;
t121 = pkin(2) * t43 + t34;
t119 = qJ(4) * t43;
t118 = qJ(4) * t78;
t111 = qJ(5) * qJDD(2);
t109 = qJD(2) * qJD(5);
t108 = qJD(4) * qJD(3);
t106 = 0.2e1 * t109;
t8 = -t63 + t129;
t105 = t76 * t8 + t78 * t9;
t104 = qJ(5) * qJD(3) * t78;
t102 = t132 * t96 - t77 * t45;
t101 = t52 - t130;
t100 = qJDD(3) * qJ(4) + t78 * t116 + t9;
t97 = -t78 * t134 - pkin(2);
t93 = t76 * t37 + t125;
t20 = -t78 * (-t79 + t131) + t128;
t21 = t76 * t48 + t78 * t52;
t11 = -qJDD(2) * pkin(2) - t80 * pkin(6) - t102;
t92 = t79 * pkin(3) - t100;
t90 = t103 + t113;
t91 = -t90 * qJ(5) - t141 - t98;
t65 = 0.2e1 * t108;
t5 = t65 - t92;
t89 = -t36 * pkin(3) + t11 + (-t103 - t90) * qJ(4);
t88 = pkin(3) * t52 + qJ(4) * t48 + t5;
t87 = 0.2e1 * qJD(4) * t115 - t89;
t86 = t91 + t129;
t84 = pkin(4) * t130 + t140 + t92;
t83 = (pkin(3) * qJD(3) - 0.2e1 * qJD(4)) * t115 + t89;
t82 = (t37 - t64) * pkin(3) + t87;
t81 = -pkin(3) * t64 + qJ(4) * t35 + t87;
t60 = -0.2e1 * t78 * t109;
t59 = t76 * t106;
t44 = (t68 - t69) * t80;
t19 = t78 * t47 + t76 * t55;
t17 = t35 * t76;
t16 = (t36 - t64) * t78;
t3 = (t104 + (-0.2e1 * qJD(5) + t33) * t76) * qJD(2) + t86;
t2 = t60 + t65 - t84;
t1 = qJ(5) * t130 - t135 + t83;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, 0, 0, 0, 0, 0, t19, -t21, 0, t76 * t9 - t78 * t8, 0, 0, 0, 0, 0, 0, t19, 0, t21, t76 * t5 - t78 * t6, 0, 0, 0, 0, 0, 0, t19, t21, 0, t76 * t2 - t78 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t102, t85, 0, 0, t17, t93, t20, t16, t18, 0, -t78 * t11 + t123, t76 * t11 - t122, t105 + t121, -pkin(2) * t11 + pkin(6) * t105, t17, t20, -t93, 0, -t18, t16, t37 * t117 + t78 * t82 + t123, t78 * (t65 + (t43 - t79) * pkin(3) + t100) + (t6 + t119) * t76 + t121, pkin(3) * t125 + t76 * t81 + t122, pkin(6) * (t78 * t5 + t76 * t6) + (-pkin(2) + t94) * t83, t17, -t93, -t20, t16, t18, 0, t76 * (qJ(4) * t37 + qJ(5) * t47) + t78 * ((-t55 - t130) * qJ(5) + (t36 + t37) * pkin(4) + t82 + t136) + t123, t76 * (t101 * qJ(5) + t135 + t81) + t78 * (-qJ(5) * t48 + t134 * t35) + t122, -t34 + (-0.2e1 * t108 + (t106 + t111) * t78 + t84) * t78 + t97 * t43 + (-qJ(5) * t103 - t119 + t59 + (-t99 + t111) * t76 - t91) * t76, (t97 - t117) * t1 + (pkin(6) - qJ(5)) * (t78 * t2 + t76 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, t44, t113, t107, t112, qJDD(3), -t8, -t9, 0, 0, -t107, t113, -t44, qJDD(3), -t112, t107, -t6 + t139, (-pkin(3) * t76 + t118) * qJDD(2), t88, -pkin(3) * t6 + qJ(4) * t5, -t107, -t44, -t113, t107, t112, qJDD(3), t141 + t59 + (-t33 * t76 - t104) * qJD(2) - t86 + t139, t101 * pkin(4) - t140 + t60 + t88, (t134 * t76 - t118) * qJDD(2), qJ(4) * t2 - t134 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t113, -t52, t6, 0, 0, 0, 0, 0, 0, -t47, -t52, -t113, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t35, -t43, -t1;];
tauJ_reg = t4;
