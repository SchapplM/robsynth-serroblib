% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:18
% EndTime: 2019-12-31 18:01:20
% DurationCPUTime: 0.99s
% Computational Cost: add. (2078->191), mult. (3124->230), div. (0->0), fcn. (1770->10), ass. (0->121)
t155 = qJD(1) - qJD(4);
t80 = sin(pkin(8));
t81 = cos(pkin(8));
t83 = sin(qJ(4));
t86 = cos(qJ(4));
t43 = t83 * t80 - t86 * t81;
t154 = t155 * t43;
t156 = t154 * t155;
t131 = t80 * qJ(2);
t88 = -pkin(1) - pkin(2);
t50 = t81 * t88 - t131;
t47 = -pkin(3) + t50;
t51 = t81 * qJ(2) + t80 * t88;
t20 = t83 * t47 + t86 * t51;
t123 = pkin(8) + qJ(4);
t113 = sin(t123);
t114 = cos(t123);
t84 = sin(qJ(1));
t87 = cos(qJ(1));
t35 = -t84 * t113 - t87 * t114;
t36 = t87 * t113 - t84 * t114;
t109 = g(1) * t35 + g(2) * t36;
t117 = -pkin(3) - t131;
t61 = t88 * qJD(1) + qJD(2);
t52 = t81 * t61;
t30 = t117 * qJD(1) + t52;
t126 = qJ(2) * qJD(1);
t34 = t81 * t126 + t80 * t61;
t14 = t83 * t30 + t86 * t34;
t10 = -pkin(7) * t155 + t14;
t82 = sin(qJ(5));
t85 = cos(qJ(5));
t7 = t85 * qJD(3) - t82 * t10;
t132 = t7 * qJD(5);
t128 = qJD(4) * t86;
t129 = qJD(4) * t83;
t124 = qJD(1) * qJD(2);
t119 = t80 * t124;
t57 = t88 * qJDD(1) + qJDD(2);
t49 = t81 * t57;
t23 = t117 * qJDD(1) - t119 + t49;
t118 = t81 * t124;
t125 = qJ(2) * qJDD(1);
t137 = t81 * t125 + t80 * t57;
t28 = t118 + t137;
t5 = t30 * t128 - t34 * t129 + t83 * t23 + t86 * t28;
t75 = qJDD(1) - qJDD(4);
t3 = -t75 * pkin(7) + t5;
t1 = t82 * qJDD(3) + t85 * t3 + t132;
t68 = t85 * qJDD(3);
t142 = t85 * t10;
t8 = t82 * qJD(3) + t142;
t2 = -t8 * qJD(5) - t82 * t3 + t68;
t94 = -(t7 * t85 + t8 * t82) * qJD(5) + t1 * t85 - t2 * t82;
t91 = t109 + t94;
t45 = t86 * t80 + t83 * t81;
t139 = t155 * t45;
t101 = -t139 * t155 + t43 * t75;
t89 = qJD(5) ^ 2;
t153 = -t45 * t89 + t101;
t152 = t75 * pkin(4);
t151 = t155 * pkin(4);
t19 = t86 * t47 - t83 * t51;
t11 = -t43 * qJD(2) + t19 * qJD(4);
t150 = t11 * t155;
t12 = t45 * qJD(2) + t20 * qJD(4);
t149 = t12 * t155;
t13 = t86 * t30 - t83 * t34;
t148 = t13 * t155;
t147 = t14 * t155;
t144 = t82 * t85;
t143 = t84 * t80;
t141 = t85 * t75;
t140 = t87 * t80;
t136 = t87 * pkin(1) + t84 * qJ(2);
t135 = g(1) * t84 - g(2) * t87;
t77 = t82 ^ 2;
t78 = t85 ^ 2;
t134 = t77 - t78;
t133 = t77 + t78;
t130 = pkin(1) * qJDD(1);
t127 = qJD(5) * t155;
t79 = qJDD(3) + g(3);
t74 = t155 ^ 2;
t122 = t74 * t144;
t121 = 0.2e1 * t124;
t120 = t133 * t75;
t116 = t34 * t128 + t30 * t129 - t86 * t23 + t83 * t28;
t66 = t81 * pkin(3) + pkin(2);
t115 = pkin(3) * t143 + t87 * t66 + t136;
t112 = qJDD(2) - t130;
t111 = -0.2e1 * t127 * t144;
t110 = g(1) * t36 - g(2) * t35;
t108 = g(1) * t87 + g(2) * t84;
t105 = t7 * t82 - t8 * t85;
t104 = (-t80 * t126 + t52) * t80 - t34 * t81;
t4 = t116 + t152;
t103 = t110 - t4;
t70 = t87 * qJ(2);
t102 = pkin(3) * t140 + t70 + (-pkin(1) - t66) * t84;
t100 = t110 - t116;
t9 = -t13 + t151;
t99 = -pkin(7) * qJDD(5) + (t13 + t9 + t151) * qJD(5);
t15 = pkin(4) - t19;
t16 = -pkin(7) + t20;
t98 = -qJDD(5) * t16 + (-t15 * t155 - t11 - t9) * qJD(5);
t97 = -0.2e1 * t154 * qJD(5) - qJDD(5) * t45;
t96 = -qJD(3) * qJD(5) + t155 * t9 - t109 - t3;
t95 = -t109 - t5;
t93 = pkin(7) * t89 - t103 + t147 + t152;
t92 = -t15 * t75 + t16 * t89 + t103 - t149;
t90 = qJD(1) ^ 2;
t55 = qJDD(5) * t85 - t89 * t82;
t54 = qJDD(5) * t82 + t89 * t85;
t46 = t87 * t81 + t143;
t44 = -t84 * t81 + t140;
t32 = t78 * t75 + t111;
t31 = -t77 * t75 + t111;
t27 = t49 + (-t124 - t125) * t80;
t24 = t134 * t127 - t82 * t141;
t6 = [0, 0, 0, 0, 0, qJDD(1), t135, t108, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t130 + t135, 0, -t108 + t121 + 0.2e1 * t125, -t112 * pkin(1) - g(1) * (-t84 * pkin(1) + t70) - g(2) * t136 + (t121 + t125) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t119 - g(1) * t44 - g(2) * t46 - t49 + (-t50 + t131) * qJDD(1), -g(1) * t46 + g(2) * t44 + t51 * qJDD(1) + 0.2e1 * t118 + t137, 0, t28 * t51 + t27 * t50 - g(1) * (t88 * t84 + t70) - g(2) * (t87 * pkin(2) + t136) - t104 * qJD(2), 0, 0, 0, 0, 0, t75, -t19 * t75 - t100 + t149, t20 * t75 + t150 - t95, 0, -g(1) * t102 - g(2) * t115 + t14 * t11 - t116 * t19 - t13 * t12 + t5 * t20, -t31, -0.2e1 * t24, -t54, t32, -t55, 0, t98 * t82 - t92 * t85, t92 * t82 + t98 * t85, -t16 * t120 - t133 * t150 - t91, t4 * t15 + t9 * t12 - g(1) * (t36 * pkin(4) + t35 * pkin(7) + t102) - g(2) * (-t35 * pkin(4) + t36 * pkin(7) + t115) - t105 * t11 + t94 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t90, -t90 * qJ(2) + t112 - t135, 0, 0, 0, 0, 0, 0, -t81 * qJDD(1) - t80 * t90, t80 * qJDD(1) - t81 * t90, 0, t104 * qJD(1) + t27 * t81 + t28 * t80 - t135, 0, 0, 0, 0, 0, 0, t101, t45 * t75 + t156, 0, t116 * t43 + t139 * t13 + t14 * t154 + t5 * t45 - t135, 0, 0, 0, 0, 0, 0, t153 * t85 + t97 * t82, -t153 * t82 + t97 * t85, -t45 * t120 - t133 * t156, -t154 * t105 - t139 * t9 + t4 * t43 + t94 * t45 - t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, 0, t55, -t54, 0, -t105 * qJD(5) + t1 * t82 + t2 * t85 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t100 - t147, t95 - t148, 0, 0, t31, 0.2e1 * t24, t54, -t32, t55, 0, t99 * t82 - t93 * t85, t93 * t82 + t99 * t85, -pkin(7) * t120 + t133 * t148 + t91, t103 * pkin(4) + t91 * pkin(7) + t105 * t13 - t9 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t134 * t74, -t82 * t75, t122, -t141, qJDD(5), g(3) * t85 + t68 + (t8 - t142) * qJD(5) + t96 * t82, t132 + (qJD(5) * t10 - t79) * t82 + t96 * t85, 0, 0;];
tau_reg = t6;
