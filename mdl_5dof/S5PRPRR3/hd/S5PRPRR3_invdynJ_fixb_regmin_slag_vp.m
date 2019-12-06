% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:37
% EndTime: 2019-12-05 15:47:43
% DurationCPUTime: 1.02s
% Computational Cost: add. (788->170), mult. (1634->244), div. (0->0), fcn. (1272->14), ass. (0->120)
t92 = cos(qJ(4));
t123 = -pkin(4) * t92 - pkin(3);
t93 = cos(qJ(2));
t132 = t93 * qJD(1);
t57 = qJD(2) * pkin(2) + t132;
t90 = sin(qJ(2));
t142 = qJD(1) * t90;
t84 = sin(pkin(9));
t60 = t84 * t142;
t86 = cos(pkin(9));
t28 = t86 * t57 - t60;
t21 = qJD(2) * t123 - t28;
t36 = t132 * t86 - t60;
t162 = t21 + t36;
t85 = sin(pkin(8));
t87 = cos(pkin(8));
t116 = g(1) * t87 + g(2) * t85;
t80 = qJ(2) + pkin(9);
t71 = sin(t80);
t72 = cos(t80);
t159 = -g(3) * t72 + t116 * t71;
t79 = qJD(4) + qJD(5);
t127 = qJD(1) * qJD(2);
t160 = t90 * qJDD(1) + t93 * t127;
t74 = t93 * qJDD(1);
t40 = qJDD(2) * pkin(2) - t127 * t90 + t74;
t17 = -t160 * t84 + t86 * t40;
t13 = -qJDD(2) * pkin(3) - t17;
t61 = t86 * t142;
t34 = t132 * t84 + t61;
t68 = pkin(2) * t84 + pkin(6);
t153 = t86 * pkin(2);
t69 = -pkin(3) - t153;
t94 = qJD(4) ^ 2;
t161 = qJD(2) * t34 - qJDD(2) * t69 - t68 * t94 - t13 + t159;
t128 = t92 * qJDD(2);
t89 = sin(qJ(4));
t130 = t89 * qJDD(2);
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t46 = t88 * t92 + t89 * t91;
t158 = t46 * t79;
t9 = qJD(2) * t158 - t128 * t91 + t130 * t88;
t45 = t88 * t89 - t91 * t92;
t104 = t45 * t79;
t29 = t57 * t84 + t61;
t119 = t29 + (pkin(6) + pkin(7)) * qJD(2);
t15 = t92 * qJD(3) - t119 * t89;
t156 = t79 * t92;
t16 = t89 * qJD(3) + t119 * t92;
t155 = g(3) * t71;
t152 = pkin(7) + t68;
t139 = qJD(2) * t92;
t124 = t91 * t139;
t140 = qJD(2) * t89;
t125 = t88 * t140;
t37 = -t124 + t125;
t39 = -t139 * t88 - t140 * t91;
t151 = t39 * t37;
t78 = qJDD(4) + qJDD(5);
t150 = t45 * t78;
t149 = t46 * t78;
t83 = qJ(4) + qJ(5);
t76 = sin(t83);
t148 = t85 * t76;
t77 = cos(t83);
t147 = t85 * t77;
t146 = t87 * t76;
t145 = t87 * t77;
t144 = t91 * t16;
t81 = t89 ^ 2;
t143 = -t92 ^ 2 + t81;
t43 = t84 * t90 - t86 * t93;
t141 = qJD(2) * t43;
t138 = qJD(4) * t89;
t136 = qJD(5) * t88;
t135 = qJD(5) * t91;
t131 = qJDD(1) - g(3);
t126 = qJD(2) * qJD(4);
t18 = t160 * t86 + t84 * t40;
t121 = t92 * t126;
t120 = qJD(4) * t152;
t14 = qJDD(2) * pkin(6) + t18;
t118 = pkin(7) * qJDD(2) + t14;
t117 = pkin(4) * t138 - t34;
t115 = -g(1) * t85 + g(2) * t87;
t12 = qJD(4) * pkin(4) + t15;
t113 = -t12 * t88 - t144;
t112 = -t104 * t79 + t149;
t41 = t152 * t89;
t42 = t152 * t92;
t111 = -t41 * t91 - t42 * t88;
t110 = -t41 * t88 + t42 * t91;
t44 = t84 * t93 + t86 * t90;
t108 = -qJDD(3) - t115;
t33 = t44 * qJD(2);
t103 = qJD(2) * t33 + qJDD(2) * t43 + t44 * t94;
t102 = 0.2e1 * qJD(4) * t141 - qJDD(4) * t44;
t101 = -g(3) * t93 + t116 * t90;
t22 = -qJD(2) * pkin(3) - t28;
t100 = -qJDD(4) * t68 + (qJD(2) * t69 + t22 + t36) * qJD(4);
t8 = qJD(5) * t124 - t79 * t125 + t88 * t128 + (t121 + t130) * t91;
t98 = -t22 * qJD(2) + t116 * t72 - t14 + t155;
t73 = t92 * qJDD(3);
t2 = qJDD(4) * pkin(4) - qJD(4) * t16 - t118 * t89 + t73;
t97 = -g(1) * (-t145 * t72 - t148) - g(2) * (-t147 * t72 + t146) + t21 * t37 + t16 * t136 + t77 * t155 + (-t16 * t79 - t2) * t88;
t3 = qJD(4) * t15 + t89 * qJDD(3) + t118 * t92;
t96 = -g(1) * (-t146 * t72 + t147) - g(2) * (-t148 * t72 - t145) + t113 * qJD(5) + t91 * t2 + t21 * t39 - t88 * t3 + t76 * t155;
t95 = qJD(2) ^ 2;
t52 = qJDD(4) * t92 - t89 * t94;
t51 = qJDD(4) * t89 + t92 * t94;
t50 = t123 - t153;
t32 = t92 * t120;
t31 = t89 * t120;
t10 = -t37 ^ 2 + t39 ^ 2;
t7 = -t158 * t79 - t150;
t6 = (t126 * t89 - t128) * pkin(4) + t13;
t5 = -t39 * t79 - t9;
t4 = t37 * t79 + t8;
t1 = [t131, 0, qJDD(2) * t93 - t90 * t95, -qJDD(2) * t90 - t93 * t95, -t141 * t29 - t17 * t43 + t18 * t44 - t28 * t33 - g(3), 0, 0, 0, 0, 0, t102 * t89 - t103 * t92, t102 * t92 + t103 * t89, 0, 0, 0, 0, 0, t33 * t37 + t43 * t9 + t141 * t158 + ((t89 * t136 + t88 * t138 - t156 * t91) * t79 - t149) * t44, -t33 * t39 + t43 * t8 - t141 * t104 + (-(-t89 * t135 - t91 * t138 - t156 * t88) * t79 + t150) * t44; 0, qJDD(2), t74 + t101, t116 * t93 - t131 * t90, t28 * t34 - t29 * t36 + (t17 * t86 + t18 * t84 + t101) * pkin(2), qJDD(2) * t81 + 0.2e1 * t121 * t89, -0.2e1 * t126 * t143 + 0.2e1 * t128 * t89, t51, t52, 0, t100 * t89 + t161 * t92, t100 * t92 - t161 * t89, t104 * t39 + t46 * t8, t104 * t37 + t158 * t39 - t45 * t8 - t46 * t9, t112, t7, 0, (-qJD(5) * t110 + t88 * t31 - t91 * t32) * t79 + t111 * t78 + t50 * t9 + t6 * t45 + t117 * t37 + t159 * t77 + t162 * t158, -(qJD(5) * t111 - t91 * t31 - t88 * t32) * t79 - t110 * t78 + t50 * t8 + t6 * t46 - t117 * t39 - t159 * t76 - t162 * t104; 0, 0, 0, 0, -t108, 0, 0, 0, 0, 0, t52, -t51, 0, 0, 0, 0, 0, t7, -t112; 0, 0, 0, 0, 0, -t89 * t95 * t92, t143 * t95, t130, t128, qJDD(4), t115 * t92 + t89 * t98 + t73, t108 * t89 + t92 * t98, -t151, t10, t4, t5, t78, -(-t15 * t88 - t144) * t79 + (-t136 * t79 - t140 * t37 + t91 * t78) * pkin(4) + t96, (-qJD(5) * t12 + t15 * t79 - t3) * t91 + (-t135 * t79 + t140 * t39 - t88 * t78) * pkin(4) + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t10, t4, t5, t78, -t113 * t79 + t96, (-t3 + (-qJD(5) + t79) * t12) * t91 + t97;];
tau_reg = t1;
