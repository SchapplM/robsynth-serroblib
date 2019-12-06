% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:24
% EndTime: 2019-12-05 18:20:27
% DurationCPUTime: 1.02s
% Computational Cost: add. (1213->201), mult. (1950->281), div. (0->0), fcn. (1204->14), ass. (0->138)
t133 = pkin(1) * qJD(2);
t120 = qJD(1) * t133;
t88 = sin(qJ(2));
t128 = qJDD(1) * t88;
t91 = cos(qJ(2));
t169 = pkin(1) * t128 + t91 * t120;
t82 = qJ(1) + qJ(2);
t72 = pkin(8) + t82;
t62 = sin(t72);
t63 = cos(t72);
t168 = g(2) * t63 + g(3) * t62;
t156 = t91 * pkin(1);
t69 = qJDD(1) * t156;
t76 = qJDD(1) + qJDD(2);
t34 = t76 * pkin(2) - t88 * t120 + t69;
t84 = sin(pkin(8));
t86 = cos(pkin(8));
t18 = t169 * t86 + t84 * t34;
t79 = qJD(1) + qJD(2);
t12 = t76 * qJ(4) + t79 * qJD(4) + t18;
t83 = sin(pkin(9));
t85 = cos(pkin(9));
t8 = -t85 * qJDD(3) + t83 * t12;
t159 = t8 * t83;
t9 = t83 * qJDD(3) + t85 * t12;
t167 = t9 * t85 + t159;
t73 = sin(t82);
t74 = cos(t82);
t166 = g(2) * t74 + g(3) * t73;
t134 = pkin(1) * qJD(1);
t124 = t91 * t134;
t125 = t88 * t134;
t54 = t84 * t125;
t39 = t86 * t124 - t54;
t130 = qJD(4) - t39;
t145 = t85 * t79;
t51 = -qJD(5) + t145;
t129 = qJD(5) + t51;
t43 = t79 * pkin(2) + t124;
t55 = t86 * t125;
t27 = t84 * t43 + t55;
t23 = t79 * qJ(4) + t27;
t19 = -t85 * qJD(3) + t83 * t23;
t154 = t19 * t83;
t20 = t83 * qJD(3) + t85 * t23;
t165 = -t129 * t20 - t79 * t154;
t164 = pkin(1) * t88;
t163 = pkin(2) * t73;
t162 = pkin(2) * t74;
t158 = t86 * pkin(2);
t89 = sin(qJ(1));
t157 = t89 * pkin(1);
t92 = cos(qJ(1));
t155 = t92 * pkin(1);
t59 = t84 * t164;
t104 = t86 * t91 * t133 - qJD(2) * t59;
t33 = qJD(4) + t104;
t153 = t33 * t79;
t87 = sin(qJ(5));
t152 = t33 * t87;
t146 = t85 * t76;
t48 = -qJDD(5) + t146;
t151 = t48 * t85;
t150 = t76 * t87;
t90 = cos(qJ(5));
t149 = t76 * t90;
t75 = t79 ^ 2;
t77 = t83 ^ 2;
t148 = t77 * t75;
t147 = t83 * t87;
t144 = t85 * t87;
t143 = t85 * t90;
t142 = t86 * t88;
t141 = t87 * t48;
t140 = t87 * t90;
t139 = t90 * t48;
t138 = t168 * t85;
t68 = pkin(2) + t156;
t137 = pkin(1) * t142 + t84 * t68;
t136 = t85 ^ 2 + t77;
t81 = t90 ^ 2;
t135 = t87 ^ 2 - t81;
t132 = qJD(5) * t87;
t131 = qJD(5) * t90;
t126 = t69 + t166;
t123 = t79 * t131;
t122 = t51 * t132;
t119 = -g(2) * t73 + g(3) * t74;
t118 = t136 * t76;
t26 = t86 * t43 - t54;
t117 = t48 - t146;
t116 = t48 + t146;
t115 = t86 * t68 - t59;
t114 = t79 * t129;
t113 = t130 * t90;
t112 = qJD(1) * (-qJD(2) + t79);
t111 = qJD(2) * (-qJD(1) - t79);
t17 = -t169 * t84 + t86 * t34;
t108 = qJD(4) - t26;
t107 = t20 * t85 + t154;
t36 = -pkin(3) - t115;
t38 = (t84 * t91 + t142) * t133;
t106 = -t36 * t76 - t38 * t79;
t37 = t84 * t124 + t55;
t64 = -pkin(3) - t158;
t105 = t37 * t79 - t64 * t76;
t103 = -t85 * pkin(4) - t83 * pkin(7) - pkin(3);
t102 = -t62 * pkin(3) + t63 * qJ(4) - t163;
t101 = qJDD(4) - t17;
t40 = t103 - t158;
t61 = t84 * pkin(2) + qJ(4);
t100 = -t61 * t144 + t90 * t40;
t15 = t103 * t79 + t108;
t28 = t62 * t144 + t63 * t90;
t30 = t63 * t144 - t62 * t90;
t7 = t103 * t76 + t101;
t99 = -g(2) * t30 - g(3) * t28 + (t87 * t7 + t90 * t9 + (t90 * t15 - t87 * t20) * qJD(5)) * t85 + t90 * t159;
t29 = t62 * t143 - t63 * t87;
t31 = -t63 * t143 - t62 * t87;
t98 = -g(2) * t31 + g(3) * t29 + t131 * t154 + t8 * t147;
t97 = g(2) * t62 - g(3) * t63 + t167;
t96 = -t63 * pkin(3) - t62 * qJ(4) - t162;
t95 = g(1) * t83 - t129 * t15 - t9;
t94 = t130 * t87 + t61 * t131;
t14 = -t76 * pkin(3) + t101;
t93 = -t51 ^ 2 - t148;
t42 = t83 * t132 * t145;
t35 = qJ(4) + t137;
t25 = (-0.2e1 * t87 * t123 + t76 * t81) * t77;
t24 = t103 - t115;
t22 = -t79 * pkin(3) + t108;
t21 = 0.2e1 * (t135 * t79 * qJD(5) - t76 * t140) * t77;
t13 = t14 * t83;
t11 = (t116 * t87 + (t51 + t145) * t131) * t83;
t10 = t42 + (-t116 * t90 + t122) * t83;
t5 = t90 * t7;
t2 = -t87 * t9 + t5 + (-t87 * t15 - t90 * t20) * qJD(5);
t1 = [qJDD(1), g(2) * t92 + g(3) * t89, -g(2) * t89 + g(3) * t92, t76, (t88 * t111 + t76 * t91) * pkin(1) + t126, ((-qJDD(1) - t76) * t88 + t91 * t111) * pkin(1) + t119, t18 * t137 + t27 * t104 + t17 * t115 - t26 * t38 - g(2) * (-t155 - t162) - g(3) * (-t157 - t163), (t106 - t14) * t85 + t138, t13 + (-t106 - t168) * t83, t35 * t118 + t136 * t153 + t97, t14 * t36 + t22 * t38 - g(2) * (t96 - t155) - g(3) * (t102 - t157) + t167 * t35 + t107 * t33, t25, t21, t10, t11, t151, -(-t24 * t132 + t90 * t38) * t51 - t24 * t139 + (-(-t35 * t131 - t152) * t51 + t35 * t141 - t2) * t85 + (t79 * t152 + (t123 + t150) * t35) * t77 + t98, (t33 * t143 + t87 * t38) * t51 + (t35 * t143 + t87 * t24) * t48 + (t35 * t76 + t153) * t90 * t77 + (t90 * t24 * t51 + (-t154 + (-t51 * t85 - t77 * t79) * t35) * t87) * qJD(5) + t99; 0, 0, 0, t76, t112 * t164 + t126, (t91 * t112 - t128) * pkin(1) + t119, t26 * t37 - t27 * t39 + (t17 * t86 + t18 * t84 + t166) * pkin(2), (t105 - t14) * t85 + t138, t13 + (-t105 - t168) * t83, t130 * t79 * t136 + t61 * t118 + t97, t14 * t64 - t22 * t37 - g(2) * t96 - g(3) * t102 + (t130 * t20 + t9 * t61) * t85 + (t130 * t19 + t8 * t61) * t83, t25, t21, t10, t11, t151, -t100 * t48 - t2 * t85 + (t40 * t132 + t90 * t37 + t94 * t85) * t51 + (t61 * t150 + t94 * t79) * t77 + t98, (t61 * t143 + t87 * t40) * t48 + (t85 * t113 - t87 * t37) * t51 + (t100 * t51 - t19 * t147) * qJD(5) + (t61 * t149 + (-t61 * t132 + t113) * t79) * t77 + t99; 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, -t8 * t85 + t9 * t83 - g(1), 0, 0, 0, 0, 0, (t117 * t87 + (t51 - t145) * t131) * t83, t42 + (t117 * t90 - t122) * t83; 0, 0, 0, 0, 0, 0, 0, -t146, t83 * t76, -t136 * t75, -t107 * t79 + t14 - t168, 0, 0, 0, 0, 0, t93 * t87 - t139, t93 * t90 + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140 * t148, -t135 * t148, (-t87 * t114 + t149) * t83, (-t90 * t114 - t150) * t83, -t48, -g(2) * t28 + g(3) * t30 + t165 * t90 + t95 * t87 + t5, -g(2) * t29 - g(3) * t31 + t95 * t90 + (-t165 - t7) * t87;];
tau_reg = t1;
