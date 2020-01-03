% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:39
% EndTime: 2020-01-03 12:07:40
% DurationCPUTime: 0.72s
% Computational Cost: add. (1247->163), mult. (2126->221), div. (0->0), fcn. (1241->16), ass. (0->122)
t133 = pkin(1) * qJD(1);
t86 = sin(qJ(2));
t123 = t86 * t133;
t90 = cos(qJ(2));
t149 = t90 * pkin(1);
t70 = qJDD(1) * t149;
t77 = qJDD(1) + qJDD(2);
t36 = t77 * pkin(2) - qJD(2) * t123 + t70;
t78 = qJD(1) + qJD(2);
t48 = t78 * pkin(2) + t90 * t133;
t85 = sin(qJ(3));
t89 = cos(qJ(3));
t131 = qJD(2) * t90;
t119 = qJD(1) * t131;
t124 = qJDD(1) * t86;
t99 = (t119 + t124) * pkin(1);
t156 = -(qJD(3) * t48 + t99) * t89 - t85 * t36;
t128 = qJD(3) * t86;
t118 = qJD(1) * t128;
t114 = pkin(1) * t118;
t51 = t85 * t114;
t15 = -t156 - t51;
t82 = sin(pkin(9));
t83 = cos(pkin(9));
t72 = qJDD(3) + t77;
t127 = qJD(3) * t89;
t121 = t86 * t127;
t129 = qJD(3) * t85;
t122 = t48 * t129;
t32 = t89 * t36;
t94 = t32 + (-t85 * t124 + (-t85 * t131 - t121) * qJD(1)) * pkin(1) - t122;
t9 = t72 * pkin(3) + t94;
t4 = -t82 * t15 + t83 * t9;
t81 = qJ(1) + qJ(2);
t76 = qJ(3) + t81;
t65 = pkin(9) + t76;
t54 = sin(t65);
t55 = cos(t65);
t154 = -t72 * pkin(4) + g(2) * t55 + g(3) * t54 - t4;
t132 = pkin(2) * qJD(3);
t145 = t83 * t85;
t142 = t86 * t89;
t109 = -t85 * t90 - t142;
t42 = t109 * t133;
t143 = t85 * t86;
t108 = t89 * t90 - t143;
t43 = t108 * t133;
t139 = -t83 * t42 + t82 * t43 - (t82 * t89 + t145) * t132;
t147 = t82 * t85;
t68 = t89 * pkin(2) + pkin(3);
t104 = -pkin(2) * t147 + t83 * t68;
t37 = -pkin(4) - t104;
t135 = pkin(2) * t145 + t82 * t68;
t38 = pkin(8) + t135;
t73 = qJD(3) + t78;
t92 = qJD(5) ^ 2;
t153 = -t139 * t73 + t37 * t72 + t38 * t92;
t152 = pkin(2) * t72;
t5 = t83 * t15 + t82 * t9;
t31 = t89 * t123 + t85 * t48;
t148 = t82 * t31;
t146 = t83 * t31;
t88 = cos(qJ(5));
t141 = t88 * t72;
t69 = pkin(2) + t149;
t53 = t89 * t69;
t41 = -pkin(1) * t143 + pkin(3) + t53;
t44 = pkin(1) * t142 + t85 * t69;
t140 = t82 * t41 + t83 * t44;
t138 = -t82 * t42 - t83 * t43 + (t83 * t89 - t147) * t132;
t66 = sin(t76);
t74 = sin(t81);
t137 = pkin(2) * t74 + pkin(3) * t66;
t67 = cos(t76);
t75 = cos(t81);
t136 = pkin(2) * t75 + pkin(3) * t67;
t84 = sin(qJ(5));
t79 = t84 ^ 2;
t134 = -t88 ^ 2 + t79;
t126 = t88 * qJD(5);
t125 = qJDD(4) - g(1);
t120 = g(2) * t74 - g(3) * t75;
t30 = -t85 * t123 + t89 * t48;
t28 = t73 * pkin(3) + t30;
t16 = t83 * t28 - t148;
t13 = -t73 * pkin(4) - t16;
t117 = t13 * t126 + t154 * t84;
t116 = qJD(1) * (-qJD(2) + t78);
t115 = qJD(2) * (-qJD(1) - t78);
t113 = g(2) * t66 - g(3) * t67 + t51;
t112 = -g(2) * t67 - g(3) * t66;
t110 = t83 * t41 - t82 * t44;
t106 = t112 + t32;
t105 = -g(2) * t75 - g(3) * t74 + t70;
t20 = -pkin(4) - t110;
t21 = pkin(8) + t140;
t24 = t69 * t127 + (t108 * qJD(2) - t85 * t128) * pkin(1);
t25 = -t69 * t129 + (t109 * qJD(2) - t121) * pkin(1);
t6 = t82 * t24 - t83 * t25;
t102 = t20 * t72 + t21 * t92 + t6 * t73;
t18 = t82 * t30 + t146;
t60 = t82 * pkin(3) + pkin(8);
t61 = -t83 * pkin(3) - pkin(4);
t101 = -t18 * t73 + t60 * t92 + t61 * t72;
t100 = -t72 * pkin(8) + g(2) * t54 - g(3) * t55 - t13 * t73 - t5;
t7 = t83 * t24 + t82 * t25;
t98 = -qJDD(5) * t21 + (t20 * t73 - t7) * qJD(5);
t19 = t83 * t30 - t148;
t97 = -qJDD(5) * t60 + (t61 * t73 + t19) * qJD(5);
t96 = -qJDD(5) * t38 + (t37 * t73 - t138) * qJD(5);
t95 = (-pkin(2) * t73 - t48) * qJD(3) - t99;
t93 = t113 + t156;
t91 = cos(qJ(1));
t87 = sin(qJ(1));
t71 = t73 ^ 2;
t50 = qJDD(5) * t88 - t92 * t84;
t49 = qJDD(5) * t84 + t92 * t88;
t33 = 0.2e1 * t84 * t73 * t126 + t79 * t72;
t26 = -0.2e1 * t134 * t73 * qJD(5) + 0.2e1 * t84 * t141;
t17 = t82 * t28 + t146;
t10 = t13 * qJD(5) * t84;
t1 = [qJDD(1), -g(2) * t91 - g(3) * t87, g(2) * t87 - g(3) * t91, t77, (t86 * t115 + t77 * t90) * pkin(1) + t105, ((-qJDD(1) - t77) * t86 + t90 * t115) * pkin(1) + t120, t72, -t122 + t25 * t73 + t53 * t72 + (-t89 * t118 + (-t119 + (-qJDD(1) - t72) * t86) * t85) * pkin(1) + t106, -t24 * t73 - t44 * t72 + t93, t5 * t140 + t17 * t7 + t4 * t110 - t16 * t6 - g(2) * (t91 * pkin(1) + t136) - g(3) * (t87 * pkin(1) + t137), t33, t26, t49, t50, 0, t10 + t98 * t84 + (-t102 - t154) * t88, t102 * t84 + t98 * t88 + t117; 0, 0, 0, t77, t86 * pkin(1) * t116 + t105, (t90 * t116 - t124) * pkin(1) + t120, t72, -t42 * t73 + (-t114 + t152) * t89 + t95 * t85 + t106, t43 * t73 + (-t36 - t152) * t85 + t95 * t89 + t113, -g(2) * t136 - g(3) * t137 + t4 * t104 + t5 * t135 + t138 * t17 + t139 * t16, t33, t26, t49, t50, 0, t10 + t96 * t84 + (-t154 - t153) * t88, t153 * t84 + t96 * t88 + t117; 0, 0, 0, 0, 0, 0, t72, t31 * t73 + t112 + t94, t30 * t73 + t93, t16 * t18 - t17 * t19 + (t4 * t83 + t5 * t82 + t112) * pkin(3), t33, t26, t49, t50, 0, t10 + t97 * t84 + (-t101 - t154) * t88, t101 * t84 + t97 * t88 + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, 0, 0, 0, 0, 0, t50, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84 * t71 * t88, t134 * t71, t84 * t72, t141, qJDD(5), t100 * t84 + t125 * t88, t100 * t88 - t125 * t84;];
tau_reg = t1;
