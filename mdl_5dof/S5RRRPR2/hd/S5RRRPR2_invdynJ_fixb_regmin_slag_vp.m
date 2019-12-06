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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:41:01
% EndTime: 2019-12-05 18:41:03
% DurationCPUTime: 0.74s
% Computational Cost: add. (1247->166), mult. (2126->223), div. (0->0), fcn. (1241->16), ass. (0->122)
t133 = pkin(1) * qJD(1);
t86 = sin(qJ(2));
t121 = t86 * t133;
t90 = cos(qJ(2));
t148 = t90 * pkin(1);
t70 = qJDD(1) * t148;
t77 = qJDD(1) + qJDD(2);
t36 = t77 * pkin(2) - qJD(2) * t121 + t70;
t78 = qJD(1) + qJD(2);
t48 = t78 * pkin(2) + t90 * t133;
t85 = sin(qJ(3));
t89 = cos(qJ(3));
t131 = qJD(2) * t90;
t118 = qJD(1) * t131;
t124 = qJDD(1) * t86;
t99 = (t118 + t124) * pkin(1);
t156 = -(qJD(3) * t48 + t99) * t89 - t85 * t36;
t81 = qJ(1) + qJ(2);
t76 = qJ(3) + t81;
t65 = pkin(9) + t76;
t54 = sin(t65);
t55 = cos(t65);
t154 = g(2) * t55 + g(3) * t54;
t66 = sin(t76);
t67 = cos(t76);
t153 = g(2) * t67 + g(3) * t66;
t132 = pkin(2) * qJD(3);
t83 = cos(pkin(9));
t143 = t83 * t85;
t140 = t86 * t89;
t106 = -t85 * t90 - t140;
t42 = t106 * t133;
t141 = t85 * t86;
t105 = t89 * t90 - t141;
t43 = t105 * t133;
t82 = sin(pkin(9));
t137 = -t83 * t42 + t82 * t43 - (t82 * t89 + t143) * t132;
t145 = t82 * t85;
t68 = t89 * pkin(2) + pkin(3);
t104 = -pkin(2) * t145 + t83 * t68;
t37 = -pkin(4) - t104;
t135 = pkin(2) * t143 + t82 * t68;
t38 = pkin(8) + t135;
t72 = qJDD(3) + t77;
t73 = qJD(3) + t78;
t92 = qJD(5) ^ 2;
t152 = -t137 * t73 + t37 * t72 + t38 * t92;
t151 = pkin(2) * t72;
t88 = cos(qJ(5));
t126 = t88 * qJD(5);
t31 = t89 * t121 + t85 * t48;
t146 = t82 * t31;
t30 = -t85 * t121 + t89 * t48;
t28 = t73 * pkin(3) + t30;
t16 = t83 * t28 - t146;
t13 = -t73 * pkin(4) - t16;
t128 = qJD(3) * t86;
t117 = qJD(1) * t128;
t114 = pkin(1) * t117;
t51 = t85 * t114;
t15 = -t156 - t51;
t129 = qJD(3) * t85;
t32 = t89 * t36;
t109 = -t48 * t129 + t32;
t127 = qJD(3) * t89;
t120 = t86 * t127;
t94 = (-t85 * t124 + (-t85 * t131 - t120) * qJD(1)) * pkin(1) + t109;
t9 = t72 * pkin(3) + t94;
t4 = -t82 * t15 + t83 * t9;
t2 = -t72 * pkin(4) - t4;
t84 = sin(qJ(5));
t147 = t13 * t126 + t2 * t84;
t5 = t83 * t15 + t82 * t9;
t144 = t83 * t31;
t139 = t88 * t72;
t69 = pkin(2) + t148;
t53 = t89 * t69;
t41 = -pkin(1) * t141 + pkin(3) + t53;
t44 = pkin(1) * t140 + t85 * t69;
t138 = t82 * t41 + t83 * t44;
t136 = -t82 * t42 - t83 * t43 + (t83 * t89 - t145) * t132;
t79 = t84 ^ 2;
t134 = -t88 ^ 2 + t79;
t125 = qJDD(4) - g(1);
t123 = t13 * qJD(5) * t84 + t154 * t88;
t74 = sin(t81);
t75 = cos(t81);
t122 = g(2) * t75 + g(3) * t74 + t70;
t119 = -g(2) * t74 + g(3) * t75;
t116 = qJD(1) * (-qJD(2) + t78);
t115 = qJD(2) * (-qJD(1) - t78);
t113 = -g(2) * t66 + g(3) * t67 + t51;
t112 = -pkin(2) * t74 - pkin(3) * t66;
t111 = -pkin(2) * t75 - pkin(3) * t67;
t107 = t83 * t41 - t82 * t44;
t20 = -pkin(4) - t107;
t21 = pkin(8) + t138;
t24 = t69 * t127 + (t105 * qJD(2) - t85 * t128) * pkin(1);
t25 = -t69 * t129 + (t106 * qJD(2) - t120) * pkin(1);
t6 = t82 * t24 - t83 * t25;
t102 = -t20 * t72 - t21 * t92 - t6 * t73;
t18 = t82 * t30 + t144;
t60 = t82 * pkin(3) + pkin(8);
t61 = -t83 * pkin(3) - pkin(4);
t101 = t18 * t73 - t60 * t92 - t61 * t72;
t100 = -t72 * pkin(8) - g(2) * t54 + g(3) * t55 - t13 * t73 - t5;
t7 = t83 * t24 + t82 * t25;
t98 = -qJDD(5) * t21 + (t20 * t73 - t7) * qJD(5);
t19 = t83 * t30 - t146;
t97 = -qJDD(5) * t60 + (t61 * t73 + t19) * qJD(5);
t96 = -qJDD(5) * t38 + (t37 * t73 - t136) * qJD(5);
t95 = (-pkin(2) * t73 - t48) * qJD(3) - t99;
t93 = t113 + t156;
t91 = cos(qJ(1));
t87 = sin(qJ(1));
t71 = t73 ^ 2;
t50 = qJDD(5) * t88 - t92 * t84;
t49 = qJDD(5) * t84 + t92 * t88;
t33 = 0.2e1 * t84 * t73 * t126 + t79 * t72;
t26 = -0.2e1 * t134 * t73 * qJD(5) + 0.2e1 * t84 * t139;
t17 = t82 * t28 + t144;
t1 = [qJDD(1), g(2) * t91 + g(3) * t87, -g(2) * t87 + g(3) * t91, t77, (t86 * t115 + t77 * t90) * pkin(1) + t122, ((-qJDD(1) - t77) * t86 + t90 * t115) * pkin(1) + t119, t72, t25 * t73 + t53 * t72 + (-t89 * t117 + (-t118 + (-qJDD(1) - t72) * t86) * t85) * pkin(1) + t109 + t153, -t24 * t73 - t44 * t72 + t93, t5 * t138 + t17 * t7 + t4 * t107 - t16 * t6 - g(2) * (-t91 * pkin(1) + t111) - g(3) * (-t87 * pkin(1) + t112), t33, t26, t49, t50, 0, t98 * t84 + (t102 - t2) * t88 + t123, t98 * t88 + (-t102 - t154) * t84 + t147; 0, 0, 0, t77, t86 * pkin(1) * t116 + t122, (t90 * t116 - t124) * pkin(1) + t119, t72, -t42 * t73 + t32 + (-t114 + t151) * t89 + t95 * t85 + t153, t43 * t73 + (-t36 - t151) * t85 + t95 * t89 + t113, -g(2) * t111 - g(3) * t112 + t4 * t104 + t5 * t135 + t136 * t17 + t137 * t16, t33, t26, t49, t50, 0, t96 * t84 + (-t152 - t2) * t88 + t123, t96 * t88 + (-t154 + t152) * t84 + t147; 0, 0, 0, 0, 0, 0, t72, t31 * t73 + t153 + t94, t30 * t73 + t93, t16 * t18 - t17 * t19 + (t4 * t83 + t5 * t82 + t153) * pkin(3), t33, t26, t49, t50, 0, t97 * t84 + (t101 - t2) * t88 + t123, t97 * t88 + (-t101 - t154) * t84 + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, 0, 0, 0, 0, 0, t50, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84 * t71 * t88, t134 * t71, t84 * t72, t139, qJDD(5), t100 * t84 + t125 * t88, t100 * t88 - t125 * t84;];
tau_reg = t1;
