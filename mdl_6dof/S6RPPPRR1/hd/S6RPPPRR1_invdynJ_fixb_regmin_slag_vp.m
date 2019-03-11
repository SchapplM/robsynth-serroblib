% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:16
% EndTime: 2019-03-09 01:30:19
% DurationCPUTime: 1.07s
% Computational Cost: add. (917->221), mult. (1543->275), div. (0->0), fcn. (934->10), ass. (0->124)
t67 = sin(pkin(9));
t43 = pkin(1) * t67 + qJ(3);
t58 = qJ(1) + pkin(9);
t51 = sin(t58);
t52 = cos(t58);
t158 = -g(1) * t51 + g(2) * t52;
t134 = pkin(1) * qJDD(1);
t35 = t43 * qJD(1);
t32 = qJD(4) + t35;
t23 = -qJD(1) * pkin(7) + t32;
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t13 = -qJD(2) * t70 + t23 * t73;
t10 = -qJD(5) * pkin(5) - t13;
t68 = cos(pkin(9));
t46 = -pkin(1) * t68 - pkin(2);
t39 = qJ(4) - t46;
t96 = pkin(5) * t70 - pkin(8) * t73;
t22 = t39 + t96;
t12 = t22 * qJD(1) - qJD(3);
t131 = qJD(5) * t13;
t137 = qJDD(1) * qJ(3) + t67 * t134;
t61 = qJD(3) * qJD(1);
t110 = -t61 - t137;
t98 = qJDD(4) - t110;
t19 = -qJDD(1) * pkin(7) + t98;
t102 = -qJDD(5) * pkin(8) - qJD(6) * t12 - qJDD(2) * t73 - t19 * t70 - t131;
t127 = qJD(5) * t73;
t116 = qJD(1) * qJD(5);
t105 = t73 * t116;
t118 = t70 * qJDD(1);
t26 = qJDD(6) + t105 + t118;
t38 = -pkin(7) + t43;
t14 = qJD(2) * t73 + t23 * t70;
t130 = qJD(5) * t14;
t4 = -qJDD(5) * pkin(5) + t70 * qJDD(2) - t73 * t19 + t130;
t40 = qJD(1) * t70 + qJD(6);
t157 = -(qJD(6) * t22 + t38 * t127) * t40 + t4 * t73 + (-qJD(3) * t40 - qJD(5) * t10 - t26 * t38 + t102) * t70;
t72 = cos(qJ(6));
t138 = t72 * t73;
t123 = t72 * qJD(5);
t124 = qJD(6) * t73;
t69 = sin(qJ(6));
t85 = t70 * t123 + t69 * t124;
t156 = t26 * t138 - t85 * t40;
t106 = t70 * t116;
t117 = t73 * qJDD(1);
t132 = qJD(1) * t73;
t30 = qJD(5) * t69 + t72 * t132;
t9 = qJD(6) * t30 - t72 * qJDD(5) + (-t106 + t117) * t69;
t133 = qJD(1) * t39;
t24 = -qJD(3) + t133;
t155 = (qJD(3) + t24 + t133) * qJD(5) + qJDD(5) * t38;
t95 = g(1) * t52 + g(2) * t51;
t97 = pkin(5) * t73 + pkin(8) * t70;
t154 = (pkin(8) * qJD(6) + t97 * qJD(1)) * t40 + t95 * t73 - g(3) * t70 + t4;
t151 = t70 * t9;
t8 = -t85 * qJD(1) + qJD(6) * t123 + t69 * qJDD(5) + t72 * t117;
t150 = t8 * t69;
t149 = t30 * t127 + t8 * t70;
t148 = t22 * t26;
t28 = t69 * t132 - t123;
t147 = t28 * t40;
t146 = t28 * t73;
t145 = t30 * t40;
t144 = t30 * t72;
t143 = t30 * t73;
t142 = t69 * t26;
t141 = t69 * t70;
t140 = t70 * t72;
t139 = t72 * t26;
t64 = t73 ^ 2;
t136 = t70 ^ 2 - t64;
t75 = qJD(5) ^ 2;
t76 = qJD(1) ^ 2;
t135 = -t75 - t76;
t129 = qJD(5) * t28;
t128 = qJD(5) * t70;
t126 = qJD(6) * t70;
t125 = qJD(6) * t72;
t65 = qJDD(2) - g(3);
t122 = qJDD(1) * t39;
t120 = qJDD(5) * t70;
t119 = qJDD(5) * t73;
t115 = qJD(4) * qJD(1);
t113 = t40 * t141;
t112 = t40 * t140;
t74 = cos(qJ(1));
t111 = t74 * pkin(1) + t52 * pkin(2) + t51 * qJ(3);
t108 = -qJDD(1) * pkin(2) - t68 * t134 + qJDD(3);
t71 = sin(qJ(1));
t107 = -t71 * pkin(1) + t52 * qJ(3);
t11 = qJD(5) * pkin(8) + t14;
t27 = t97 * qJD(5) + qJD(4);
t59 = qJDD(1) * qJ(4);
t99 = -t59 + t108;
t6 = t27 * qJD(1) + t96 * qJDD(1) - t99;
t103 = qJD(6) * t11 - t6;
t101 = -0.2e1 * t105;
t100 = qJD(1) + t126;
t94 = g(1) * t71 - g(2) * t74;
t93 = -qJD(5) * t23 - t65;
t90 = g(3) * t73 + t102;
t89 = t40 * t125 + t142;
t88 = -qJD(6) * t69 * t40 + t139;
t86 = t108 + t158;
t84 = qJD(1) * t24 + t95;
t83 = -t59 + t86;
t82 = qJDD(1) * t43 + t137 + 0.2e1 * t61 - t95;
t81 = -pkin(8) * t26 + (t10 + t13) * t40;
t80 = qJD(2) * qJD(5) - t19 + t84;
t20 = -t99 + t115;
t79 = -t38 * t75 + t115 + t122 - t158 + t20;
t34 = -t70 * t75 + t119;
t33 = -t73 * t75 - t120;
t25 = qJD(5) * t113;
t18 = t52 * t140 - t51 * t69;
t17 = -t52 * t141 - t51 * t72;
t16 = -t51 * t140 - t52 * t69;
t15 = t51 * t141 - t52 * t72;
t5 = t72 * t6;
t2 = t11 * t72 + t12 * t69;
t1 = -t11 * t69 + t12 * t72;
t3 = [qJDD(1), t94, g(1) * t74 + g(2) * t71 (t94 + (t67 ^ 2 + t68 ^ 2) * t134) * pkin(1), qJDD(1) * t46 + t86, t82, -t110 * t43 + t35 * qJD(3) + t108 * t46 - g(1) * (-pkin(2) * t51 + t107) - g(2) * t111, qJDD(4) + t82, 0.2e1 * t115 - t83 + t122, t20 * t39 + t24 * qJD(4) + t98 * t43 + t32 * qJD(3) - g(1) * ((-pkin(2) - qJ(4)) * t51 + t107) - g(2) * (qJ(4) * t52 + t111) qJDD(1) * t64 + t70 * t101, 0.2e1 * t136 * t116 - 0.2e1 * t70 * t117, t34, t33, 0, t155 * t73 + t79 * t70, -t155 * t70 + t79 * t73, t8 * t138 - t85 * t30 (t28 * t72 + t30 * t69) * t128 + (-t150 - t72 * t9 + (t28 * t69 - t144) * qJD(6)) * t73, t149 + t156, -t151 + t25 + (-t89 - t129) * t73, t40 * t127 + t26 * t70, -g(1) * t16 - g(2) * t18 + (t38 * t129 + t5) * t70 + (-qJD(3) * t28 + qJD(5) * t1 - t38 * t9) * t73 + (t148 + t27 * t40 + (t10 * t73 + (-t38 * t40 - t11) * t70) * qJD(6)) * t72 + t157 * t69, t38 * t30 * t128 - g(1) * t15 - g(2) * t17 + (-qJD(3) * t30 - qJD(5) * t2 - t38 * t8) * t73 + (-(-t38 * t126 + t27) * t40 - t148 + t103 * t70 - t10 * t124) * t69 + t157 * t72; 0, 0, 0, t65, 0, 0, t65, 0, 0, t65, 0, 0, 0, 0, 0, t33, -t34, 0, 0, 0, 0, 0, t151 + t25 + (-t89 + t129) * t73, t149 - t156; 0, 0, 0, 0, qJDD(1), -t76, -qJD(1) * t35 + t86, -t76, -qJDD(1) (-qJD(4) - t32) * qJD(1) + t83, 0, 0, 0, 0, 0, t101 - t118, 0.2e1 * t106 - t117, 0, 0, 0, 0, 0 (t113 + t146) * qJD(1) - t88 (t112 + t143) * qJD(1) + t89; 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t76, -t84 + t98, 0, 0, 0, 0, 0, t135 * t70 + t119, t135 * t73 - t120, 0, 0, 0, 0, 0, -t73 * t9 + (t129 - t142) * t70 + (-t100 * t72 - t69 * t127) * t40, -t73 * t8 + (qJD(5) * t30 - t139) * t70 + (t100 * t69 - t73 * t123) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t76 * t70, -t136 * t76, t117, -t118, qJDD(5), t93 * t70 - t80 * t73 + t130, t80 * t70 + t93 * t73 + t131, t40 * t144 + t150 (t8 - t147) * t72 + (-t9 - t145) * t69 (t112 - t143) * qJD(1) + t89 (-t113 + t146) * qJD(1) + t88, -t40 * t132, -pkin(5) * t9 - t1 * t132 - t14 * t28 - t154 * t72 + t81 * t69, -pkin(5) * t8 + t2 * t132 - t14 * t30 + t154 * t69 + t81 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t28, -t28 ^ 2 + t30 ^ 2, t8 + t147, t145 - t9, t26, -g(1) * t17 + g(2) * t15 - t10 * t30 - t11 * t125 + t2 * t40 + t90 * t69 + t5, g(1) * t18 - g(2) * t16 + t1 * t40 + t10 * t28 + t103 * t69 + t90 * t72;];
tau_reg  = t3;
