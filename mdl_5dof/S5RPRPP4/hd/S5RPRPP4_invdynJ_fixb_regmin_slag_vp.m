% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:01
% EndTime: 2021-01-15 11:24:06
% DurationCPUTime: 1.10s
% Computational Cost: add. (1483->229), mult. (2716->268), div. (0->0), fcn. (1637->8), ass. (0->122)
t100 = cos(qJ(1));
t98 = sin(qJ(1));
t165 = g(1) * t98 - g(2) * t100;
t171 = qJDD(2) - t165;
t101 = -pkin(1) - pkin(6);
t62 = t101 * qJD(1) + qJD(2);
t127 = -qJ(4) * qJD(1) + t62;
t134 = qJD(1) * qJD(3);
t99 = cos(qJ(3));
t129 = t99 * t134;
t97 = sin(qJ(3));
t136 = t97 * qJDD(1);
t170 = t129 + t136;
t123 = g(1) * t100 + g(2) * t98;
t91 = (qJD(1) * qJD(2));
t114 = -t123 + (2 * t91);
t90 = qJDD(1) * qJ(2);
t169 = 0.2e1 * t90 + t114;
t144 = qJDD(1) * pkin(1);
t168 = t144 - t171;
t95 = sin(pkin(7));
t96 = cos(pkin(7));
t52 = t95 * t99 + t96 * t97;
t164 = t52 * qJD(1);
t167 = t164 * qJD(3);
t148 = qJD(1) * t97;
t131 = t95 * t148;
t147 = qJD(1) * t99;
t48 = t96 * t147 - t131;
t43 = t48 ^ 2;
t166 = -t164 ^ 2 - t43;
t163 = qJ(4) - t101;
t146 = qJD(3) * t97;
t113 = -t99 * qJD(4) + t146 * t163;
t128 = t163 * t99;
t39 = -qJD(3) * t128 - t97 * qJD(4);
t21 = -t96 * t113 + t95 * t39;
t22 = t95 * t113 + t96 * t39;
t130 = t97 * t134;
t135 = t99 * qJDD(1);
t132 = t95 * t135 + t170 * t96;
t28 = t95 * t130 - t132;
t122 = -t96 * t135 + t95 * t136;
t29 = t122 + t167;
t57 = t163 * t97;
t31 = t96 * t128 - t95 * t57;
t32 = -t95 * t128 - t96 * t57;
t162 = -t164 * t22 + t21 * t48 + t32 * t28 - t31 * t29;
t88 = qJ(3) + pkin(7);
t82 = cos(t88);
t161 = t22 * qJD(3) + t32 * qJDD(3) + t123 * t82;
t133 = qJD(1) * qJD(4);
t61 = t101 * qJDD(1) + qJDD(2);
t53 = t99 * t61;
t17 = -t99 * t133 - t62 * t146 + qJDD(3) * pkin(3) + t53 + (t130 - t135) * qJ(4);
t145 = qJD(3) * t99;
t23 = t127 * t145 + (-qJ(4) * qJDD(1) - t133 + t61) * t97;
t5 = t95 * t17 + t96 * t23;
t81 = sin(t88);
t160 = g(3) * t82 + t165 * t81 - t5;
t157 = g(3) * t97;
t156 = t97 * pkin(3);
t55 = pkin(3) * t148 + qJD(1) * qJ(2) + qJD(4);
t20 = pkin(4) * t164 - t48 * qJ(5) + t55;
t155 = t20 * t48;
t152 = t163 * t98;
t40 = t127 * t97;
t151 = t95 * t40;
t36 = t96 * t40;
t4 = t96 * t17 - t95 * t23;
t41 = t127 * t99;
t38 = qJD(3) * pkin(3) + t41;
t19 = t95 * t38 + t36;
t94 = t99 ^ 2;
t149 = t97 ^ 2 - t94;
t103 = qJD(1) ^ 2;
t143 = t103 * qJ(2);
t142 = t55 * qJD(1);
t25 = t96 * t41 - t151;
t140 = qJD(5) - t25;
t102 = qJD(3) ^ 2;
t139 = -t102 - t103;
t66 = pkin(3) * t145 + qJD(2);
t137 = qJDD(3) * t97;
t76 = qJ(2) + t156;
t18 = t96 * t38 - t151;
t35 = t170 * pkin(3) + qJDD(4) + t90 + t91;
t56 = pkin(4) * t96 + qJ(5) * t95 + pkin(3);
t58 = -t95 * pkin(4) + qJ(5) * t96;
t119 = t56 * t97 - t58 * t99 + qJ(2);
t3 = -qJDD(3) * pkin(4) + qJDD(5) - t4;
t46 = -t96 * t145 + t95 * t146;
t47 = -t95 * t145 - t96 * t146;
t51 = -t95 * t97 + t96 * t99;
t117 = t164 * t46 + t52 * t28 + t51 * t29 - t47 * t48;
t116 = 0.2e1 * qJ(2) * t134 + qJDD(3) * t101;
t115 = -t165 - t143;
t112 = -qJD(1) * t164 + t47 * qJD(3) + t51 * qJDD(3);
t111 = qJD(1) * t48 - t46 * qJD(3) + t52 * qJDD(3);
t24 = t95 * t41 + t36;
t74 = g(3) * t81;
t110 = t24 * qJD(3) - t165 * t82 + t4 + t74;
t109 = (t48 - t131) * qJD(3) + t132;
t108 = -t28 * pkin(4) + t29 * qJ(5) + t35;
t12 = -qJD(3) * pkin(4) + qJD(5) - t18;
t13 = qJD(3) * qJ(5) + t19;
t89 = qJDD(3) * qJ(5);
t2 = qJD(3) * qJD(5) + t5 + t89;
t107 = -t12 * t47 - t13 * t46 + t2 * t52 - t3 * t51 - t165;
t106 = t18 * t47 - t19 * t46 + t4 * t51 + t5 * t52 - t165;
t105 = -t21 * qJD(3) - t31 * qJDD(3) - t123 * t81;
t104 = -t101 * t102 + t169;
t83 = qJDD(3) * t99;
t77 = -t96 * pkin(3) - pkin(4);
t73 = t95 * pkin(3) + qJ(5);
t72 = t163 * t100;
t27 = t52 * pkin(4) - t51 * qJ(5) + t76;
t26 = pkin(3) * t147 + t48 * pkin(4) + qJ(5) * t164;
t9 = t122 + 0.2e1 * t167;
t8 = -t46 * pkin(4) - t47 * qJ(5) - t51 * qJD(5) + t66;
t1 = -t48 * qJD(5) + t108;
t6 = [qJDD(1), t165, t123, -0.2e1 * t144 + t171, t169, t168 * pkin(1) + (t114 + t90) * qJ(2), t94 * qJDD(1) - 0.2e1 * t97 * t129, 0.2e1 * t149 * t134 - 0.2e1 * t97 * t135, -t102 * t97 + t83, -t102 * t99 - t137, 0, t104 * t97 + t116 * t99, t104 * t99 - t116 * t97, t164 * t66 - t76 * t28 + t35 * t52 - t55 * t46 + t105, -t76 * t29 + t35 * t51 + t55 * t47 + t66 * t48 - t161, -t106 + t162, t5 * t32 + t19 * t22 - t4 * t31 - t18 * t21 + t35 * t76 + t55 * t66 - g(1) * (t76 * t100 - t152) - g(2) * (t76 * t98 + t72), t1 * t52 + t164 * t8 - t20 * t46 - t27 * t28 + t105, -t107 + t162, -t1 * t51 - t20 * t47 + t27 * t29 - t8 * t48 + t161, t2 * t32 + t13 * t22 + t1 * t27 + t20 * t8 + t3 * t31 + t12 * t21 - g(1) * (t119 * t100 - t152) - g(2) * (t119 * t98 + t72); 0, 0, 0, qJDD(1), -t103, -t143 - t168, 0, 0, 0, 0, 0, t139 * t97 + t83, t139 * t99 - t137, t112, -t111, t117, t106 - t142, t112, t117, t111, -t20 * qJD(1) + t107; 0, 0, 0, 0, 0, 0, t99 * t103 * t97, -t149 * t103, t135, -t136, qJDD(3), t115 * t99 + t157 + t53, g(3) * t99 + (-t115 - t61) * t97, -t55 * t48 + (qJDD(3) * t96 - t147 * t164) * pkin(3) + t110, t25 * qJD(3) + t55 * t164 + (-qJDD(3) * t95 - t48 * t147) * pkin(3) + t160, (t19 - t24) * t48 + (-t18 + t25) * t164 + (t28 * t95 + t29 * t96) * pkin(3), t18 * t24 - t19 * t25 + (t157 + t4 * t96 + t5 * t95 + (-t165 - t142) * t99) * pkin(3), -t155 - t26 * t164 - qJDD(5) + (pkin(4) - t77) * qJDD(3) + t110, t73 * t28 - t77 * t29 + (t13 - t24) * t48 + (t12 - t140) * t164, t73 * qJDD(3) - t20 * t164 + t26 * t48 + t89 + (0.2e1 * qJD(5) - t25) * qJD(3) - t160, t2 * t73 + t3 * t77 - t20 * t26 - t12 * t24 - g(3) * (-t81 * pkin(4) + t82 * qJ(5) - t156) + t140 * t13 - t165 * (t56 * t99 + t58 * t97); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t9, t166, t164 * t19 + t18 * t48 - t123 + t35, t109, t166, t9, t13 * t164 + (-qJD(5) - t12) * t48 + t108 - t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164 * t48 - qJDD(3), -t122, -t43 - t102, -t13 * qJD(3) + t165 * t51 + t155 + t3 - t74;];
tau_reg = t6;
