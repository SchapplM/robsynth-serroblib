% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR1
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
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:42
% EndTime: 2022-01-20 09:51:45
% DurationCPUTime: 0.74s
% Computational Cost: add. (1043->175), mult. (1683->226), div. (0->0), fcn. (1120->16), ass. (0->124)
t113 = cos(qJ(2));
t153 = pkin(1) * qJD(2);
t139 = qJD(1) * t153;
t110 = sin(qJ(2));
t144 = qJDD(1) * t110;
t168 = pkin(1) * t144 + t113 * t139;
t107 = cos(pkin(9));
t112 = cos(qJ(5));
t147 = t112 * t107;
t105 = sin(pkin(9));
t109 = sin(qJ(5));
t150 = t109 * t105;
t51 = -t147 + t150;
t52 = t112 * t105 + t109 * t107;
t104 = qJ(1) + qJ(2);
t93 = pkin(8) + t104;
t79 = cos(t93);
t163 = g(2) * t79;
t78 = sin(t93);
t164 = g(1) * t78;
t167 = t163 - t164;
t103 = qJD(1) + qJD(2);
t38 = t52 * t103;
t94 = sin(t104);
t95 = cos(t104);
t166 = g(1) * t94 - g(2) * t95;
t108 = cos(pkin(8));
t154 = pkin(1) * qJD(1);
t140 = t113 * t154;
t106 = sin(pkin(8));
t141 = t110 * t154;
t68 = t106 * t141;
t45 = t108 * t140 - t68;
t146 = qJD(4) - t45;
t145 = t105 ^ 2 + t107 ^ 2;
t165 = pkin(2) * t94;
t161 = pkin(1) * t110;
t160 = t107 * pkin(4);
t159 = t108 * pkin(2);
t111 = sin(qJ(1));
t158 = t111 * pkin(1);
t157 = t113 * pkin(1);
t86 = qJDD(1) * t157;
t99 = qJDD(1) + qJDD(2);
t40 = t99 * pkin(2) - t110 * t139 + t86;
t21 = t106 * t40 + t168 * t108;
t14 = t99 * qJ(4) + t103 * qJD(4) + t21;
t9 = t105 * qJDD(3) + t107 * t14;
t55 = t103 * pkin(2) + t140;
t69 = t108 * t141;
t33 = t106 * t55 + t69;
t151 = t108 * t110;
t85 = pkin(2) + t157;
t156 = pkin(1) * t151 + t106 * t85;
t155 = g(1) * t95 + g(2) * t94;
t152 = t107 * t99;
t136 = t103 * t147;
t143 = qJD(5) * t136 + t52 * t99;
t84 = pkin(2) * t95;
t142 = t79 * pkin(3) + t78 * qJ(4) + t84;
t74 = t106 * t161;
t138 = t103 * t150;
t135 = -pkin(3) - t160;
t20 = -t168 * t106 + t108 * t40;
t118 = qJDD(4) - t20;
t15 = -t99 * pkin(3) + t118;
t134 = -t15 - t163;
t32 = t108 * t55 - t68;
t133 = t108 * t85 - t74;
t132 = qJD(1) * (-qJD(2) + t103);
t131 = qJD(2) * (-qJD(1) - t103);
t130 = t86 + t166;
t42 = -pkin(3) - t133;
t128 = g(1) * t79 + g(2) * t78;
t127 = t51 * t99;
t126 = qJD(4) - t32;
t41 = qJ(4) + t156;
t30 = (-pkin(7) - t41) * t105;
t96 = t107 * pkin(7);
t31 = t107 * t41 + t96;
t124 = t109 * t31 - t112 * t30;
t123 = t109 * t30 + t112 * t31;
t76 = t106 * pkin(2) + qJ(4);
t48 = (-pkin(7) - t76) * t105;
t49 = t107 * t76 + t96;
t122 = t109 * t49 - t112 * t48;
t121 = t109 * t48 + t112 * t49;
t120 = -t78 * pkin(3) + t79 * qJ(4) - t165;
t119 = t108 * t113 * t153 - qJD(2) * t74;
t10 = t135 * t99 + t118;
t24 = t135 * t103 + t126;
t46 = t51 * qJD(5);
t102 = pkin(9) + qJ(5);
t91 = sin(t102);
t117 = t10 * t52 + t167 * t91 - t24 * t46;
t47 = t52 * qJD(5);
t92 = cos(t102);
t116 = t10 * t51 - t167 * t92 + t24 * t47;
t88 = t107 * qJDD(3);
t8 = -t105 * t14 + t88;
t115 = -t8 * t105 + t9 * t107 - t128;
t114 = cos(qJ(1));
t97 = t114 * pkin(1);
t80 = -pkin(3) - t159;
t58 = t107 * t164;
t57 = t135 - t159;
t44 = (t106 * t113 + t151) * t153;
t43 = t106 * t140 + t69;
t39 = qJD(4) + t119;
t36 = -t136 + t138;
t34 = t42 - t160;
t29 = t103 * qJ(4) + t33;
t28 = -t103 * pkin(3) + t126;
t26 = -t47 * qJD(5) - t51 * qJDD(5);
t25 = -t46 * qJD(5) + t52 * qJDD(5);
t23 = t105 * qJD(3) + t107 * t29;
t22 = t107 * qJD(3) - t105 * t29;
t19 = t103 * t47 + t127;
t18 = -qJD(5) * t138 + t143;
t4 = pkin(7) * t152 + t9;
t3 = t88 + (-pkin(7) * t99 - t14) * t105;
t2 = t18 * t52 - t38 * t46;
t1 = -t18 * t51 - t52 * t19 + t46 * t36 - t38 * t47;
t5 = [qJDD(1), g(1) * t111 - g(2) * t114, g(1) * t114 + g(2) * t111, t99, (t110 * t131 + t113 * t99) * pkin(1) + t130, ((-qJDD(1) - t99) * t110 + t113 * t131) * pkin(1) + t155, t21 * t156 + t33 * t119 + t20 * t133 - t32 * t44 - g(1) * (-t158 - t165) - g(2) * (t84 + t97), t58 + (-t103 * t44 - t42 * t99 + t134) * t107, t115 + t145 * (t103 * t39 + t41 * t99), t15 * t42 + t28 * t44 - g(1) * (t120 - t158) - g(2) * (t97 + t142) + (t23 * t39 + t9 * t41) * t107 + (-t22 * t39 - t8 * t41) * t105, t2, t1, t25, t26, 0, t44 * t36 + t34 * t19 - t124 * qJDD(5) + (-t123 * qJD(5) - t52 * t39) * qJD(5) + t116, t44 * t38 + t34 * t18 - t123 * qJDD(5) + (t124 * qJD(5) + t51 * t39) * qJD(5) + t117; 0, 0, 0, t99, t132 * t161 + t130, (t113 * t132 - t144) * pkin(1) + t155, t32 * t43 - t33 * t45 + (t106 * t21 + t108 * t20 + t166) * pkin(2), t58 + (t103 * t43 - t80 * t99 + t134) * t107, t115 + (t146 * t103 + t99 * t76) * t145, t15 * t80 - t28 * t43 - g(1) * t120 - g(2) * t142 + (t146 * t23 + t9 * t76) * t107 + (-t146 * t22 - t8 * t76) * t105, t2, t1, t25, t26, 0, t57 * t19 - t122 * qJDD(5) - t43 * t36 + (-t121 * qJD(5) - t146 * t52) * qJD(5) + t116, t57 * t18 - t121 * qJDD(5) - t43 * t38 + (t122 * qJD(5) + t146 * t51) * qJD(5) + t117; 0, 0, 0, 0, 0, 0, qJDD(3) - g(3), 0, 0, t9 * t105 + t8 * t107 - g(3), 0, 0, 0, 0, 0, t26, -t25; 0, 0, 0, 0, 0, 0, 0, -t152, -t145 * t103 ^ 2, -t164 + (t22 * t105 - t23 * t107) * t103 - t134, 0, 0, 0, 0, 0, 0.2e1 * t38 * qJD(5) + t127, (-t36 - t138) * qJD(5) + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t36, -t36 ^ 2 + t38 ^ 2, (t36 - t138) * qJD(5) + t143, -t127, qJDD(5), -g(3) * t92 - t109 * t4 + t112 * t3 + t128 * t91 - t24 * t38, g(3) * t91 - t109 * t3 - t112 * t4 + t128 * t92 + t24 * t36;];
tau_reg = t5;
