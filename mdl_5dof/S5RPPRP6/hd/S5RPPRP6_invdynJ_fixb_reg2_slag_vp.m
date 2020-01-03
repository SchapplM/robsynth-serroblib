% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:21
% EndTime: 2019-12-31 17:55:23
% DurationCPUTime: 0.91s
% Computational Cost: add. (1501->214), mult. (2818->234), div. (0->0), fcn. (1861->8), ass. (0->120)
t96 = sin(qJ(1));
t98 = cos(qJ(1));
t171 = g(1) * t96 - g(2) * t98;
t91 = sin(pkin(7));
t144 = qJD(1) * t91;
t95 = sin(qJ(4));
t133 = t95 * t144;
t92 = cos(pkin(7));
t97 = cos(qJ(4));
t54 = t97 * t91 + t95 * t92;
t109 = -qJD(4) * t133 + t54 * qJDD(1);
t151 = t97 * t92;
t132 = qJD(1) * t151;
t46 = t132 - t133;
t172 = t109 + (-t46 + t132) * qJD(4);
t85 = t91 ^ 2;
t86 = t92 ^ 2;
t147 = t85 + t86;
t88 = qJDD(1) * qJ(2);
t89 = qJD(1) * qJD(2);
t170 = t88 + t89;
t94 = -pkin(1) - qJ(3);
t63 = t94 * qJD(1) + qJD(2);
t169 = t147 * t63;
t124 = g(1) * t98 + g(2) * t96;
t60 = qJDD(3) + t170;
t110 = -t124 + t60;
t48 = t54 * qJD(4);
t55 = -t95 * t91 + t151;
t24 = -t48 * qJD(4) + t55 * qJDD(4);
t44 = t54 * qJD(1);
t168 = -qJD(1) * t44 + t24;
t23 = qJD(4) * t132 + t109;
t142 = qJD(4) * t97;
t143 = qJD(4) * t95;
t49 = t92 * t142 - t91 * t143;
t117 = t23 * t54 + t44 * t49;
t160 = t46 ^ 2;
t41 = t44 ^ 2;
t167 = -t41 - t160;
t166 = -t41 + t160;
t158 = -pkin(6) + t94;
t57 = t158 * t91;
t58 = t158 * t92;
t26 = t95 * t57 - t97 * t58;
t13 = -t54 * qJD(3) - t26 * qJD(4);
t27 = t97 * t57 + t95 * t58;
t14 = t55 * qJD(3) + t27 * qJD(4);
t137 = t92 * qJDD(1);
t138 = t91 * qJDD(1);
t120 = -t97 * t137 + t95 * t138;
t22 = qJD(1) * t48 + t120;
t164 = -t13 * t44 + t14 * t46 - t26 * t22 - t27 * t23;
t87 = pkin(7) + qJ(4);
t78 = cos(t87);
t163 = -t13 * qJD(4) - t27 * qJDD(4) - t124 * t78;
t161 = -qJD(1) * qJD(3) + qJDD(1) * t94;
t56 = qJDD(2) + t161;
t126 = -pkin(6) * qJDD(1) + t56;
t35 = t126 * t91;
t36 = t126 * t92;
t128 = -pkin(6) * qJD(1) + t63;
t38 = t128 * t92;
t134 = t38 * t142 + t97 * t35 + t95 * t36;
t77 = sin(t87);
t162 = -g(3) * t78 - t171 * t77 + t134;
t159 = 0.2e1 * t89;
t79 = t91 * pkin(3);
t154 = t46 * t44;
t37 = t128 * t91;
t153 = t95 * t37;
t150 = t49 * qJD(4) + t54 * qJDD(4);
t149 = t98 * pkin(1) + t96 * qJ(2);
t70 = qJ(2) + t79;
t146 = pkin(1) * qJDD(1);
t141 = qJDD(4) * pkin(4);
t19 = t97 * t37 + t95 * t38;
t140 = t19 * qJD(4);
t76 = qJD(1) * qJ(2) + qJD(3);
t18 = t97 * t38 - t153;
t139 = qJD(5) - t18;
t135 = qJDD(4) * qJ(5);
t59 = pkin(3) * t144 + t76;
t81 = t98 * qJ(2);
t131 = -t96 * pkin(1) + t81;
t130 = t147 * t56;
t129 = t18 + t153;
t127 = t37 * t142 + t38 * t143 + t95 * t35 - t97 * t36;
t73 = pkin(3) * t138;
t52 = t73 + t60;
t125 = qJDD(2) - t146;
t121 = qJD(1) * t46 + t150;
t118 = t77 * pkin(4) - t78 * qJ(5);
t6 = -t55 * t22 - t48 * t46;
t93 = -pkin(6) - qJ(3);
t116 = t98 * t79 + t96 * t93 + t131;
t114 = t96 * t79 - t98 * t93 + t149;
t112 = -t6 - t117;
t108 = g(3) * t77 - t171 * t78 - t127;
t107 = t23 * pkin(4) + t22 * qJ(5) + t52;
t106 = t22 * t54 - t55 * t23 + t48 * t44 - t46 * t49;
t1 = t135 + (qJD(5) - t153) * qJD(4) + t134;
t12 = -qJD(4) * pkin(4) + t139;
t16 = qJD(4) * qJ(5) + t19;
t3 = qJDD(5) + t127 - t141;
t105 = t1 * t54 + t12 * t48 + t16 * t49 - t3 * t55 - t171;
t4 = -t37 * t143 + t134;
t104 = -t127 * t55 - t18 * t48 + t19 * t49 + t4 * t54 - t171;
t103 = t110 + t170;
t102 = -t14 * qJD(4) - t26 * qJDD(4) - t124 * t77;
t17 = t44 * pkin(4) - t46 * qJ(5) + t59;
t101 = t17 * t46 + qJDD(5) - t108;
t100 = (t46 + t132) * qJD(4) + t109;
t99 = qJD(1) ^ 2;
t21 = t46 * pkin(4) + t44 * qJ(5);
t20 = t54 * pkin(4) - t55 * qJ(5) + t70;
t10 = 0.2e1 * t44 * qJD(4) + t120;
t9 = t49 * pkin(4) + t48 * qJ(5) - t55 * qJD(5) + qJD(2);
t2 = -t46 * qJD(5) + t107;
t5 = [0, 0, 0, 0, 0, qJDD(1), t171, t124, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - 0.2e1 * t146 - t171, -t124 + 0.2e1 * t88 + t159, -t125 * pkin(1) - g(1) * t131 - g(2) * t149 + (t88 + t159) * qJ(2), t86 * qJDD(1), -0.2e1 * t91 * t137, 0, t85 * qJDD(1), 0, 0, t103 * t91, t103 * t92, t171 + t147 * (-t161 - t56), t60 * qJ(2) + t76 * qJD(2) - g(1) * (t94 * t96 + t81) - g(2) * (t98 * qJ(3) + t149) + t94 * t130 - qJD(3) * t169, t6, t106, t24, t117, -t150, 0, qJD(2) * t44 + t70 * t23 + t59 * t49 + t52 * t54 + t102, qJD(2) * t46 - t70 * t22 - t59 * t48 + t52 * t55 + t163, -t104 + t164, -g(1) * t116 - g(2) * t114 + t59 * qJD(2) + t127 * t26 + t19 * t13 - t18 * t14 + t4 * t27 + t52 * t70, t6, t24, -t106, 0, t150, t117, t17 * t49 + t2 * t54 + t20 * t23 + t9 * t44 + t102, -t105 + t164, t17 * t48 - t2 * t55 + t20 * t22 - t9 * t46 - t163, t1 * t27 + t16 * t13 + t2 * t20 + t17 * t9 + t3 * t26 + t12 * t14 - g(1) * (t118 * t98 + t116) - g(2) * (t118 * t96 + t114); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t99, -t99 * qJ(2) + t125 - t171, 0, 0, 0, 0, 0, 0, -t99 * t91, -t99 * t92, -t147 * qJDD(1), -t76 * qJD(1) + t130 - t171, 0, 0, 0, 0, 0, 0, t168, -t121, t112, -t59 * qJD(1) + t104, 0, 0, 0, 0, 0, 0, t168, t112, t121, -t17 * qJD(1) + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, t137, -t147 * t99, qJD(1) * t169 + t110, 0, 0, 0, 0, 0, 0, t100, -t10, t167, t18 * t46 + t19 * t44 + t110 + t73, 0, 0, 0, 0, 0, 0, t100, t167, t10, t16 * t44 + (-qJD(5) - t12) * t46 + t107 - t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t166, -t120, -t154, -t172, qJDD(4), -t59 * t46 + t108 + t140, t129 * qJD(4) + t59 * t44 - t162, 0, 0, t154, -t120, -t166, qJDD(4), t172, -t154, -t21 * t44 - t101 + t140 + 0.2e1 * t141, pkin(4) * t22 - t23 * qJ(5) + (t16 - t19) * t46 + (t12 - t139) * t44, 0.2e1 * t135 - t17 * t44 + t21 * t46 + (0.2e1 * qJD(5) - t129) * qJD(4) + t162, -t3 * pkin(4) + g(3) * t118 + t1 * qJ(5) - t12 * t19 + t139 * t16 - t17 * t21 - t171 * (pkin(4) * t78 + qJ(5) * t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t154, -t120, -qJD(4) ^ 2 - t160, -t16 * qJD(4) + t101 - t141;];
tau_reg = t5;
