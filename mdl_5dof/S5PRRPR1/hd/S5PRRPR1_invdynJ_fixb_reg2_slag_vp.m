% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPR1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:57
% EndTime: 2019-12-05 16:15:59
% DurationCPUTime: 1.06s
% Computational Cost: add. (1690->197), mult. (2438->244), div. (0->0), fcn. (1629->12), ass. (0->131)
t118 = sin(qJ(3));
t163 = pkin(2) * qJD(2);
t149 = t118 * t163;
t120 = cos(qJ(3));
t172 = t120 * pkin(2);
t164 = -qJD(3) * t149 + qJDD(2) * t172;
t146 = qJDD(4) - t164;
t107 = qJDD(2) + qJDD(3);
t173 = t107 * pkin(3);
t48 = t146 - t173;
t111 = pkin(8) + qJ(2);
t104 = qJ(3) + t111;
t91 = cos(t104);
t82 = g(2) * t91;
t180 = t48 + t82;
t90 = sin(t104);
t83 = g(1) * t90;
t178 = t82 - t83;
t115 = cos(pkin(9));
t119 = cos(qJ(5));
t157 = t119 * t115;
t114 = sin(pkin(9));
t117 = sin(qJ(5));
t160 = t117 * t114;
t58 = -t157 + t160;
t59 = t119 * t114 + t117 * t115;
t155 = qJD(2) * t120;
t132 = -pkin(2) * t155 + qJD(4);
t179 = g(1) * t91 + g(2) * t90;
t112 = qJD(2) + qJD(3);
t51 = t59 * t112;
t152 = qJDD(2) * t118;
t153 = qJD(3) * t120;
t40 = t107 * qJ(4) + t112 * qJD(4) + (qJD(2) * t153 + t152) * pkin(2);
t97 = t115 * qJDD(1);
t30 = -t114 * t40 + t97;
t31 = t114 * qJDD(1) + t115 * t40;
t177 = -t30 * t114 + t31 * t115;
t108 = t114 ^ 2;
t109 = t115 ^ 2;
t156 = t108 + t109;
t176 = t112 * t156;
t175 = t51 ^ 2;
t101 = sin(t111);
t174 = pkin(2) * t101;
t105 = t115 * pkin(7);
t147 = t112 * t160;
t49 = -t112 * t157 + t147;
t171 = t51 * t49;
t116 = -pkin(7) - qJ(4);
t70 = t116 * t114;
t71 = t115 * qJ(4) + t105;
t38 = -t117 * t71 + t119 * t70;
t170 = t38 * qJD(5) - t132 * t58;
t39 = t117 * t70 + t119 * t71;
t169 = -t39 * qJD(5) - t132 * t59;
t134 = t58 * t107;
t54 = t59 * qJD(5);
t19 = t112 * t54 + t134;
t141 = qJD(5) * t157;
t142 = qJD(5) * t160;
t53 = -t141 + t142;
t168 = -t59 * t19 + t53 * t49;
t167 = t180 * t114;
t64 = t112 * qJ(4) + t149;
t45 = t114 * qJD(1) + t115 * t64;
t166 = t91 * pkin(3) + t90 * qJ(4);
t161 = t115 * t107;
t154 = qJD(3) * t118;
t151 = t59 * t107 + t112 * t141;
t150 = pkin(2) * t154;
t93 = t115 * pkin(4) + pkin(3);
t145 = -t90 * pkin(3) + t91 * qJ(4);
t143 = t112 * t154;
t140 = -t90 * t116 + t91 * t93;
t20 = t97 + (-pkin(7) * t107 - t40) * t114;
t21 = pkin(7) * t161 + t31;
t139 = -t117 * t21 + t119 * t20;
t138 = -t164 + t178;
t137 = t156 * t107;
t136 = t112 * t149;
t99 = t115 * qJD(1);
t34 = t99 + (-pkin(7) * t112 - t64) * t114;
t35 = t112 * t105 + t45;
t12 = -t117 * t35 + t119 * t34;
t13 = t117 * t34 + t119 * t35;
t127 = t117 * t20 + t119 * t21;
t4 = t12 * qJD(5) + t127;
t5 = -t13 * qJD(5) + t139;
t133 = t12 * t53 - t13 * t54 - t4 * t58 - t5 * t59 - t179;
t103 = cos(t111);
t131 = g(1) * t101 - g(2) * t103;
t18 = t112 * t142 - t151;
t130 = -t58 * t18 + t51 * t54;
t129 = -t91 * t116 - t90 * t93;
t44 = -t114 * t64 + t99;
t128 = t44 * t114 - t45 * t115;
t92 = t118 * pkin(2) + qJ(4);
t56 = (-pkin(7) - t92) * t114;
t57 = t115 * t92 + t105;
t27 = -t117 * t57 + t119 * t56;
t28 = t117 * t56 + t119 * t57;
t126 = -t179 + t177;
t110 = pkin(9) + qJ(5);
t100 = sin(t110);
t36 = -t93 * t107 + t146;
t47 = -t93 * t112 + t132;
t125 = t178 * t100 + t36 * t59 - t47 * t53;
t102 = cos(t110);
t124 = -t178 * t102 + t36 * t58 + t47 * t54;
t123 = -t136 - t173;
t94 = -pkin(3) - t172;
t122 = pkin(2) * t143 + t107 * t94;
t113 = qJDD(1) - g(3);
t88 = pkin(2) * t103;
t86 = t109 * t107;
t85 = t108 * t107;
t78 = pkin(2) * t153 + qJD(4);
t69 = -t93 - t172;
t68 = t115 * t83;
t65 = 0.2e1 * t114 * t161;
t61 = -t112 * pkin(3) + t132;
t46 = t49 ^ 2;
t33 = -t54 * qJD(5) - t58 * qJDD(5);
t32 = -t53 * qJD(5) + t59 * qJDD(5);
t11 = -t28 * qJD(5) - t59 * t78;
t10 = t27 * qJD(5) - t58 * t78;
t7 = t19 * t58 + t49 * t54;
t6 = -t18 * t59 - t51 * t53;
t1 = -t130 + t168;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * t114 + t30 * t115 - g(3), 0, 0, 0, 0, 0, 0, t33, -t32, t130 + t168, -t12 * t54 - t13 * t53 + t4 * t59 - t5 * t58 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t131, g(1) * t103 + g(2) * t101, 0, 0, 0, 0, 0, 0, 0, t107, (t107 * t120 - t143) * pkin(2) - t138, ((-qJDD(2) - t107) * t118 + (-qJD(2) - t112) * t153) * pkin(2) + t179, 0, (t131 + (t118 ^ 2 + t120 ^ 2) * qJDD(2) * pkin(2)) * pkin(2), t85, t65, 0, t86, 0, 0, t68 + (-t122 - t180) * t115, (t122 - t83) * t114 + t167, t92 * t137 + t176 * t78 + t126, t48 * t94 + t61 * t150 - g(1) * (t145 - t174) - g(2) * (t88 + t166) + (t31 * t92 + t45 * t78) * t115 + (-t30 * t92 - t44 * t78) * t114, t6, t1, t32, t7, t33, 0, t11 * qJD(5) + t27 * qJDD(5) + t49 * t150 + t69 * t19 + t124, -t10 * qJD(5) - t28 * qJDD(5) + t51 * t150 - t69 * t18 + t125, -t10 * t49 - t11 * t51 + t27 * t18 - t28 * t19 + t133, t4 * t28 + t13 * t10 + t5 * t27 + t12 * t11 + t36 * t69 + t47 * t150 - g(1) * (t129 - t174) - g(2) * (t140 + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t136 - t138, (-t152 + (-qJD(3) + t112) * t155) * pkin(2) + t179, 0, 0, t85, t65, 0, t86, 0, 0, t68 + (-t123 - t180) * t115, (t123 - t83) * t114 + t167, qJ(4) * t137 + t132 * t176 + t126, -t48 * pkin(3) - g(1) * t145 - g(2) * t166 - t128 * qJD(4) + t177 * qJ(4) + (-t118 * t61 + t128 * t120) * t163, t6, t1, t32, t7, t33, 0, t169 * qJD(5) + t38 * qJDD(5) - t49 * t149 - t93 * t19 + t124, -t170 * qJD(5) - t39 * qJDD(5) - t51 * t149 + t93 * t18 + t125, -t169 * t51 - t170 * t49 + t38 * t18 - t39 * t19 + t133, -g(1) * t129 - g(2) * t140 + t169 * t12 + t170 * t13 - t47 * t149 - t36 * t93 + t5 * t38 + t4 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, t114 * t107, -t156 * t112 ^ 2, t128 * t112 + t178 + t48, 0, 0, 0, 0, 0, 0, 0.2e1 * t51 * qJD(5) + t134, (-t49 - t147) * qJD(5) + t151, -t46 - t175, t12 * t51 + t13 * t49 + t178 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, -t46 + t175, (t49 - t147) * qJD(5) + t151, -t171, -t134, qJDD(5), -g(3) * t102 + t100 * t179 - t47 * t51 + t139, g(3) * t100 + t102 * t179 + t47 * t49 - t127, 0, 0;];
tau_reg = t2;
