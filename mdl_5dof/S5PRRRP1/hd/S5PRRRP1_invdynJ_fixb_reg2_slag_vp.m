% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:22
% EndTime: 2019-12-05 16:40:26
% DurationCPUTime: 1.18s
% Computational Cost: add. (1239->234), mult. (1822->280), div. (0->0), fcn. (944->8), ass. (0->142)
t89 = sin(qJ(3));
t135 = qJDD(2) * t89;
t91 = cos(qJ(3));
t142 = qJD(3) * t91;
t81 = qJDD(2) + qJDD(3);
t23 = t81 * pkin(7) + (qJD(2) * t142 + t135) * pkin(2);
t169 = qJD(1) * qJD(4) + t23;
t147 = pkin(2) * qJD(2);
t128 = t89 * t147;
t83 = qJD(2) + qJD(3);
t37 = t83 * pkin(7) + t128;
t90 = cos(qJ(4));
t77 = t90 * qJD(1);
t88 = sin(qJ(4));
t20 = -t88 * t37 + t77;
t137 = t88 * qJD(1);
t21 = t90 * t37 + t137;
t140 = qJD(4) * t88;
t5 = t88 * qJDD(1) - t37 * t140 + t169 * t90;
t74 = t90 * qJDD(1);
t6 = -t21 * qJD(4) - t88 * t23 + t74;
t95 = t5 * t90 - t6 * t88 + (-t20 * t90 - t21 * t88) * qJD(4);
t82 = pkin(8) + qJ(2);
t75 = qJ(3) + t82;
t64 = sin(t75);
t57 = g(2) * t64;
t65 = cos(t75);
t59 = g(1) * t65;
t151 = t57 + t59;
t119 = qJ(5) * t83 + t37;
t106 = t119 * t90;
t125 = t83 * t140;
t68 = t90 * pkin(4) + pkin(3);
t156 = t68 * t81;
t168 = -pkin(4) * t125 + t156;
t109 = t20 * t88 - t21 * t90;
t144 = qJD(2) * t91;
t127 = pkin(2) * t144;
t38 = -t83 * pkin(3) - t127;
t167 = t109 * t91 - t38 * t89;
t84 = t88 ^ 2;
t166 = pkin(4) * t84;
t58 = g(1) * t64;
t71 = sin(t82);
t165 = g(1) * t71;
t164 = g(2) * t65;
t163 = g(3) * t90;
t162 = t81 * pkin(3);
t161 = t91 * pkin(2);
t13 = -t119 * t88 + t77;
t146 = qJD(4) * pkin(4);
t9 = t13 + t146;
t160 = -t13 + t9;
t158 = t64 * t90;
t157 = t65 * t88;
t80 = t83 ^ 2;
t155 = t80 * t90;
t154 = t83 * t88;
t153 = t83 * t90;
t60 = t88 * t81;
t61 = t90 * t81;
t87 = -qJ(5) - pkin(7);
t152 = t65 * pkin(3) + t64 * pkin(7);
t150 = -qJD(3) * t128 + qJDD(2) * t161;
t85 = t90 ^ 2;
t149 = t84 - t85;
t148 = t84 + t85;
t78 = t90 * qJ(5);
t67 = t89 * pkin(2) + pkin(7);
t145 = -qJ(5) - t67;
t143 = qJD(3) * t89;
t141 = qJD(4) * t83;
t139 = qJD(4) * t90;
t138 = qJDD(4) * pkin(4);
t76 = t90 * qJD(5);
t19 = -t68 * t83 + qJD(5) - t127;
t136 = -qJD(5) - t19;
t3 = t83 * t76 + (-t125 + t61) * qJ(5) + t5;
t133 = t3 * t90 - t151;
t43 = g(2) * t157;
t8 = qJDD(5) - t150 - t168;
t132 = t19 * t139 + t8 * t88 + t43;
t22 = -t150 - t162;
t131 = t38 * t139 + t22 * t88 + t43;
t114 = qJD(4) * t127;
t116 = t83 * t128;
t44 = g(1) * t158;
t130 = t88 * t114 + t90 * t116 + t44;
t129 = pkin(2) * t142;
t126 = t83 * t143;
t15 = t19 * t140;
t124 = t83 * t139;
t123 = -t8 - t164;
t122 = -t22 - t164;
t121 = t148 * t81;
t120 = -t64 * t87 + t65 * t68;
t118 = qJD(4) * t87;
t117 = qJD(4) * t145;
t115 = t88 * t124;
t113 = g(1) * (-t64 * pkin(3) + t65 * pkin(7));
t92 = qJD(4) ^ 2;
t112 = -pkin(7) * t92 + t162;
t72 = cos(t82);
t111 = -g(2) * t72 + t165;
t14 = t137 + t106;
t110 = t14 * t90 - t88 * t9;
t35 = pkin(2) * t143 + pkin(4) * t140;
t47 = -t68 - t161;
t108 = t35 * t83 + t47 * t81;
t107 = -t64 * t68 - t65 * t87;
t105 = -t150 - t58 + t164;
t104 = g(1) * t157 + t88 * t57 - t163 + t74;
t103 = t148 * t127;
t102 = -pkin(3) * t141 - pkin(7) * qJDD(4);
t101 = g(2) * t158 + g(3) * t88 + t90 * t59 - t5;
t100 = -qJ(5) * t81 - t169;
t99 = -t116 - t58;
t69 = -pkin(3) - t161;
t98 = pkin(2) * t126 + t67 * t92 + t69 * t81;
t96 = -qJDD(4) * t67 + (t69 * t83 - t129) * qJD(4);
t94 = -t151 + t95;
t86 = qJDD(1) - g(3);
t63 = pkin(2) * t72;
t52 = t88 * t155;
t51 = t90 * pkin(7) + t78;
t50 = t87 * t88;
t49 = t90 * t114;
t40 = qJDD(4) * t90 - t92 * t88;
t39 = qJDD(4) * t88 + t92 * t90;
t33 = t90 * t67 + t78;
t32 = t149 * t80;
t31 = t145 * t88;
t28 = t38 * t140;
t27 = -t88 * qJD(5) + t90 * t118;
t26 = t88 * t118 + t76;
t25 = t85 * t81 - 0.2e1 * t115;
t24 = t84 * t81 + 0.2e1 * t115;
t12 = -0.2e1 * t149 * t141 + 0.2e1 * t88 * t61;
t11 = (-qJD(5) - t129) * t88 + t90 * t117;
t10 = t88 * t117 + t90 * t129 + t76;
t2 = t138 + t74 - qJD(4) * t106 + (-qJD(5) * t83 + t100) * t88;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, t40, -t39, 0, -t109 * qJD(4) + t5 * t88 + t6 * t90 - g(3), 0, 0, 0, 0, 0, 0, t40, -t39, 0, qJD(4) * t110 + t2 * t90 + t3 * t88 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t111, g(1) * t72 + g(2) * t71, 0, 0, 0, 0, 0, 0, 0, t81, (t81 * t91 - t126) * pkin(2) - t105, ((-qJDD(2) - t81) * t89 + (-qJD(2) - t83) * t142) * pkin(2) + t151, 0, (t111 + (t89 ^ 2 + t91 ^ 2) * qJDD(2) * pkin(2)) * pkin(2), t24, t12, t39, t25, t40, 0, t28 + t44 + t96 * t88 + (t122 - t98) * t90, t96 * t90 + (t98 - t58) * t88 + t131, t148 * t83 * t129 + t67 * t121 + t94, t22 * t69 - t113 - g(2) * (t63 + t152) + (-t167 * qJD(3) + t165) * pkin(2) + t95 * t67, t24, t12, t39, t25, t40, 0, t31 * qJDD(4) + t15 + t44 + (t47 * t154 + t11) * qJD(4) + (-t108 + t123) * t90, -t33 * qJDD(4) + (t47 * t153 - t10) * qJD(4) + (t108 - t58) * t88 + t132, (t10 * t83 + t33 * t81 + (-t31 * t83 - t9) * qJD(4)) * t90 + (-t11 * t83 - t31 * t81 - t2 + (-t33 * t83 - t14) * qJD(4)) * t88 + t133, t3 * t33 + t14 * t10 + t2 * t31 + t9 * t11 + t8 * t47 + t19 * t35 - g(1) * (-pkin(2) * t71 + t107) - g(2) * (t120 + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t105 + t116, (-t135 + (-qJD(3) + t83) * t144) * pkin(2) + t151, 0, 0, t24, t12, t39, t25, t40, 0, t28 + t102 * t88 + (t112 + t122) * t90 + t130, t49 + t102 * t90 + (-t112 + t99) * t88 + t131, pkin(7) * t121 - t83 * t103 + t94, -t22 * pkin(3) + t95 * pkin(7) - g(2) * t152 + t167 * t147 - t113, t24, t12, t39, t25, t40, 0, t50 * qJDD(4) + t15 + (-t68 * t154 + t27) * qJD(4) + (t123 + t168) * t90 + t130, -t51 * qJDD(4) + t49 + (t99 - t156) * t88 + (-t26 + (-t68 * t90 + t166) * t83) * qJD(4) + t132, (-qJD(4) * t9 + t51 * t81) * t90 + (-t14 * qJD(4) - t50 * t81 - t2) * t88 + (t26 * t90 - t27 * t88 + (-t50 * t90 - t51 * t88) * qJD(4) - t103) * t83 + t133, t3 * t51 + t14 * t26 + t2 * t50 + t9 * t27 - t8 * t68 + pkin(4) * t15 - g(1) * t107 - g(2) * t120 + (-t110 * t91 - t19 * t89) * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t32, t60, t52, t61, qJDD(4), (-t38 * t83 - t23) * t88 + t104, t20 * qJD(4) - t38 * t153 + t101, 0, 0, -t52, t32, t60, t52, t61, qJDD(4), 0.2e1 * t138 + (t14 - t106) * qJD(4) + (pkin(4) * t155 + t136 * t83 + t100) * t88 + t104, -t80 * t166 - t81 * t78 + t13 * qJD(4) + (qJ(5) * t140 + t136 * t90) * t83 + t101, -pkin(4) * t60 + (-t146 + t160) * t153, t160 * t14 + (-t163 + t2 + (-t19 * t83 + t151) * t88) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 + 0.2e1 * t125, t60 + 0.2e1 * t124, -t148 * t80, -t110 * t83 - t123 - t58;];
tau_reg = t1;
