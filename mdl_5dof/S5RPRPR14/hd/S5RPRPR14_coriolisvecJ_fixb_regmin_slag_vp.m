% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR14_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:16
% EndTime: 2019-12-31 18:35:19
% DurationCPUTime: 0.95s
% Computational Cost: add. (1167->175), mult. (2585->260), div. (0->0), fcn. (1732->6), ass. (0->107)
t118 = cos(pkin(8));
t74 = sin(pkin(8));
t76 = sin(qJ(3));
t78 = cos(qJ(3));
t88 = t118 * t76 + t74 * t78;
t142 = t88 * qJD(1);
t144 = qJD(5) + t142;
t77 = cos(qJ(5));
t110 = t77 * qJD(3);
t117 = qJD(1) * t76;
t98 = t118 * t78;
t45 = qJD(1) * t98 - t74 * t117;
t75 = sin(qJ(5));
t26 = t75 * t45 - t110;
t147 = t144 * t26;
t79 = -pkin(1) - pkin(6);
t54 = t79 * qJD(1) + qJD(2);
t146 = -qJ(4) * qJD(1) + t54;
t36 = t45 * qJD(3);
t128 = t75 * t36;
t96 = t144 * t77;
t145 = -t144 * t96 - t128;
t141 = -qJD(5) + t144;
t104 = qJ(4) * qJD(3);
t111 = t76 * qJD(4);
t115 = qJD(3) * t78;
t25 = t54 * t115 + (-t78 * t104 - t111) * qJD(1);
t109 = t78 * qJD(4);
t116 = qJD(3) * t76;
t82 = -t54 * t116 + (t76 * t104 - t109) * qJD(1);
t3 = -t118 * t82 + t74 * t25;
t63 = t74 * pkin(3) + pkin(7);
t140 = (t78 * qJD(1) * pkin(3) + t45 * pkin(4) + pkin(7) * t142 + qJD(5) * t63) * t144 + t3;
t51 = pkin(3) * t117 + qJD(1) * qJ(2) + qJD(4);
t13 = pkin(4) * t142 - t45 * pkin(7) + t51;
t119 = qJ(4) - t79;
t52 = t119 * t76;
t99 = t119 * t78;
t24 = -t118 * t52 - t74 * t99;
t133 = t24 * t36;
t48 = -t74 * t76 + t98;
t137 = t3 * t48;
t38 = -qJD(3) * t99 - t111;
t86 = t119 * t116 - t109;
t15 = t118 * t38 + t74 * t86;
t120 = t76 * pkin(3) + qJ(2);
t21 = pkin(4) * t88 - t48 * pkin(7) + t120;
t4 = t118 * t25 + t74 * t82;
t93 = qJD(3) * t118;
t44 = -t74 * t115 - t76 * t93;
t39 = t146 * t76;
t129 = t74 * t39;
t40 = t146 * t78;
t35 = qJD(3) * pkin(3) + t40;
t11 = t118 * t35 - t129;
t7 = -qJD(3) * pkin(4) - t11;
t139 = -(qJD(5) * t21 + t15) * t144 + t7 * t44 - (qJD(5) * t13 + t4) * t88 - t133 + t137;
t70 = qJD(1) * qJD(2);
t138 = 0.2e1 * t70;
t113 = qJD(5) * t75;
t37 = qJD(3) * t142;
t9 = qJD(5) * t110 - t45 * t113 - t77 * t37;
t136 = t48 * t9;
t135 = t9 * t75;
t134 = t21 * t36;
t28 = t75 * qJD(3) + t77 * t45;
t132 = t28 * t45;
t131 = t45 * t26;
t130 = t88 * t36;
t127 = t75 * t37;
t31 = t77 * t36;
t80 = qJD(3) ^ 2;
t126 = t80 * t76;
t125 = t80 * t78;
t33 = t118 * t39;
t12 = t74 * t35 + t33;
t103 = qJD(1) * qJD(3);
t100 = t78 * t103;
t124 = pkin(3) * t100 + t70;
t123 = t76 ^ 2 - t78 ^ 2;
t81 = qJD(1) ^ 2;
t122 = -t80 - t81;
t121 = t81 * qJ(2);
t114 = qJD(5) * t48;
t112 = t51 * qJD(1);
t108 = pkin(3) * t115 + qJD(2);
t106 = qJ(2) * qJD(3);
t102 = 0.2e1 * qJD(1);
t92 = -qJD(5) * t88 - qJD(1);
t8 = qJD(3) * pkin(7) + t12;
t1 = t77 * t13 - t75 * t8;
t2 = t75 * t13 + t77 * t8;
t90 = t31 + (-t142 * t75 - t113) * t144;
t89 = -t48 * t113 + t77 * t44;
t17 = t118 * t40 - t129;
t84 = -t63 * t36 + (t17 + t7) * t144;
t43 = t74 * t116 - t78 * t93;
t83 = t11 * t44 - t12 * t43 + t4 * t88 - t137;
t64 = -t118 * pkin(3) - pkin(4);
t23 = t119 * t98 - t74 * t52;
t18 = -t43 * pkin(4) - t44 * pkin(7) + t108;
t16 = t74 * t40 + t33;
t14 = -t118 * t86 + t74 * t38;
t10 = t28 * qJD(5) - t127;
t6 = t36 * pkin(4) + t37 * pkin(7) + t124;
t5 = t77 * t6;
t19 = [0, 0, 0, 0, t138, qJ(2) * t138, -0.2e1 * t76 * t100, 0.2e1 * t123 * t103, -t126, -t125, 0, -t79 * t126 + (qJD(2) * t76 + t78 * t106) * t102, -t79 * t125 + (qJD(2) * t78 - t76 * t106) * t102, t14 * t45 - t142 * t15 - t23 * t37 - t133 - t83, t51 * t108 - t11 * t14 + t12 * t15 + t124 * t120 + t3 * t23 + t4 * t24, t77 * t136 + t89 * t28, (-t26 * t77 - t28 * t75) * t44 + (-t10 * t77 - t135 + (t26 * t75 - t28 * t77) * qJD(5)) * t48, t144 * t89 - t28 * t43 + t48 * t31 + t88 * t9, -t48 * t128 - t10 * t88 + t26 * t43 + (-t77 * t114 - t75 * t44) * t144, -t144 * t43 + t130, -t1 * t43 + t23 * t10 + t14 * t26 + t5 * t88 + (t18 * t144 + t134 + (-t144 * t24 + t7 * t48 - t8 * t88) * qJD(5)) * t77 + t139 * t75, t14 * t28 + t2 * t43 + t23 * t9 + (-(-qJD(5) * t24 + t18) * t144 - t134 - (-qJD(5) * t8 + t6) * t88 - t7 * t114) * t75 + t139 * t77; 0, 0, 0, 0, -t81, -t121, 0, 0, 0, 0, 0, t122 * t76, t122 * t78, t142 * t43 + t48 * t37 - t44 * t45 - t130, t83 - t112, 0, 0, 0, 0, 0, -t88 * t128 - t48 * t10 - t44 * t26 + (t43 * t75 + t92 * t77) * t144, -t88 * t31 - t44 * t28 - t136 + (t43 * t77 - t92 * t75) * t144; 0, 0, 0, 0, 0, 0, t78 * t81 * t76, -t123 * t81, 0, 0, 0, -t78 * t121, t76 * t121, (t12 - t16) * t45 - (t11 - t17) * t142 + (t118 * t37 - t36 * t74) * pkin(3), t11 * t16 - t12 * t17 + (-t78 * t112 - t118 * t3 + t4 * t74) * pkin(3), t28 * t96 + t135, (t9 - t147) * t77 + (-t144 * t28 - t10) * t75, -t132 - t145, t90 + t131, -t144 * t45, -t1 * t45 + t64 * t10 - t140 * t77 - t16 * t26 + t84 * t75, t140 * t75 - t16 * t28 + t2 * t45 + t64 * t9 + t84 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142 ^ 2 - t45 ^ 2, t11 * t45 + t12 * t142 + t124, 0, 0, 0, 0, 0, t90 - t131, -t132 + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t26, -t26 ^ 2 + t28 ^ 2, t9 + t147, t141 * t28 + t127, t36, t141 * t2 - t7 * t28 - t75 * t4 + t5, t141 * t1 + t7 * t26 - t77 * t4 - t75 * t6;];
tauc_reg = t19;
