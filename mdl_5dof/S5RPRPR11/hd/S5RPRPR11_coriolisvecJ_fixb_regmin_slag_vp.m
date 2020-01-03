% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:55
% EndTime: 2019-12-31 18:27:58
% DurationCPUTime: 0.91s
% Computational Cost: add. (1073->182), mult. (2904->235), div. (0->0), fcn. (2162->6), ass. (0->111)
t150 = cos(qJ(3));
t94 = cos(pkin(8));
t126 = t150 * t94;
t116 = qJD(1) * t126;
t93 = sin(pkin(8));
t96 = sin(qJ(3));
t143 = t96 * t93;
t128 = qJD(1) * t143;
t62 = -t116 + t128;
t73 = t150 * t93 + t94 * t96;
t64 = t73 * qJD(1);
t95 = sin(qJ(5));
t97 = cos(qJ(5));
t155 = t62 * t95 + t64 * t97;
t90 = qJD(3) - qJD(5);
t148 = t155 * t90;
t82 = qJD(3) * t116;
t52 = qJD(3) * t128 - t82;
t67 = t73 * qJD(3);
t53 = qJD(1) * t67;
t2 = qJD(5) * t155 - t95 * t52 - t53 * t97;
t163 = t2 + t148;
t127 = t94 * pkin(2) + pkin(1);
t76 = -qJD(1) * t127 + qJD(2);
t162 = -t64 * qJ(4) + t76;
t140 = pkin(6) + qJ(2);
t78 = t140 * t94;
t75 = qJD(1) * t78;
t59 = t96 * t75;
t77 = t140 * t93;
t74 = qJD(1) * t77;
t40 = -t150 * t74 - t59;
t130 = qJD(4) - t40;
t110 = -t97 * t62 + t95 * t64;
t119 = qJD(5) * t110 + t97 * t52 - t95 * t53;
t147 = t110 * t90;
t161 = t119 + t147;
t160 = -t110 ^ 2 + t155 ^ 2;
t129 = qJD(1) * qJD(2);
t121 = qJD(2) * t150;
t113 = qJD(1) * t121;
t120 = qJD(3) * t150;
t139 = t113 * t94 - t120 * t74;
t91 = qJD(3) * qJD(4);
t14 = t91 + (-qJD(3) * t75 - t129 * t93) * t96 + t139;
t4 = t53 * pkin(7) + t14;
t124 = t96 * t129;
t134 = qJD(3) * t96;
t15 = t113 * t93 + t120 * t75 + t124 * t94 - t134 * t74;
t7 = pkin(7) * t52 + t15;
t98 = -pkin(3) - pkin(4);
t9 = t62 * t98 - t162;
t159 = t9 * t155 + t95 * t4 - t97 * t7;
t43 = t150 * t78 - t77 * t96;
t157 = t155 * t110;
t156 = qJD(5) + t90;
t131 = -t64 * pkin(7) + t130;
t154 = t110 * t9 - t97 * t4 - t95 * t7;
t153 = (t40 + t59) * qJD(3) + t93 * t124 - t139;
t152 = t64 ^ 2;
t151 = t53 * pkin(3);
t22 = t62 * pkin(3) + t162;
t149 = t22 * t64;
t146 = t64 * t62;
t41 = t150 * t75 - t74 * t96;
t138 = t93 ^ 2 + t94 ^ 2;
t137 = t52 * qJ(4);
t136 = t62 * qJ(4);
t23 = -t77 * t120 + t94 * t121 + (-qJD(2) * t93 - qJD(3) * t78) * t96;
t133 = t23 * qJD(3);
t24 = qJD(2) * t73 + qJD(3) * t43;
t132 = t24 * qJD(3);
t125 = t138 * qJD(1) ^ 2;
t118 = t90 ^ 2;
t21 = pkin(7) * t62 + t41;
t42 = t150 * t77 + t78 * t96;
t13 = qJD(3) * t98 + t131;
t92 = qJD(3) * qJ(4);
t16 = t21 + t92;
t112 = t97 * t13 - t95 * t16;
t111 = -t95 * t13 - t97 * t16;
t72 = -t126 + t143;
t108 = t72 * t97 - t73 * t95;
t39 = t72 * t95 + t73 * t97;
t107 = t64 * qJD(4) - t137;
t66 = -t120 * t94 + t134 * t93;
t106 = -t66 * qJ(4) + t73 * qJD(4);
t105 = 0.2e1 * t138 * t129;
t104 = t73 * qJ(4) + t127;
t103 = qJD(3) * t41 - t15;
t100 = 0.2e1 * t64 * qJD(3);
t57 = t62 ^ 2;
t37 = pkin(3) * t72 - t104;
t36 = pkin(3) * t64 + t136;
t35 = t92 + t41;
t34 = t82 + (t62 - t128) * qJD(3);
t33 = -t82 + (t62 + t128) * qJD(3);
t32 = -qJD(3) * pkin(3) + t130;
t26 = pkin(7) * t72 + t43;
t25 = -pkin(7) * t73 + t42;
t19 = t72 * t98 + t104;
t18 = pkin(3) * t67 - t106;
t17 = t64 * t98 - t136;
t12 = -t107 + t151;
t11 = t66 * pkin(7) + t24;
t10 = t67 * pkin(7) + t23;
t8 = t67 * t98 + t106;
t6 = qJD(5) * t39 - t95 * t66 - t67 * t97;
t5 = qJD(5) * t108 - t97 * t66 + t95 * t67;
t3 = t53 * t98 + t107;
t1 = [0, 0, 0, 0, 0, t105, qJ(2) * t105, -t52 * t73 - t64 * t66, t52 * t72 - t53 * t73 + t62 * t66 - t64 * t67, -t66 * qJD(3), -t67 * qJD(3), 0, -t127 * t53 + t67 * t76 - t132, t127 * t52 - t66 * t76 - t133, t12 * t72 + t18 * t62 + t22 * t67 + t37 * t53 - t132, -t14 * t72 + t15 * t73 - t23 * t62 + t24 * t64 - t32 * t66 - t35 * t67 - t42 * t52 - t43 * t53, -t12 * t73 - t18 * t64 + t22 * t66 + t37 * t52 + t133, t12 * t37 + t14 * t43 + t15 * t42 + t18 * t22 + t23 * t35 + t24 * t32, -t119 * t39 + t155 * t5, -t108 * t119 - t110 * t5 - t155 * t6 - t2 * t39, -t5 * t90, t6 * t90, 0, t8 * t110 + t19 * t2 - t3 * t108 + t9 * t6 - (-t95 * t10 + t97 * t11 + (-t25 * t95 - t26 * t97) * qJD(5)) * t90, t8 * t155 - t19 * t119 + t3 * t39 + t9 * t5 + (t97 * t10 + t95 * t11 + (t25 * t97 - t26 * t95) * qJD(5)) * t90; 0, 0, 0, 0, 0, -t125, -qJ(2) * t125, 0, 0, 0, 0, 0, t100, -t33, t100, -t57 - t152, t33, t151 + t137 + t35 * t62 + (-qJD(4) - t32) * t64, 0, 0, 0, 0, 0, -t2 + t148, t119 - t147; 0, 0, 0, 0, 0, 0, 0, t146, -t57 + t152, t34, 0, 0, -t64 * t76 + t103, t76 * t62 + t153, -t36 * t62 + t103 - t149, pkin(3) * t52 - t53 * qJ(4) + (t35 - t41) * t64 + (t32 - t130) * t62, -t22 * t62 + t36 * t64 - t153 + 0.2e1 * t91, -t15 * pkin(3) + t14 * qJ(4) + t130 * t35 - t22 * t36 - t32 * t41, -t157, -t160, t161, t163, 0, -t17 * t110 + (t131 * t95 + t97 * t21) * t90 + (-(-qJ(4) * t97 - t95 * t98) * t90 - t111) * qJD(5) + t159, -t17 * t155 + (t131 * t97 - t95 * t21) * t90 + ((-qJ(4) * t95 + t97 * t98) * t90 + t112) * qJD(5) - t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t34, -qJD(3) ^ 2 - t152, -qJD(3) * t35 + t149 + t15, 0, 0, 0, 0, 0, -t110 * t64 - t118 * t95, -t118 * t97 - t155 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t160, -t161, -t163, 0, t111 * t156 - t159, -t112 * t156 + t154;];
tauc_reg = t1;
