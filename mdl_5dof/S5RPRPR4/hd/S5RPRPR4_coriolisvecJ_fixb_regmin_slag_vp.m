% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:23:10
% EndTime: 2022-01-23 09:23:13
% DurationCPUTime: 0.87s
% Computational Cost: add. (1120->162), mult. (2798->240), div. (0->0), fcn. (1960->8), ass. (0->107)
t100 = sin(qJ(5));
t101 = sin(qJ(3));
t103 = cos(qJ(3));
t96 = sin(pkin(9));
t98 = cos(pkin(9));
t81 = t98 * t101 + t96 * t103;
t130 = qJD(1) * t81;
t102 = cos(qJ(5));
t125 = qJD(1) * t103;
t118 = t98 * t125;
t126 = qJD(1) * t101;
t72 = t96 * t126 - t118;
t63 = t102 * t72;
t29 = -t100 * t130 - t63;
t93 = qJD(3) + qJD(5);
t135 = t29 * t93;
t124 = qJD(5) * t100;
t74 = t81 * qJD(3);
t66 = qJD(1) * t74;
t122 = qJD(1) * qJD(3);
t115 = t103 * t122;
t116 = t101 * t122;
t85 = t96 * t116;
t67 = t98 * t115 - t85;
t4 = -qJD(5) * t63 - t100 * t66 + t102 * t67 - t124 * t130;
t148 = t4 - t135;
t109 = t100 * t72 - t102 * t130;
t106 = qJD(5) * t109 - t100 * t67 - t102 * t66;
t136 = t109 * t93;
t147 = t106 - t136;
t146 = t109 * t29;
t89 = sin(pkin(8)) * pkin(1) + pkin(6);
t131 = qJ(4) + t89;
t145 = t109 ^ 2 - t29 ^ 2;
t91 = -cos(pkin(8)) * pkin(1) - pkin(2);
t82 = -t103 * pkin(3) + t91;
t129 = qJD(1) * t82;
t71 = qJD(4) + t129;
t38 = t72 * pkin(4) + t71;
t140 = t72 * pkin(7);
t111 = t131 * qJD(1);
t54 = t101 * qJD(2) + t103 * t111;
t134 = t98 * t54;
t132 = qJD(3) * pkin(3);
t53 = t103 * qJD(2) - t111 * t101;
t48 = t53 + t132;
t19 = t96 * t48 + t134;
t9 = t19 - t140;
t144 = t9 * t124 - t38 * t29;
t121 = qJD(1) * qJD(4);
t36 = t53 * qJD(3) + t103 * t121;
t37 = -t54 * qJD(3) - t101 * t121;
t10 = -t96 * t36 + t98 * t37;
t2 = -t67 * pkin(7) + t10;
t11 = t98 * t36 + t96 * t37;
t3 = -t66 * pkin(7) + t11;
t143 = -t100 * t3 + t102 * t2 + t38 * t109;
t142 = qJD(5) - t93;
t141 = pkin(3) * t96;
t139 = t130 * pkin(7);
t138 = pkin(3) * t101;
t80 = t96 * t101 - t98 * t103;
t39 = t100 * t81 + t102 * t80;
t77 = t80 * qJD(3);
t12 = -qJD(5) * t39 - t100 * t74 - t102 * t77;
t137 = t12 * t93;
t44 = t96 * t54;
t21 = t98 * t53 - t44;
t112 = qJD(3) * t131;
t57 = t103 * qJD(4) - t101 * t112;
t58 = -t101 * qJD(4) - t103 * t112;
t23 = t98 * t57 + t96 * t58;
t78 = t131 * t101;
t79 = t131 * t103;
t35 = -t96 * t78 + t98 * t79;
t133 = t101 ^ 2 - t103 ^ 2;
t84 = qJD(1) * t91;
t104 = qJD(3) ^ 2;
t128 = t104 * t101;
t127 = t104 * t103;
t120 = t101 * t132;
t119 = pkin(3) * t126;
t18 = t98 * t48 - t44;
t20 = -t96 * t53 - t134;
t22 = -t96 * t57 + t98 * t58;
t34 = -t98 * t78 - t96 * t79;
t113 = 0.2e1 * t130;
t7 = qJD(3) * pkin(4) - t139 + t18;
t110 = -t100 * t7 - t102 * t9;
t40 = -t100 * t80 + t102 * t81;
t108 = 0.2e1 * qJD(3) * t84;
t105 = qJD(1) ^ 2;
t90 = t98 * pkin(3) + pkin(4);
t87 = pkin(3) * t116;
t56 = t74 * pkin(4) + t120;
t55 = pkin(4) * t130 + t119;
t52 = t80 * pkin(4) + t82;
t43 = t66 * pkin(4) + t87;
t25 = -t80 * pkin(7) + t35;
t24 = -t81 * pkin(7) + t34;
t17 = -t74 * pkin(7) + t23;
t16 = t77 * pkin(7) + t22;
t15 = t21 - t139;
t14 = t20 + t140;
t13 = qJD(5) * t40 - t100 * t77 + t102 * t74;
t8 = t13 * t93;
t1 = [0, 0, 0, 0, 0.2e1 * t101 * t115, -0.2e1 * t133 * t122, t127, -t128, 0, t101 * t108 - t89 * t127, t103 * t108 + t89 * t128, t82 * t66 + t71 * t74 + (t22 + (qJD(1) * t80 + t72) * t138) * qJD(3), t82 * t67 - t71 * t77 + (t113 * t138 - t23) * qJD(3), -t10 * t81 - t11 * t80 - t130 * t22 + t18 * t77 - t19 * t74 - t23 * t72 - t34 * t67 - t35 * t66, t10 * t34 + t11 * t35 + t18 * t22 + t19 * t23 + (t71 + t129) * t120, -t109 * t12 + t4 * t40, t106 * t40 + t109 * t13 + t12 * t29 - t4 * t39, t137, -t8, 0, -t56 * t29 - t52 * t106 + t43 * t39 + t38 * t13 + (-t100 * t17 + t102 * t16 + (-t100 * t24 - t102 * t25) * qJD(5)) * t93, -t56 * t109 + t52 * t4 + t43 * t40 + t38 * t12 - (t100 * t16 + t102 * t17 + (-t100 * t25 + t102 * t24) * qJD(5)) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t127, -t74 * qJD(3), t77 * qJD(3), t130 * t74 - t81 * t66 + t80 * t67 + t77 * t72, -t10 * t80 + t11 * t81 - t18 * t74 - t19 * t77, 0, 0, 0, 0, 0, -t8, -t137; 0, 0, 0, 0, -t101 * t105 * t103, t133 * t105, 0, 0, 0, -t84 * t126, -t84 * t125, -t20 * qJD(3) - t72 * t119 - t130 * t71 + t10, t21 * qJD(3) - t119 * t130 + t71 * t72 - t11, (t19 + t20) * t130 + (-t18 + t21) * t72 + (-t66 * t96 - t67 * t98) * pkin(3), -t18 * t20 - t19 * t21 + (t10 * t98 + t11 * t96 - t71 * t126) * pkin(3), t146, t145, t148, t147, 0, t55 * t29 - (-t100 * t15 + t102 * t14) * t93 + ((-t100 * t90 - t102 * t141) * t93 + t110) * qJD(5) + t143, -t102 * t3 - t100 * t2 + t55 * t109 + (t100 * t14 + t102 * t15) * t93 + (-(-t100 * t141 + t102 * t90) * t93 - t102 * t7) * qJD(5) + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113 * qJD(3), -t85 + (-t72 + t118) * qJD(3), -t130 ^ 2 - t72 ^ 2, t130 * t18 + t19 * t72 + t87, 0, 0, 0, 0, 0, -t106 - t136, t4 + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t145, t148, t147, 0, t142 * t110 + t143, (-t9 * t93 - t2) * t100 + (-t142 * t7 - t3) * t102 + t144;];
tauc_reg = t1;
