% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:42
% EndTime: 2021-01-15 16:33:48
% DurationCPUTime: 0.87s
% Computational Cost: add. (1253->171), mult. (3183->241), div. (0->0), fcn. (2190->6), ass. (0->120)
t79 = sin(qJ(2));
t121 = t79 * qJD(1);
t64 = qJD(2) * pkin(6) + t121;
t106 = pkin(7) * qJD(2) + t64;
t81 = cos(qJ(3));
t40 = t106 * t81;
t77 = sin(qJ(4));
t32 = t77 * t40;
t78 = sin(qJ(3));
t39 = t106 * t78;
t35 = qJD(3) * pkin(3) - t39;
t80 = cos(qJ(4));
t109 = t80 * t35 - t32;
t56 = t77 * t81 + t80 * t78;
t47 = t56 * qJD(2);
t128 = t47 * qJ(5);
t8 = t109 - t128;
t74 = qJD(3) + qJD(4);
t156 = t74 * t56;
t22 = t156 * qJD(2);
t123 = qJD(4) * t77;
t82 = cos(qJ(2));
t120 = t82 * qJD(1);
t125 = qJD(3) * t78;
t27 = -t64 * t125 + (-pkin(7) * t125 + t81 * t120) * qJD(2);
t124 = qJD(3) * t81;
t28 = -t64 * t124 + (-pkin(7) * t124 - t78 * t120) * qJD(2);
t92 = -(qJD(4) * t35 + t27) * t80 + t40 * t123 - t77 * t28;
t155 = -t22 * qJ(5) - t92;
t118 = qJD(2) * qJD(3);
t154 = -0.2e1 * t118;
t139 = t80 * t81;
t140 = t77 * t78;
t55 = -t139 + t140;
t42 = t55 * t79;
t149 = pkin(6) + pkin(7);
t113 = qJD(3) * t149;
t57 = t78 * t113;
t58 = t81 * t113;
t59 = t149 * t78;
t60 = t149 * t81;
t99 = t77 * t59 - t80 * t60;
t153 = t99 * qJD(4) + t56 * t120 + t77 * t57 - t80 * t58;
t122 = qJD(4) * t80;
t96 = t55 * t82;
t152 = -qJD(1) * t96 + t59 * t122 + t60 * t123 + t80 * t57 + t77 * t58;
t83 = qJD(3) ^ 2;
t84 = qJD(2) ^ 2;
t151 = (t83 + t84) * t79;
t97 = pkin(3) * t125 - t121;
t150 = t47 ^ 2;
t7 = t74 * pkin(4) + t8;
t148 = t7 - t8;
t147 = pkin(3) * t74;
t101 = t74 * t140;
t25 = -t81 * t122 - t80 * t124 + t101;
t146 = -t25 * qJ(5) + t56 * qJD(5) - t153;
t145 = qJ(5) * t156 + t55 * qJD(5) + t152;
t114 = qJD(2) * t139;
t127 = qJD(2) * t78;
t115 = t77 * t127;
t45 = -t114 + t115;
t144 = t45 * t74;
t143 = t47 * t45;
t142 = t47 * t74;
t73 = -t81 * pkin(3) - pkin(2);
t52 = t73 * qJD(2) - t120;
t141 = t52 * t47;
t34 = t80 * t40;
t138 = t83 * t78;
t137 = t83 * t81;
t136 = -t80 * t39 - t32;
t112 = t81 * t118;
t135 = -qJD(4) * t114 - t80 * t112;
t116 = pkin(3) * t127;
t53 = qJD(2) * t121 + qJD(3) * t116;
t134 = t78 ^ 2 - t81 ^ 2;
t132 = qJD(2) * pkin(2);
t21 = qJD(2) * t101 + t135;
t131 = t21 * qJ(5);
t129 = t45 * qJ(5);
t126 = qJD(2) * t79;
t107 = t45 * pkin(4) + qJD(5);
t24 = t107 + t52;
t119 = qJD(5) + t24;
t111 = -t77 * t27 + t80 * t28;
t108 = t77 * t39 - t34;
t105 = pkin(4) * t156 + t97;
t102 = t82 * t154;
t16 = t22 * pkin(4) + t53;
t100 = -t77 * t35 - t34;
t98 = qJD(2) * t132;
t94 = -0.2e1 * qJD(3) * t132;
t93 = -t74 * t115 - t135;
t91 = t100 * qJD(4) + t111;
t89 = t52 * t45 + t92;
t88 = t91 + t131;
t87 = (-t34 + (-t35 - t147) * t77) * qJD(4) + t111;
t85 = t119 * t45 - t155;
t72 = t80 * pkin(3) + pkin(4);
t61 = t122 * t147;
t44 = t45 ^ 2;
t41 = t56 * t79;
t38 = t55 * pkin(4) + t73;
t30 = t47 * pkin(4) + t116;
t19 = -t55 * qJ(5) - t99;
t18 = -t56 * qJ(5) - t80 * t59 - t77 * t60;
t17 = -t44 + t150;
t15 = t142 - t22;
t14 = t93 + t144;
t13 = t74 * t42 - t82 * t47;
t12 = -qJD(2) * t96 - t156 * t79;
t11 = -t128 + t136;
t10 = t108 + t129;
t9 = -t100 - t129;
t4 = t45 * t126 + t13 * t74 - t82 * t22;
t3 = -t12 * t74 + t47 * t126 + t82 * t21;
t2 = -t47 * qJD(5) + t88;
t1 = -t45 * qJD(5) + t155;
t5 = [0, 0, -t84 * t79, -t84 * t82, 0, 0, 0, 0, 0, t78 * t102 - t81 * t151, t81 * t102 + t78 * t151, 0, 0, 0, 0, 0, t4, t3, t4, t3, -t12 * t45 - t13 * t47 - t41 * t21 + t42 * t22, -t1 * t42 + t9 * t12 + t24 * t126 + t7 * t13 - t16 * t82 - t2 * t41; 0, 0, 0, 0, 0.2e1 * t78 * t112, t134 * t154, t137, -t138, 0, -pkin(6) * t137 + t78 * t94, pkin(6) * t138 + t81 * t94, -t21 * t56 - t47 * t25, -t156 * t47 + t21 * t55 - t56 * t22 + t25 * t45, -t25 * t74, -t156 * t74, 0, t153 * t74 + t156 * t52 + t73 * t22 + t97 * t45 + t53 * t55, t152 * t74 - t73 * t21 - t52 * t25 + t97 * t47 + t53 * t56, t105 * t45 - t146 * t74 + t156 * t24 + t16 * t55 + t38 * t22, t105 * t47 + t145 * t74 + t16 * t56 - t38 * t21 - t24 * t25, -t1 * t55 + t145 * t45 + t146 * t47 - t156 * t9 + t18 * t21 - t19 * t22 - t2 * t56 + t7 * t25, t1 * t19 + t105 * t24 - t145 * t9 - t146 * t7 + t16 * t38 + t2 * t18; 0, 0, 0, 0, -t78 * t84 * t81, t134 * t84, 0, 0, 0, t78 * t98, t81 * t98, t143, t17, t14, t15, 0, -t108 * t74 - t45 * t116 - t141 + t87, -t47 * t116 + t136 * t74 - t61 + t89, -t10 * t74 - t119 * t47 - t30 * t45 + t131 + t87, t11 * t74 - t30 * t47 - t61 + t85, t72 * t21 + (t10 + t9) * t47 + (t11 - t7) * t45 + (-t22 * t77 + (-t45 * t80 + t47 * t77) * qJD(4)) * pkin(3), -t7 * t10 - t9 * t11 + t2 * t72 - t24 * t30 + (t1 * t77 + (-t7 * t77 + t80 * t9) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t17, t14, t15, 0, -t100 * t74 - t141 + t91, t109 * t74 + t89, t9 * t74 + (-t107 - t24) * t47 + t88, -t150 * pkin(4) + t8 * t74 + t85, t21 * pkin(4) - t148 * t45, t148 * t9 + (-t24 * t47 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 + t142, t93 - t144, -t44 - t150, t9 * t45 + t7 * t47 + t16;];
tauc_reg = t5;
