% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:42
% EndTime: 2021-01-15 16:04:49
% DurationCPUTime: 1.06s
% Computational Cost: add. (1414->219), mult. (3885->336), div. (0->0), fcn. (2978->10), ass. (0->132)
t81 = cos(pkin(5));
t137 = qJD(1) * t81;
t83 = sin(qJ(3));
t122 = t83 * t137;
t80 = sin(pkin(5));
t138 = qJD(1) * t80;
t84 = sin(qJ(2));
t124 = t84 * t138;
t63 = qJD(2) * pkin(7) + t124;
t86 = cos(qJ(3));
t101 = -t86 * t63 - t122;
t87 = cos(qJ(2));
t123 = t87 * t138;
t106 = qJD(4) + t123;
t131 = qJ(4) * qJD(3);
t121 = t86 * t131;
t166 = t101 * qJD(3) + (-t106 * t83 - t121) * qJD(2);
t139 = cos(pkin(10));
t118 = t139 * t83;
t79 = sin(pkin(10));
t61 = t79 * t86 + t118;
t56 = t61 * qJD(2);
t136 = qJD(2) * t83;
t117 = t139 * t86;
t68 = qJD(2) * t117;
t53 = t79 * t136 - t68;
t50 = qJD(5) + t53;
t165 = -qJD(5) + t50;
t140 = qJD(3) * pkin(3);
t164 = t83 * t140 - t124;
t127 = pkin(3) * t136;
t69 = t86 * t137;
t17 = (-t83 * t63 + t69) * qJD(3) + (t106 * t86 - t83 * t131) * qJD(2);
t3 = -t139 * t166 + t79 * t17;
t72 = t79 * pkin(3) + pkin(8);
t163 = (t56 * pkin(4) + t53 * pkin(8) + qJD(5) * t72 + t127) * t50 + t3;
t145 = qJ(4) + pkin(7);
t153 = t79 * t83;
t65 = t145 * t86;
t39 = t139 * t65 - t145 * t153;
t55 = t61 * qJD(3);
t47 = qJD(2) * t55;
t107 = t3 * t61 - t39 * t47;
t115 = qJD(3) * t145;
t51 = t86 * qJD(4) - t83 * t115;
t95 = -t83 * qJD(4) - t86 * t115;
t98 = t117 - t153;
t143 = -t98 * t123 + t139 * t51 + t79 * t95;
t75 = -t86 * pkin(3) - pkin(2);
t49 = t75 * qJD(2) + qJD(4) - t123;
t16 = t53 * pkin(4) - t56 * pkin(8) + t49;
t28 = -pkin(4) * t98 - t61 * pkin(8) + t75;
t119 = t139 * t17;
t4 = t166 * t79 + t119;
t112 = qJ(4) * qJD(2) + t63;
t35 = t112 * t86 + t122;
t154 = t79 * t35;
t34 = -t112 * t83 + t69;
t32 = t34 + t140;
t7 = t139 * t32 - t154;
t5 = -qJD(3) * pkin(4) - t7;
t58 = t98 * qJD(3);
t162 = (qJD(5) * t16 + t4) * t98 + t5 * t58 + (-qJD(5) * t28 - t143) * t50 + t107;
t82 = sin(qJ(5));
t85 = cos(qJ(5));
t42 = t82 * qJD(3) + t85 * t56;
t130 = qJD(2) * qJD(3);
t120 = t83 * t130;
t66 = t79 * t120;
t48 = qJD(3) * t68 - t66;
t19 = t42 * qJD(5) + t82 * t48;
t161 = pkin(3) * t83;
t132 = t85 * qJD(3);
t133 = qJD(5) * t82;
t18 = qJD(5) * t132 - t56 * t133 + t85 * t48;
t160 = t18 * t82;
t159 = t28 * t47;
t40 = t82 * t56 - t132;
t158 = t40 * t50;
t157 = t42 * t50;
t156 = t42 * t56;
t155 = t56 * t40;
t152 = t80 * t84;
t151 = t80 * t87;
t89 = qJD(2) ^ 2;
t150 = t80 * t89;
t149 = t82 * t47;
t43 = t85 * t47;
t88 = qJD(3) ^ 2;
t147 = t88 * t83;
t146 = t88 * t86;
t144 = -t61 * t123 - t139 * t95 + t79 * t51;
t30 = t139 * t35;
t8 = t79 * t32 + t30;
t52 = pkin(3) * t120 + qJD(2) * t124;
t142 = t83 ^ 2 - t86 ^ 2;
t141 = qJD(2) * pkin(2);
t135 = qJD(2) * t84;
t134 = qJD(5) * t61;
t129 = t84 * t150;
t126 = t80 * t135;
t125 = qJD(2) * t151;
t114 = t50 * t85;
t111 = t83 * t125;
t110 = t86 * t125;
t109 = t55 * pkin(4) - t58 * pkin(8) + t164;
t6 = qJD(3) * pkin(8) + t8;
t1 = t85 * t16 - t82 * t6;
t2 = t82 * t16 + t85 * t6;
t105 = t43 + (-t53 * t82 - t133) * t50;
t102 = -t83 * t152 + t81 * t86;
t59 = t86 * t152 + t81 * t83;
t25 = t79 * t102 + t139 * t59;
t104 = -t85 * t151 - t82 * t25;
t103 = t82 * t151 - t85 * t25;
t100 = -t61 * t133 + t85 * t58;
t99 = t141 * qJD(2);
t97 = t59 * qJD(3);
t94 = -0.2e1 * qJD(3) * t141;
t12 = t139 * t34 - t154;
t93 = -t72 * t47 + (t12 + t5) * t50;
t91 = -t97 - t111;
t73 = -t139 * pkin(3) - pkin(4);
t38 = t145 * t118 + t79 * t65;
t33 = t102 * qJD(3) + t110;
t24 = -t139 * t102 + t79 * t59;
t14 = t47 * pkin(4) - t48 * pkin(8) + t52;
t13 = t85 * t14;
t11 = t139 * t33 + t79 * t91;
t10 = t79 * t34 + t30;
t9 = -t139 * t91 + t79 * t33;
t15 = [0, 0, -t129, -t87 * t150, 0, 0, 0, 0, 0, -t86 * t129 + (-t97 - 0.2e1 * t111) * qJD(3), t83 * t129 + (-t33 - t110) * qJD(3), -t9 * qJD(3) + (t53 * t135 - t47 * t87) * t80, -t11 * qJD(3) + (t56 * t135 - t48 * t87) * t80, -t11 * t53 + t24 * t48 - t25 * t47 + t9 * t56, t8 * t11 + t3 * t24 + t4 * t25 - t7 * t9 + (t49 * t135 - t52 * t87) * t80, 0, 0, 0, 0, 0, (t103 * qJD(5) - t82 * t11 + t85 * t126) * t50 + t104 * t47 + t9 * t40 + t24 * t19, -(t104 * qJD(5) + t85 * t11 + t82 * t126) * t50 + t103 * t47 + t9 * t42 + t24 * t18; 0, 0, 0, 0, 0.2e1 * t86 * t120, -0.2e1 * t142 * t130, t146, -t147, 0, -pkin(7) * t146 + t83 * t94, pkin(7) * t147 + t86 * t94, -t53 * t124 + t75 * t47 + t49 * t55 - t52 * t98 + (t53 * t161 - t144) * qJD(3), -t56 * t124 + t75 * t48 + t49 * t58 + t52 * t61 + (t56 * t161 - t143) * qJD(3), -t143 * t53 + t144 * t56 + t38 * t48 + t4 * t98 - t8 * t55 - t7 * t58 + t107, t143 * t8 - t144 * t7 + t164 * t49 + t3 * t38 + t4 * t39 + t52 * t75, t18 * t85 * t61 + t100 * t42, (-t40 * t85 - t42 * t82) * t58 + (-t160 - t19 * t85 + (t40 * t82 - t42 * t85) * qJD(5)) * t61, t100 * t50 - t18 * t98 + t42 * t55 + t61 * t43, -t61 * t149 + t19 * t98 - t40 * t55 + (-t85 * t134 - t82 * t58) * t50, -t47 * t98 + t50 * t55, t1 * t55 - t13 * t98 + t38 * t19 + t144 * t40 + (t159 + t109 * t50 + (-t39 * t50 + t5 * t61 + t6 * t98) * qJD(5)) * t85 + t162 * t82, t38 * t18 - t2 * t55 + t144 * t42 + (-t159 + (-qJD(5) * t6 + t14) * t98 - t5 * t134 + (qJD(5) * t39 - t109) * t50) * t82 + t162 * t85; 0, 0, 0, 0, -t83 * t89 * t86, t142 * t89, 0, 0, 0, t83 * t99, t86 * t99, t10 * qJD(3) - t53 * t127 - t49 * t56 - t3, -t119 + t49 * t53 + (-t79 * t101 + t12) * qJD(3) + (t79 * t121 + (-pkin(3) * t56 + t79 * t106) * t83) * qJD(2), (-t10 + t8) * t56 + (t12 - t7) * t53 + (-t139 * t48 - t47 * t79) * pkin(3), t7 * t10 - t8 * t12 + (-t49 * t136 - t139 * t3 + t4 * t79) * pkin(3), t42 * t114 + t160, (t18 - t158) * t85 + (-t19 - t157) * t82, t50 * t114 + t149 - t156, t105 + t155, -t50 * t56, -t1 * t56 - t10 * t40 - t163 * t85 + t73 * t19 + t93 * t82, -t10 * t42 + t163 * t82 + t73 * t18 + t2 * t56 + t93 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t56 * qJD(3), -t66 + (t68 - t53) * qJD(3), -t53 ^ 2 - t56 ^ 2, t8 * t53 + t7 * t56 + t52, 0, 0, 0, 0, 0, t105 - t155, -t50 ^ 2 * t85 - t149 - t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t40, -t40 ^ 2 + t42 ^ 2, t18 + t158, t157 - t19, t47, t165 * t2 - t82 * t4 - t5 * t42 + t13, t165 * t1 - t82 * t14 - t85 * t4 + t5 * t40;];
tauc_reg = t15;
