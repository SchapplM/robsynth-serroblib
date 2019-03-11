% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:34
% EndTime: 2019-03-09 11:50:40
% DurationCPUTime: 1.86s
% Computational Cost: add. (3569->218), mult. (7940->396), div. (0->0), fcn. (7766->8), ass. (0->125)
t104 = sin(qJ(4));
t106 = cos(qJ(4));
t141 = qJD(4) * t106;
t101 = sin(pkin(10));
t102 = cos(pkin(10));
t105 = sin(qJ(2));
t107 = cos(qJ(2));
t83 = t101 * t107 + t102 * t105;
t132 = t83 * t141;
t143 = qJD(2) * t107;
t144 = qJD(2) * t105;
t78 = -t101 * t144 + t102 * t143;
t165 = t104 * t78 + t132;
t158 = cos(qJ(5));
t124 = t158 * qJD(5);
t164 = qJD(4) * t158 + t124;
t129 = t158 * t106;
t103 = sin(qJ(5));
t146 = t103 * t104;
t85 = -t129 + t146;
t100 = t106 ^ 2;
t147 = -t104 ^ 2 + t100;
t118 = t147 * qJD(4);
t163 = qJD(4) + qJD(5);
t162 = pkin(9) * t83;
t130 = t158 * t104;
t86 = t103 * t106 + t130;
t54 = t163 * t86;
t161 = t54 * pkin(5);
t82 = t101 * t105 - t102 * t107;
t160 = t82 * pkin(4);
t94 = pkin(2) * t101 + pkin(8);
t159 = pkin(9) + t94;
t53 = -t106 * t164 + t146 * t163;
t157 = t86 * t53;
t156 = -qJ(3) - pkin(7);
t131 = -pkin(2) * t107 - pkin(1);
t50 = pkin(3) * t82 - pkin(8) * t83 + t131;
t47 = t106 * t50;
t88 = t156 * t105;
t89 = t156 * t107;
t56 = t101 * t88 - t102 * t89;
t23 = -t104 * t56 - t106 * t162 + t160 + t47;
t149 = t104 * t83;
t52 = t106 * t56;
t152 = t104 * t50 + t52;
t25 = -pkin(9) * t149 + t152;
t24 = t158 * t25;
t155 = t103 * t23 + t24;
t142 = qJD(4) * t104;
t133 = t83 * t142;
t134 = t83 * t146;
t17 = t78 * t130 - t103 * t133 - qJD(5) * t134 + (t103 * t78 + t164 * t83) * t106;
t44 = t86 * t83;
t154 = -t17 * t86 + t44 * t53;
t79 = t159 * t104;
t80 = t159 * t106;
t151 = -t103 * t79 + t158 * t80;
t148 = t106 * t78;
t145 = t104 * t106;
t140 = qJD(5) * t103;
t139 = -0.2e1 * pkin(1) * qJD(2);
t95 = -pkin(2) * t102 - pkin(3);
t138 = 0.2e1 * qJD(4) * t95;
t137 = t158 * pkin(4);
t98 = pkin(2) * t144;
t136 = pkin(4) * t142;
t135 = pkin(4) * t140;
t128 = t104 * t141;
t120 = qJD(2) * t156;
t75 = t107 * qJD(3) + t105 * t120;
t76 = -t105 * qJD(3) + t107 * t120;
t41 = t101 * t76 + t102 * t75;
t77 = t83 * qJD(2);
t42 = pkin(3) * t77 - pkin(8) * t78 + t98;
t13 = -t104 * t42 - t106 * t41 - t141 * t50 + t142 * t56;
t12 = -pkin(9) * t165 - t13;
t121 = -t104 * t41 + t106 * t42;
t9 = -pkin(9) * t148 + t77 * pkin(4) + (-t52 + (-t50 + t162) * t104) * qJD(4) + t121;
t127 = -t103 * t12 + t158 * t9;
t126 = qJD(4) * t159;
t40 = t101 * t75 - t102 * t76;
t55 = -t101 * t89 - t102 * t88;
t123 = -t103 * t25 + t158 * t23;
t122 = -t103 * t80 - t158 * t79;
t119 = -0.4e1 * t83 * t145;
t117 = pkin(4) * t124;
t37 = pkin(4) * t149 + t55;
t16 = t54 * t83 + t78 * t85;
t45 = t129 * t83 - t134;
t116 = -t16 * t85 + t45 * t54;
t115 = t40 * t83 + t55 * t78;
t114 = t53 * t82 - t77 * t86;
t113 = t77 * t83 + t78 * t82;
t112 = -t77 * t94 + t78 * t95;
t111 = t82 * t94 - t83 * t95;
t27 = pkin(4) * t165 + t40;
t87 = -pkin(4) * t106 + t95;
t110 = t104 * t77 + t141 * t82;
t3 = -t103 * t9 - t12 * t158 - t124 * t23 + t140 * t25;
t73 = t104 * t126;
t74 = t106 * t126;
t28 = t103 * t74 + t124 * t79 + t140 * t80 + t158 * t73;
t108 = -t56 * t77 + t115;
t4 = -qJD(5) * t155 + t127;
t29 = -qJD(5) * t151 + t103 * t73 - t158 * t74;
t97 = t137 + pkin(5);
t81 = t83 ^ 2;
t57 = pkin(5) * t85 + t87;
t51 = 0.2e1 * t82 * t77;
t49 = t136 + t161;
t48 = t106 * t77 - t142 * t82;
t33 = -qJ(6) * t85 + t151;
t32 = -qJ(6) * t86 + t122;
t30 = -t54 * t82 - t77 * t85;
t26 = pkin(5) * t44 + t37;
t19 = t53 * qJ(6) - t86 * qJD(6) + t29;
t18 = -qJ(6) * t54 - qJD(6) * t85 - t28;
t14 = -qJD(4) * t152 + t121;
t10 = pkin(5) * t17 + t27;
t6 = -qJ(6) * t44 + t155;
t5 = pkin(5) * t82 - qJ(6) * t45 + t123;
t2 = -qJ(6) * t17 - qJD(6) * t44 - t3;
t1 = t77 * pkin(5) + t16 * qJ(6) - t45 * qJD(6) + t4;
t7 = [0, 0, 0, 0.2e1 * t105 * t143, 0.2e1 * (-t105 ^ 2 + t107 ^ 2) * qJD(2), 0, 0, 0, t105 * t139, t107 * t139, -0.2e1 * t41 * t82 + 0.2e1 * t108, 0.2e1 * t131 * t98 + 0.2e1 * t55 * t40 + 0.2e1 * t56 * t41, 0.2e1 * t100 * t78 * t83 - 0.2e1 * t128 * t81, -0.2e1 * t118 * t81 + t119 * t78, 0.2e1 * t106 * t113 - 0.2e1 * t133 * t82, -0.2e1 * t104 * t113 - 0.2e1 * t132 * t82, t51, 0.2e1 * t104 * t108 + 0.2e1 * t132 * t55 + 0.2e1 * t14 * t82 + 0.2e1 * t47 * t77, 0.2e1 * t106 * t115 + 0.2e1 * t13 * t82 - 0.2e1 * t133 * t55 - 0.2e1 * t152 * t77, -0.2e1 * t45 * t16, 0.2e1 * t16 * t44 - 0.2e1 * t17 * t45, -0.2e1 * t16 * t82 + 0.2e1 * t45 * t77, -0.2e1 * t17 * t82 - 0.2e1 * t44 * t77, t51, 0.2e1 * t123 * t77 + 0.2e1 * t17 * t37 + 0.2e1 * t27 * t44 + 0.2e1 * t4 * t82, -0.2e1 * t155 * t77 - 0.2e1 * t16 * t37 + 0.2e1 * t27 * t45 + 0.2e1 * t3 * t82, -0.2e1 * t1 * t45 + 0.2e1 * t16 * t5 - 0.2e1 * t17 * t6 - 0.2e1 * t2 * t44, 0.2e1 * t1 * t5 + 0.2e1 * t10 * t26 + 0.2e1 * t2 * t6; 0, 0, 0, 0, 0, t143, -t144, 0, -pkin(7) * t143, pkin(7) * t144 (-t101 * t77 - t102 * t78) * pkin(2) (t101 * t41 - t102 * t40) * pkin(2), t118 * t83 + t145 * t78, qJD(4) * t119 + t147 * t78, t110, t48, 0, -t40 * t106 + t112 * t104 + (t104 * t55 - t106 * t111) * qJD(4), t40 * t104 + t112 * t106 + (t104 * t111 + t106 * t55) * qJD(4), -t16 * t86 - t45 * t53, -t116 + t154, -t114, t30, 0, t122 * t77 + t136 * t44 + t17 * t87 + t27 * t85 + t29 * t82 + t37 * t54, t136 * t45 - t151 * t77 - t16 * t87 + t27 * t86 + t28 * t82 - t37 * t53, -t1 * t86 + t16 * t32 - t17 * t33 - t18 * t44 - t19 * t45 - t2 * t85 + t5 * t53 - t54 * t6, t1 * t32 + t10 * t57 + t18 * t6 + t19 * t5 + t2 * t33 + t26 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t128, 0.2e1 * t118, 0, 0, 0, t104 * t138, t106 * t138, -0.2e1 * t157, 0.2e1 * t53 * t85 - 0.2e1 * t54 * t86, 0, 0, 0, 0.2e1 * t136 * t85 + 0.2e1 * t54 * t87, 0.2e1 * t136 * t86 - 0.2e1 * t53 * t87, -0.2e1 * t18 * t85 - 0.2e1 * t19 * t86 + 0.2e1 * t32 * t53 - 0.2e1 * t33 * t54, 0.2e1 * t18 * t33 + 0.2e1 * t19 * t32 + 0.2e1 * t49 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, 0, 0, 0, 0, t48, -t110, 0, 0, 0, 0, 0, t30, t114, t116 + t154, -t1 * t85 + t2 * t86 - t5 * t54 - t53 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t86 - t19 * t85 - t32 * t54 - t33 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t54 * t85 - 0.2e1 * t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133 + t148, -t165, t77, t14, t13, 0, 0, -t16, -t17, t77, t77 * t137 + (-t24 + (-t23 - t160) * t103) * qJD(5) + t127 (-t103 * t77 - t124 * t82) * pkin(4) + t3, t97 * t16 + (-t103 * t17 + (t103 * t45 - t158 * t44) * qJD(5)) * pkin(4), t1 * t97 + (t103 * t2 + (-t103 * t5 + t158 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, -t142, 0, -t94 * t141, t94 * t142, 0, 0, -t53, -t54, 0, t29, t28, t97 * t53 + (-t103 * t54 + (t103 * t86 - t158 * t85) * qJD(5)) * pkin(4), t19 * t97 + (t103 * t18 + (-t103 * t32 + t158 * t33) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, -t141, 0, 0, 0, 0, 0, -t54, t53, 0, -t54 * t97 + (-t103 * t53 + (t103 * t85 + t158 * t86) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t135, -0.2e1 * t117, 0, 0.2e1 * (t137 - t97) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, t77, t4, t3, pkin(5) * t16, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t54, 0, t29, t28, pkin(5) * t53, t19 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t53, 0, -t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t117, 0, -pkin(5) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
