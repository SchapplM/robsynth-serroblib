% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:08
% EndTime: 2019-12-31 19:13:12
% DurationCPUTime: 1.23s
% Computational Cost: add. (1628->204), mult. (3536->295), div. (0->0), fcn. (2410->6), ass. (0->125)
t91 = sin(qJ(4));
t92 = sin(qJ(3));
t94 = cos(qJ(4));
t95 = cos(qJ(3));
t57 = t91 * t95 + t94 * t92;
t52 = t57 * qJD(1);
t162 = qJD(5) + t52;
t135 = qJD(1) * t95;
t122 = t94 * t135;
t136 = qJD(1) * t92;
t123 = t91 * t136;
t51 = -t122 + t123;
t85 = qJD(3) + qJD(4);
t90 = sin(qJ(5));
t93 = cos(qJ(5));
t40 = -t90 * t51 - t93 * t85;
t165 = t162 * t40;
t104 = t93 * t51 - t90 * t85;
t164 = t104 * t162;
t163 = qJD(5) - t162;
t161 = t52 * t85;
t96 = -pkin(1) - pkin(6);
t66 = t96 * qJD(1) + qJD(2);
t49 = -pkin(7) * t135 + t95 * t66;
t45 = qJD(3) * pkin(3) + t49;
t117 = pkin(7) * qJD(1) - t66;
t133 = qJD(3) * t95;
t47 = t117 * t133;
t160 = (qJD(4) * t45 - t47) * t94;
t14 = -t104 * qJD(5) - t161 * t90;
t157 = pkin(7) - t96;
t60 = t157 * t92;
t61 = t157 * t95;
t38 = -t91 * t60 + t94 * t61;
t134 = qJD(3) * t92;
t54 = t157 * t134;
t55 = qJD(3) * t61;
t15 = -t38 * qJD(4) + t91 * t54 - t94 * t55;
t48 = -pkin(7) * t136 + t92 * t66;
t146 = t91 * t48;
t24 = t94 * t45 - t146;
t21 = -t85 * pkin(4) - t24;
t62 = pkin(3) * t136 + qJD(1) * qJ(2);
t26 = t52 * pkin(4) + t51 * pkin(8) + t62;
t143 = t94 * t95;
t107 = t85 * t143;
t140 = t85 * t123;
t32 = qJD(1) * t107 - t140;
t56 = t91 * t92 - t143;
t77 = t92 * pkin(3) + qJ(2);
t34 = t57 * pkin(4) + t56 * pkin(8) + t77;
t131 = qJD(4) * t94;
t132 = qJD(4) * t91;
t36 = -t92 * t131 - t95 * t132 - t91 * t133 - t94 * t134;
t39 = -t94 * t60 - t91 * t61;
t46 = t117 * t134;
t118 = -t48 * t132 + t91 * t46;
t8 = t118 + t160;
t119 = -t94 * t46 - t91 * t47;
t144 = t94 * t48;
t25 = t91 * t45 + t144;
t9 = t25 * qJD(4) + t119;
t159 = -(qJD(5) * t34 + t15) * t162 - t39 * t32 - (qJD(5) * t26 + t8) * t57 - t9 * t56 + t21 * t36;
t86 = qJD(1) * qJD(2);
t158 = 0.2e1 * t86;
t129 = qJD(5) * t93;
t130 = qJD(5) * t90;
t13 = t85 * t129 + t51 * t130 - t161 * t93;
t156 = t13 * t90;
t155 = t21 * t52;
t154 = t21 * t56;
t153 = t34 * t32;
t37 = -t92 * t132 - t91 * t134 + t107;
t152 = t37 * t85;
t151 = t162 * t51;
t150 = t51 * t52;
t149 = t56 * t13;
t147 = t90 * t32;
t145 = t93 * t32;
t97 = qJD(3) ^ 2;
t142 = t97 * t92;
t141 = t97 * t95;
t126 = qJD(1) * qJD(3);
t120 = t95 * t126;
t59 = pkin(3) * t120 + t86;
t139 = t92 ^ 2 - t95 ^ 2;
t98 = qJD(1) ^ 2;
t138 = -t97 - t98;
t137 = t98 * qJ(2);
t67 = pkin(3) * t133 + qJD(2);
t127 = qJ(2) * qJD(3);
t125 = 0.2e1 * qJD(1);
t124 = pkin(3) * t135;
t121 = -pkin(3) * t85 - t45;
t115 = t93 * t162;
t33 = -t51 * pkin(4) + t52 * pkin(8);
t79 = t91 * pkin(3) + pkin(8);
t112 = qJD(5) * t79 + t124 + t33;
t111 = qJD(5) * t57 + qJD(1);
t22 = t85 * pkin(8) + t25;
t6 = t93 * t22 + t90 * t26;
t110 = t21 * t129 - t6 * t51 + t9 * t90;
t27 = t91 * t49 + t144;
t109 = pkin(3) * t132 - t27;
t28 = t94 * t49 - t146;
t108 = -pkin(3) * t131 + t28;
t106 = -t79 * t32 + t155;
t105 = t90 * t22 - t93 * t26;
t103 = -t105 * t51 + t21 * t130 - t9 * t93;
t102 = t62 * t51 - t119;
t101 = t62 * t52 - t118;
t100 = t56 * t130 + t93 * t36;
t80 = -t94 * pkin(3) - pkin(4);
t35 = t36 * t85;
t23 = t51 ^ 2 - t52 ^ 2;
t20 = t140 + (-t122 - t51) * t85;
t16 = t39 * qJD(4) - t94 * t54 - t91 * t55;
t12 = t37 * pkin(4) - t36 * pkin(8) + t67;
t11 = t32 * pkin(4) + pkin(8) * t161 + t59;
t10 = t93 * t11;
t4 = -t104 * t51 + t115 * t162 + t147;
t3 = -t162 ^ 2 * t90 - t40 * t51 + t145;
t2 = -t104 * t115 + t156;
t1 = (t13 - t165) * t93 + (-t14 + t164) * t90;
t5 = [0, 0, 0, 0, t158, qJ(2) * t158, -0.2e1 * t92 * t120, 0.2e1 * t139 * t126, -t142, -t141, 0, -t96 * t142 + (qJD(2) * t92 + t95 * t127) * t125, -t96 * t141 + (qJD(2) * t95 - t92 * t127) * t125, t161 * t56 - t51 * t36, t161 * t57 + t56 * t32 - t36 * t52 + t51 * t37, t35, -t152, 0, -t16 * t85 + t77 * t32 + t62 * t37 + t67 * t52 + t59 * t57, -t15 * t85 - t161 * t77 + t62 * t36 - t67 * t51 - t59 * t56, -t100 * t104 - t93 * t149, (t104 * t90 - t40 * t93) * t36 + (t156 + t14 * t93 + (-t104 * t93 - t40 * t90) * qJD(5)) * t56, t100 * t162 - t104 * t37 + t13 * t57 - t56 * t145, t56 * t147 - t14 * t57 - t40 * t37 + (t56 * t129 - t90 * t36) * t162, t162 * t37 + t32 * t57, t10 * t57 + t38 * t14 + t16 * t40 - t105 * t37 + (t12 * t162 + t153 + (-t162 * t39 - t22 * t57 - t154) * qJD(5)) * t93 + t159 * t90, t38 * t13 - t16 * t104 - t6 * t37 + (-(-qJD(5) * t39 + t12) * t162 - t153 - (-qJD(5) * t22 + t11) * t57 + qJD(5) * t154) * t90 + t159 * t93; 0, 0, 0, 0, -t98, -t137, 0, 0, 0, 0, 0, t138 * t92, t138 * t95, 0, 0, 0, 0, 0, -qJD(1) * t52 + t35, qJD(1) * t51 - t152, 0, 0, 0, 0, 0, -t57 * t147 + t56 * t14 - t36 * t40 + (-t111 * t93 - t37 * t90) * t162, -t57 * t145 + t149 + t36 * t104 + (t111 * t90 - t37 * t93) * t162; 0, 0, 0, 0, 0, 0, t95 * t98 * t92, -t139 * t98, 0, 0, 0, -t95 * t137, t92 * t137, -t150, t23, 0, t20, 0, -t52 * t124 + t27 * t85 + (t121 * t91 - t144) * qJD(4) + t102, t51 * t124 + t28 * t85 + (t121 * qJD(4) + t47) * t94 + t101, t2, t1, t4, t3, t151, t80 * t14 + t106 * t90 + t109 * t40 + (t108 * t90 - t112 * t93) * t162 + t103, t80 * t13 + t106 * t93 - t109 * t104 + (t108 * t93 + t112 * t90) * t162 + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, t23, 0, t20, 0, t102 + (-qJD(4) + t85) * t25, t24 * t85 + t101 - t160, t2, t1, t4, t3, t151, -pkin(4) * t14 - (-t90 * t24 + t93 * t33) * t162 - t25 * t40 + t90 * t155 + (-t129 * t162 - t147) * pkin(8) + t103, -pkin(4) * t13 + (t93 * t24 + t90 * t33) * t162 + t25 * t104 + t93 * t155 + (t130 * t162 - t145) * pkin(8) + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104 * t40, t104 ^ 2 - t40 ^ 2, t13 + t165, -t14 - t164, t32, t21 * t104 - t163 * t6 - t90 * t8 + t10, t163 * t105 - t90 * t11 + t21 * t40 - t93 * t8;];
tauc_reg = t5;
