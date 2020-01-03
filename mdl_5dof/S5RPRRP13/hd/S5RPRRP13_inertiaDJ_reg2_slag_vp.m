% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP13_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:37
% EndTime: 2019-12-31 18:59:41
% DurationCPUTime: 1.18s
% Computational Cost: add. (738->140), mult. (1646->251), div. (0->0), fcn. (1137->4), ass. (0->100)
t51 = sin(qJ(4));
t47 = t51 ^ 2;
t53 = cos(qJ(4));
t49 = t53 ^ 2;
t103 = t47 + t49;
t54 = cos(qJ(3));
t44 = t54 * qJD(3);
t22 = t103 * t44;
t111 = t54 * pkin(7);
t52 = sin(qJ(3));
t72 = t52 * pkin(3) - t111;
t31 = qJ(2) + t72;
t55 = -pkin(1) - pkin(6);
t109 = t52 * t55;
t35 = t53 * t109;
t117 = t51 * t31 + t35;
t104 = t47 - t49;
t48 = t52 ^ 2;
t50 = t54 ^ 2;
t76 = (t48 - t50) * qJD(3);
t33 = t104 * qJD(4);
t67 = pkin(4) * t51 - qJ(5) * t53;
t19 = t67 * qJD(4) - t51 * qJD(5);
t107 = t54 * t19;
t62 = -t55 + t67;
t18 = t62 * t54;
t68 = t53 * pkin(4) + t51 * qJ(5);
t32 = -pkin(3) - t68;
t116 = qJD(3) * (t32 * t52 + t111) - qJD(4) * t18 - t107;
t112 = pkin(7) * t52;
t73 = pkin(3) * t54 + t112;
t60 = qJD(3) * t73 + qJD(2);
t59 = -qJD(4) * t117 + t53 * t60;
t115 = t68 * qJD(4) - t53 * qJD(5);
t114 = 0.2e1 * qJD(2);
t113 = 0.2e1 * qJD(5);
t110 = t51 * t55;
t108 = t53 * t31;
t106 = t54 * t55;
t105 = pkin(7) * t22;
t101 = t48 + t50;
t100 = qJD(3) * t18;
t99 = qJD(4) * t51;
t45 = qJD(4) * t53;
t98 = qJD(4) * t54;
t97 = qJD(4) * t55;
t96 = t52 * qJD(3);
t94 = qJ(2) * qJD(3);
t93 = qJ(5) * qJD(3);
t92 = -0.2e1 * pkin(3) * qJD(4);
t91 = t51 * t109;
t79 = t55 * t44;
t90 = t31 * t45 + t51 * t60 + t53 * t79;
t88 = pkin(7) * t99;
t87 = pkin(7) * t45;
t86 = t51 * t98;
t85 = t51 * t97;
t84 = t53 * t98;
t83 = t51 * t45;
t82 = t55 * t96;
t81 = t53 * t96;
t80 = t52 * t44;
t78 = t54 * t93;
t77 = -pkin(4) + t110;
t40 = 0.2e1 * t80;
t75 = t51 * t81;
t74 = t50 * t83;
t10 = t52 * t77 - t108;
t9 = t52 * qJ(5) + t117;
t70 = t10 * t53 - t51 * t9;
t69 = t10 * t51 + t53 * t9;
t16 = -t91 + t108;
t66 = t117 * t51 + t16 * t53;
t65 = -t117 * t53 + t16 * t51;
t5 = t115 * t54 - t62 * t96;
t61 = -t5 + (t32 * t54 - t112) * qJD(4);
t1 = t78 + (qJD(5) - t85) * t52 + t90;
t2 = t44 * t77 - t59;
t57 = t70 * qJD(4) + t1 * t53 + t2 * t51;
t3 = t52 * t85 - t90;
t4 = -t51 * t79 + t59;
t56 = -t66 * qJD(4) - t3 * t53 - t4 * t51;
t46 = qJ(2) * t114;
t39 = -0.2e1 * t83;
t38 = 0.2e1 * t83;
t28 = -t51 * t96 + t84;
t27 = t44 * t51 + t45 * t52;
t26 = t101 * t45;
t25 = -t81 - t86;
t24 = -t44 * t53 + t52 * t99;
t23 = t101 * t99;
t15 = -0.2e1 * t49 * t80 - 0.2e1 * t74;
t14 = -0.2e1 * t47 * t80 + 0.2e1 * t74;
t13 = t104 * t98 + t75;
t12 = -t51 * t76 + t52 * t84;
t11 = -t104 * t96 + 0.4e1 * t54 * t83;
t8 = -0.2e1 * t52 * t86 - 0.2e1 * t53 * t76;
t7 = t50 * t33 + 0.2e1 * t54 * t75;
t6 = (-0.1e1 + t103) * t40;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t46, -0.2e1 * t80, 0.2e1 * t76, 0, t40, 0, 0, 0.2e1 * qJD(2) * t52 + 0.2e1 * t54 * t94, 0.2e1 * qJD(2) * t54 - 0.2e1 * t52 * t94, 0, t46, t15, 0.2e1 * t7, t8, t14, -0.2e1 * t12, t40, -0.2e1 * t50 * t53 * t97 + 0.2e1 * t4 * t52 + 0.2e1 * (t16 + 0.2e1 * t91) * t44, 0.2e1 * t50 * t85 + 0.2e1 * t3 * t52 + 0.2e1 * (-t117 + 0.2e1 * t35) * t44, 0.2e1 * t66 * t96 + 0.2e1 * (t65 * qJD(4) + t3 * t51 - t4 * t53) * t54, -0.2e1 * t55 ^ 2 * t80 - 0.2e1 * t117 * t3 + 0.2e1 * t16 * t4, t15, t8, -0.2e1 * t7, t40, 0.2e1 * t12, t14, 0.2e1 * (-t51 * t100 - t2) * t52 + 0.2e1 * (-qJD(3) * t10 + t18 * t45 + t5 * t51) * t54, -0.2e1 * t70 * t96 + 0.2e1 * (-t69 * qJD(4) - t1 * t51 + t2 * t53) * t54, 0.2e1 * (t53 * t100 + t1) * t52 + 0.2e1 * (qJD(3) * t9 + t18 * t99 - t5 * t53) * t54, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t18 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t23, 0, -t65 * t44 + (t56 - 0.2e1 * t79) * t52, 0, 0, 0, 0, 0, 0, -t26, 0, -t23, (qJD(3) * t69 - t5) * t54 + (t57 + t100) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, 0, -t44, 0, -t82, -t79, 0, 0, -t13, -t11, t27, t13, -t24, 0, (-t51 * t106 - t53 * t73) * qJD(4) + (t51 * t72 - t35) * qJD(3), (-t53 * t106 + t51 * t73) * qJD(4) + (-t53 * t111 + (pkin(3) * t53 + t110) * t52) * qJD(3), t56, -pkin(3) * t82 + pkin(7) * t56, -t13, t27, t11, 0, t24, t13, -t116 * t51 + t61 * t53, t57, t116 * t53 + t61 * t51, pkin(7) * t57 + t18 * t19 + t5 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t44, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t28, t22, -pkin(3) * t96 + t105, 0, 0, 0, 0, 0, 0, t25, t22, t28, t32 * t96 + t105 - t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -0.2e1 * t33, 0, t39, 0, 0, t51 * t92, t53 * t92, 0, 0, t38, 0, 0.2e1 * t33, 0, 0, t39, -0.2e1 * t19 * t53 + 0.2e1 * t32 * t99, 0, -0.2e1 * t19 * t51 - 0.2e1 * t32 * t45, 0.2e1 * t32 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t28, t44, t4, t3, 0, 0, 0, t25, 0, t44, t28, 0, (0.2e1 * pkin(4) - t110) * t44 + t59, (pkin(4) * t96 - qJ(5) * t98) * t53 + (t52 * t93 + (pkin(4) * qJD(4) - qJD(5)) * t54) * t51, 0.2e1 * t78 + (t113 - t85) * t52 + t90, -t2 * pkin(4) + t1 * qJ(5) + t9 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t24, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, -t24, -t115 * t52 - t44 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, -t99, 0, -t87, t88, 0, 0, 0, t45, 0, 0, t99, 0, -t87, -t115, -t88, -t115 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, qJ(5) * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t25, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
