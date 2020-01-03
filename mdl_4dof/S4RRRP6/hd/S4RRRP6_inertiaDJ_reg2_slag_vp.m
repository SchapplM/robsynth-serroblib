% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:10
% EndTime: 2019-12-31 17:19:13
% DurationCPUTime: 0.76s
% Computational Cost: add. (480->122), mult. (1389->245), div. (0->0), fcn. (948->4), ass. (0->94)
t49 = sin(qJ(2));
t105 = -0.4e1 * t49;
t48 = sin(qJ(3));
t44 = t48 ^ 2;
t50 = cos(qJ(3));
t46 = t50 ^ 2;
t89 = t44 - t46;
t51 = cos(qJ(2));
t104 = t49 * qJD(4) + (pkin(5) * qJD(3) + qJ(4) * qJD(2)) * t51;
t103 = 0.2e1 * qJD(3);
t102 = pkin(5) * t48;
t101 = t50 * pkin(3);
t43 = qJD(3) * t50;
t81 = t51 * qJD(2);
t21 = t49 * t43 + t48 * t81;
t75 = pkin(5) * t81;
t14 = t21 * pkin(3) + t75;
t100 = t14 * t48;
t99 = t14 * t50;
t29 = (pkin(3) * t48 + pkin(5)) * t49;
t98 = t29 * t48;
t92 = -qJ(4) - pkin(6);
t31 = t92 * t48;
t97 = t31 * t49;
t32 = t92 * t50;
t96 = t32 * t49;
t95 = t48 * t51;
t94 = t49 * t50;
t93 = t50 * t51;
t61 = -t51 * pkin(2) - t49 * pkin(6);
t57 = -pkin(1) + t61;
t54 = qJD(3) * t57;
t60 = pkin(2) * t49 - pkin(6) * t51;
t55 = t60 * qJD(2);
t91 = -t48 * t55 - t50 * t54;
t41 = t49 * qJD(2);
t70 = t48 * t41;
t90 = pkin(5) * t70 + t50 * t55;
t39 = pkin(5) * t93;
t16 = t48 * t57 + t39;
t45 = t49 ^ 2;
t88 = -t51 ^ 2 + t45;
t87 = qJ(4) * t49;
t86 = qJD(2) * t50;
t85 = qJD(3) * t48;
t84 = qJD(3) * t49;
t83 = qJD(3) * t51;
t80 = pkin(5) * t95;
t79 = -0.2e1 * pkin(1) * qJD(2);
t78 = -0.2e1 * t85;
t77 = pkin(3) * t41;
t76 = pkin(3) * t85;
t74 = t48 * t84;
t73 = t48 * t83;
t72 = t50 * t83;
t71 = t29 * t85;
t69 = t48 * t43;
t68 = t49 * t81;
t67 = t50 * t81;
t40 = -pkin(2) - t101;
t66 = -t40 + t101;
t65 = t88 * qJD(2);
t64 = 0.2e1 * t68;
t63 = t48 * t67;
t62 = t45 * t69;
t59 = pkin(3) * t44 + t40 * t50;
t27 = t50 * t57;
t15 = t27 - t80;
t58 = -t15 * t50 - t16 * t48;
t19 = t67 - t74;
t20 = t50 * t41 + t73;
t3 = t20 * pkin(5) + t91;
t4 = -t16 * qJD(3) + t90;
t53 = t58 * qJD(3) - t3 * t50 - t4 * t48;
t52 = qJ(4) * t74 - t104 * t50 - t48 * t54 + t90;
t37 = -0.2e1 * t68;
t36 = -0.2e1 * t69;
t35 = 0.2e1 * t69;
t28 = -0.2e1 * t89 * qJD(3);
t22 = t70 - t72;
t18 = -t48 * qJD(4) + t92 * t43;
t17 = -t50 * qJD(4) - t92 * t85;
t13 = 0.2e1 * t46 * t68 - 0.2e1 * t62;
t12 = 0.2e1 * t44 * t68 + 0.2e1 * t62;
t11 = t89 * t84 - t63;
t10 = t69 * t105 - t89 * t81;
t9 = -0.2e1 * t48 * t65 + 0.2e1 * t49 * t72;
t8 = 0.2e1 * t49 * t73 + 0.2e1 * t88 * t86;
t7 = -t48 * t87 + t16;
t6 = t89 * t45 * t103 + t63 * t105;
t5 = -t50 * t87 + t27 + (-pkin(3) - t102) * t51;
t2 = (pkin(5) * qJD(2) + qJ(4) * qJD(3)) * t94 + t104 * t48 + t91;
t1 = t52 + t77;
t23 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -0.2e1 * t65, 0, t37, 0, 0, t49 * t79, t51 * t79, 0, 0, t13, t6, t8, t12, t9, t37, 0.2e1 * t15 * t41 - 0.2e1 * t4 * t51 + 0.2e1 * (t45 * t43 + t48 * t64) * pkin(5), -0.2e1 * t16 * t41 - 0.2e1 * t3 * t51 + 0.2e1 * (-t45 * t85 + t50 * t64) * pkin(5), 0.2e1 * t58 * t81 + 0.2e1 * (t3 * t48 - t4 * t50 + (t15 * t48 - t16 * t50) * qJD(3)) * t49, 0.2e1 * pkin(5) ^ 2 * t68 + 0.2e1 * t15 * t4 - 0.2e1 * t16 * t3, t13, t6, t8, t12, t9, t37, 0.2e1 * (qJD(2) * t98 - t1) * t51 + 0.2e1 * (qJD(2) * t5 + t29 * t43 + t100) * t49, 0.2e1 * (t29 * t86 - t2) * t51 + 0.2e1 * (-qJD(2) * t7 - t71 + t99) * t49, 0.2e1 * (-t48 * t7 - t5 * t50) * t81 + 0.2e1 * (-t1 * t50 + t2 * t48 + (t48 * t5 - t50 * t7) * qJD(3)) * t49, 0.2e1 * t5 * t1 + 0.2e1 * t29 * t14 - 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, -t41, 0, -t75, pkin(5) * t41, 0, 0, -t11, t10, t22, t11, t20, 0, (pkin(6) * t93 + (-pkin(2) * t50 + t102) * t49) * qJD(3) + (t61 * t48 - t39) * qJD(2), (pkin(5) * t94 + t60 * t48) * qJD(3) + (t61 * t50 + t80) * qJD(2), t53, -pkin(2) * t75 + t53 * pkin(6), -t11, t10, t22, t11, t20, 0, -t99 - t18 * t51 + (t40 * t95 + t97) * qJD(2) + (t59 * t49 + t98) * qJD(3), t100 - t17 * t51 + (t40 * t93 + t96) * qJD(2) + (t66 * t49 * t48 + t29 * t50) * qJD(3), (-t31 * t81 - t18 * t49 - t2 + (-t5 + t96) * qJD(3)) * t50 + (t32 * t81 + t17 * t49 - t1 + (-t7 + t97) * qJD(3)) * t48, pkin(3) * t71 + t1 * t31 + t14 * t40 - t7 * t17 + t5 * t18 + t2 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t28, 0, t36, 0, 0, pkin(2) * t78, -0.2e1 * pkin(2) * t43, 0, 0, t35, t28, 0, t36, 0, 0, t66 * t78, t59 * t103, -0.2e1 * t17 * t50 - 0.2e1 * t18 * t48 + 0.2e1 * (-t31 * t50 + t32 * t48) * qJD(3), 0.2e1 * t32 * t17 + 0.2e1 * t31 * t18 + 0.2e1 * t40 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t21, t41, t4, t3, 0, 0, 0, 0, t19, 0, -t21, t41, t52 + 0.2e1 * t77, t2, -t19 * pkin(3), t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, -t85, 0, -pkin(6) * t43, pkin(6) * t85, 0, 0, 0, 0, t43, 0, -t85, 0, t18, t17, -pkin(3) * t43, t18 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t19, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t43, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t23;
