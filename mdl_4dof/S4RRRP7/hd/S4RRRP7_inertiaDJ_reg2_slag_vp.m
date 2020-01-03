% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRP7
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:04
% EndTime: 2019-12-31 17:21:07
% DurationCPUTime: 0.74s
% Computational Cost: add. (470->110), mult. (1326->226), div. (0->0), fcn. (900->4), ass. (0->85)
t47 = cos(qJ(2));
t45 = sin(qJ(2));
t94 = t45 * pkin(6);
t62 = -t47 * pkin(2) - t94;
t29 = -pkin(1) + t62;
t46 = cos(qJ(3));
t96 = pkin(5) * t47;
t37 = t46 * t96;
t44 = sin(qJ(3));
t101 = t44 * t29 + t37;
t40 = t44 ^ 2;
t42 = t46 ^ 2;
t90 = t40 - t42;
t28 = t90 * qJD(3);
t38 = t45 * qJD(2);
t85 = qJD(3) * t47;
t73 = t44 * t85;
t20 = t46 * t38 + t73;
t39 = qJD(3) * t46;
t95 = pkin(6) * t47;
t61 = pkin(2) * t45 - t95;
t52 = t61 * qJD(2);
t4 = t20 * pkin(5) - t29 * t39 - t44 * t52;
t58 = pkin(3) * t44 - qJ(4) * t46;
t53 = pkin(5) + t58;
t17 = t53 * t45;
t18 = t58 * qJD(3) - t44 * qJD(4);
t59 = t46 * pkin(3) + t44 * qJ(4);
t27 = -pkin(2) - t59;
t100 = (-t27 * t47 + t94) * qJD(2) - qJD(3) * t17 - t18 * t45;
t71 = t44 * t38;
t5 = pkin(5) * t71 - qJD(3) * t101 + t46 * t52;
t99 = t59 * qJD(3) - t46 * qJD(4);
t98 = 0.2e1 * qJD(4);
t97 = pkin(5) * t44;
t93 = t46 * t29;
t41 = t45 ^ 2;
t89 = -t47 ^ 2 + t41;
t88 = qJD(2) * t46;
t87 = qJD(3) * t44;
t86 = qJD(3) * t45;
t83 = t47 * qJD(2);
t82 = t47 * qJD(4);
t81 = qJ(4) * qJD(2);
t80 = t44 * t96;
t79 = -0.2e1 * pkin(1) * qJD(2);
t78 = -0.2e1 * pkin(2) * qJD(3);
t77 = pkin(3) * t38;
t76 = pkin(6) * t87;
t75 = pkin(6) * t39;
t74 = pkin(5) * t83;
t72 = t46 * t85;
t70 = t44 * t39;
t69 = t45 * t83;
t68 = t46 * t83;
t67 = t45 * t81;
t66 = t89 * qJD(2);
t65 = 0.2e1 * t69;
t64 = t44 * t68;
t63 = t41 * t70;
t10 = -t47 * qJ(4) + t101;
t12 = -t93 + (pkin(3) + t97) * t47;
t57 = -t10 * t44 + t12 * t46;
t15 = -t80 + t93;
t56 = -t101 * t44 - t15 * t46;
t1 = -t4 + t67 - t82;
t2 = -t5 - t77;
t49 = t57 * qJD(3) + t1 * t46 + t2 * t44;
t48 = t56 * qJD(3) - t4 * t46 - t5 * t44;
t35 = -0.2e1 * t69;
t34 = -0.2e1 * t70;
t33 = 0.2e1 * t70;
t32 = pkin(6) * t72;
t22 = t71 - t72;
t21 = t45 * t39 + t44 * t83;
t19 = -t44 * t86 + t68;
t14 = 0.2e1 * t42 * t69 - 0.2e1 * t63;
t13 = 0.2e1 * t40 * t69 + 0.2e1 * t63;
t11 = t90 * t86 - t64;
t9 = -t44 * t66 + t45 * t72;
t8 = 0.4e1 * t45 * t70 + t90 * t83;
t7 = 0.2e1 * t45 * t73 + 0.2e1 * t89 * t88;
t6 = t41 * t28 - 0.2e1 * t45 * t64;
t3 = t99 * t45 + t53 * t83;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -0.2e1 * t66, 0, t35, 0, 0, t45 * t79, t47 * t79, 0, 0, t14, 0.2e1 * t6, t7, t13, 0.2e1 * t9, t35, 0.2e1 * t15 * t38 - 0.2e1 * t5 * t47 + 0.2e1 * (t41 * t39 + t44 * t65) * pkin(5), -0.2e1 * t101 * t38 - 0.2e1 * t4 * t47 + 0.2e1 * (-t41 * t87 + t46 * t65) * pkin(5), 0.2e1 * t56 * t83 + 0.2e1 * (t4 * t44 - t46 * t5 + (-t101 * t46 + t15 * t44) * qJD(3)) * t45, 0.2e1 * pkin(5) ^ 2 * t69 - 0.2e1 * t101 * t4 + 0.2e1 * t15 * t5, t14, t7, -0.2e1 * t6, t35, -0.2e1 * t9, t13, 0.2e1 * (qJD(2) * t17 * t44 + t2) * t47 + 0.2e1 * (-qJD(2) * t12 + t17 * t39 + t3 * t44) * t45, 0.2e1 * t57 * t83 + 0.2e1 * (-t1 * t44 + t2 * t46 + (-t10 * t46 - t12 * t44) * qJD(3)) * t45, 0.2e1 * (-t17 * t88 - t1) * t47 + 0.2e1 * (qJD(2) * t10 + t17 * t87 - t3 * t46) * t45, 0.2e1 * t10 * t1 + 0.2e1 * t12 * t2 + 0.2e1 * t17 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, -t38, 0, -t74, pkin(5) * t38, 0, 0, -t11, -t8, t22, t11, t20, 0, t32 + (-pkin(2) * t46 + t97) * t86 + (t62 * t44 - t37) * qJD(2), (pkin(5) * t45 * t46 + t61 * t44) * qJD(3) + (t62 * t46 + t80) * qJD(2), t48, -pkin(2) * t74 + t48 * pkin(6), -t11, t22, t8, 0, -t20, t11, t32 + (t27 * t86 - t3) * t46 - t100 * t44, t49, (-t3 + (t27 * t45 + t95) * qJD(3)) * t44 + t100 * t46, t49 * pkin(6) + t17 * t18 + t3 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -0.2e1 * t28, 0, t34, 0, 0, t44 * t78, t46 * t78, 0, 0, t33, 0, 0.2e1 * t28, 0, 0, t34, -0.2e1 * t18 * t46 + 0.2e1 * t27 * t87, 0, -0.2e1 * t18 * t44 - 0.2e1 * t27 * t39, 0.2e1 * t27 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t21, t38, t5, t4, 0, 0, 0, t19, 0, t38, t21, 0, t5 + 0.2e1 * t77, (-pkin(3) * t83 - qJ(4) * t86) * t46 + (-t47 * t81 + (pkin(3) * qJD(3) - qJD(4)) * t45) * t44, -t4 + 0.2e1 * t67 - 0.2e1 * t82, -t2 * pkin(3) + t1 * qJ(4) + t10 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t87, 0, -t75, t76, 0, 0, 0, t39, 0, 0, t87, 0, -t75, -t99, -t76, -t99 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, qJ(4) * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t19, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t16;
