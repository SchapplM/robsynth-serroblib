% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:15
% EndTime: 2019-12-31 17:26:17
% DurationCPUTime: 0.74s
% Computational Cost: add. (1122->116), mult. (2693->233), div. (0->0), fcn. (2295->6), ass. (0->86)
t48 = sin(qJ(4));
t46 = t48 ^ 2;
t51 = cos(qJ(4));
t47 = t51 ^ 2;
t85 = t46 - t47;
t102 = t85 * qJD(4);
t100 = -pkin(6) - pkin(5);
t50 = sin(qJ(2));
t34 = t100 * t50;
t52 = cos(qJ(2));
t35 = t100 * t52;
t49 = sin(qJ(3));
t97 = cos(qJ(3));
t22 = t49 * t34 - t35 * t97;
t70 = t97 * t52;
t91 = t49 * t50;
t28 = -t70 + t91;
t90 = t49 * t52;
t29 = t50 * t97 + t90;
t44 = -pkin(2) * t52 - pkin(1);
t57 = -pkin(3) * t28 + pkin(7) * t29 - t44;
t54 = t51 * t57;
t7 = -t48 * t22 - t54;
t8 = t51 * t22 - t48 * t57;
t103 = -t48 * t7 + t51 * t8;
t101 = qJD(2) + qJD(3);
t69 = qJD(2) * t100;
t31 = t50 * t69;
t65 = qJD(2) * t70;
t12 = qJD(3) * t22 - t100 * t65 + t49 * t31;
t96 = t12 * t29;
t10 = t12 * t48;
t21 = -t34 * t97 - t49 * t35;
t95 = t21 * t12;
t94 = t21 * t49;
t68 = t97 * qJD(3);
t19 = t101 * t91 - t52 * t68 - t65;
t93 = t29 * t19;
t20 = t101 * t29;
t92 = t48 * t20;
t89 = t51 * t19;
t88 = t51 * t20;
t45 = qJD(4) * t51;
t87 = t21 * t45 + t10;
t43 = -pkin(2) * t97 - pkin(3);
t84 = pkin(2) * qJD(3);
t73 = t49 * t84;
t86 = t43 * t45 + t48 * t73;
t82 = qJD(4) * t48;
t81 = t50 * qJD(2);
t80 = t52 * qJD(2);
t79 = 0.2e1 * t28 * t20;
t78 = -0.2e1 * pkin(1) * qJD(2);
t77 = t48 * t89;
t76 = pkin(3) * t82;
t75 = pkin(3) * t45;
t74 = pkin(2) * t81;
t72 = t48 * t45;
t71 = t50 * t80;
t26 = t29 ^ 2;
t67 = t26 * t72;
t66 = pkin(2) * t68;
t64 = t48 * t8 + t51 * t7;
t42 = pkin(2) * t49 + pkin(7);
t63 = t28 * t42 - t29 * t43;
t62 = t43 * t82 - t51 * t73;
t61 = -t19 * t48 + t29 * t45;
t60 = t29 * t82 + t89;
t59 = t28 * t82 - t88;
t58 = (t46 + t47) * t97;
t56 = (-t28 * t97 + t29 * t49) * qJD(3);
t55 = pkin(3) * t20 + pkin(7) * t19 + t74;
t11 = qJD(3) * t21 - t31 * t97 - t69 * t90;
t2 = qJD(4) * t54 + t51 * t11 + t22 * t82 - t48 * t55;
t3 = -qJD(4) * t8 + t48 * t11 + t51 * t55;
t1 = -qJD(4) * t64 - t2 * t51 - t3 * t48;
t53 = pkin(2) * t56 - t19 * t43 - t20 * t42;
t38 = -0.2e1 * t72;
t37 = 0.2e1 * t72;
t27 = -0.2e1 * t102;
t24 = t58 * t84;
t17 = t21 * t82;
t14 = t28 * t45 + t92;
t6 = t102 * t29 + t77;
t4 = t19 * t85 - 0.4e1 * t29 * t72;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t71, 0.2e1 * (-t50 ^ 2 + t52 ^ 2) * qJD(2), 0, -0.2e1 * t71, 0, 0, t50 * t78, t52 * t78, 0, 0, -0.2e1 * t93, 0.2e1 * t19 * t28 - 0.2e1 * t20 * t29, 0, t79, 0, 0, 0.2e1 * t20 * t44 + 0.2e1 * t28 * t74, -0.2e1 * t19 * t44 + 0.2e1 * t29 * t74, 0.2e1 * t11 * t28 - 0.2e1 * t19 * t21 - 0.2e1 * t22 * t20 + 0.2e1 * t96, -0.2e1 * t11 * t22 + 0.2e1 * t44 * t74 + 0.2e1 * t95, -0.2e1 * t47 * t93 - 0.2e1 * t67, 0.2e1 * t102 * t26 + 0.4e1 * t29 * t77, -0.2e1 * t28 * t60 + 0.2e1 * t29 * t88, -0.2e1 * t46 * t93 + 0.2e1 * t67, -0.2e1 * t28 * t61 - 0.2e1 * t29 * t92, t79, 0.2e1 * t10 * t29 + 0.2e1 * t7 * t20 + 0.2e1 * t21 * t61 + 0.2e1 * t3 * t28, 0.2e1 * t2 * t28 - 0.2e1 * t8 * t20 - 0.2e1 * t21 * t60 + 0.2e1 * t51 * t96, 0.2e1 * t64 * t19 + 0.2e1 * (-qJD(4) * t103 + t2 * t48 - t3 * t51) * t29, -0.2e1 * t2 * t8 + 0.2e1 * t3 * t7 + 0.2e1 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, -t81, 0, -pkin(5) * t80, pkin(5) * t81, 0, 0, 0, 0, -t19, 0, -t20, 0, -t12, t11, (t19 * t97 - t20 * t49 + t56) * pkin(2), (-t97 * t12 - t11 * t49 + (t22 * t97 + t94) * qJD(3)) * pkin(2), -t6, t4, t14, t6, -t59, 0, t17 + (-qJD(4) * t63 - t12) * t51 + t53 * t48, t51 * t53 + t63 * t82 + t87, t1, t12 * t43 + (t103 * t97 + t94) * t84 + t1 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t73, -0.2e1 * t66, 0, 0, t37, t27, 0, t38, 0, 0, 0.2e1 * t62, 0.2e1 * t86, 0.2e1 * t24, 0.2e1 * (t42 * t58 + t43 * t49) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, -t20, 0, -t12, t11, 0, 0, -t6, t4, t14, t6, -t59, 0, t17 + (pkin(3) * t19 - pkin(7) * t20) * t48 + (-t12 + (-pkin(3) * t29 - pkin(7) * t28) * qJD(4)) * t51, pkin(3) * t60 + pkin(7) * t59 + t87, t1, -t12 * pkin(3) + pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t66, 0, 0, t37, t27, 0, t38, 0, 0, t62 - t76, -t75 + t86, t24, (-pkin(3) * t49 + pkin(7) * t58) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t27, 0, t38, 0, 0, -0.2e1 * t76, -0.2e1 * t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, -t61, t20, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, -t82, 0, -t42 * t45 - t48 * t66, t42 * t82 - t51 * t66, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, -t82, 0, -pkin(7) * t45, pkin(7) * t82, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
