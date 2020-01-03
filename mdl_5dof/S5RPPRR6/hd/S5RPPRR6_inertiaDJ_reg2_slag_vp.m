% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:07
% EndTime: 2019-12-31 17:58:10
% DurationCPUTime: 0.72s
% Computational Cost: add. (1024->100), mult. (2168->185), div. (0->0), fcn. (2040->8), ass. (0->74)
t41 = cos(pkin(9));
t34 = sin(pkin(8)) * pkin(1) + qJ(3);
t87 = pkin(6) + t34;
t28 = t87 * t41;
t40 = sin(pkin(9));
t86 = cos(qJ(4));
t61 = t87 * t86;
t65 = t86 * qJD(3);
t43 = sin(qJ(4));
t73 = t43 * qJD(3);
t76 = qJD(4) * t43;
t45 = (-qJD(4) * t61 - t73) * t40 - t28 * t76 + t41 * t65;
t66 = qJD(4) * t86;
t26 = t40 * t76 - t41 * t66;
t31 = t86 * t40 + t43 * t41;
t27 = t31 * qJD(4);
t30 = t43 * t40 - t86 * t41;
t95 = t30 * t26 - t31 * t27;
t32 = -cos(pkin(8)) * pkin(1) - pkin(2) - t41 * pkin(3);
t47 = t30 * pkin(4) - t31 * pkin(7) + t32;
t94 = -qJD(5) * t47 - t45;
t42 = sin(qJ(5));
t38 = t42 ^ 2;
t44 = cos(qJ(5));
t39 = t44 ^ 2;
t64 = qJD(5) * (t38 - t39);
t67 = t43 * t87;
t16 = t86 * t28 - t40 * t67;
t90 = t27 * pkin(4);
t91 = t26 * pkin(7);
t60 = t90 + t91;
t75 = qJD(5) * t42;
t1 = t16 * t75 - t42 * t60 + t94 * t44;
t74 = qJD(5) * t44;
t2 = -t16 * t74 + t94 * t42 + t44 * t60;
t3 = -t42 * t16 + t44 * t47;
t4 = t44 * t16 + t42 * t47;
t55 = t3 * t42 - t4 * t44;
t93 = qJD(5) * t55 + t1 * t42 - t2 * t44;
t15 = t43 * t28 + t40 * t61;
t7 = t28 * t66 + t41 * t73 + (-qJD(4) * t67 + t65) * t40;
t92 = t15 * t7;
t89 = t7 * t42;
t88 = t7 * t44;
t85 = t30 * t27;
t84 = t31 * t26;
t83 = t38 * t26;
t23 = t39 * t26;
t82 = t42 * t27;
t81 = t44 * t26;
t80 = t44 * t27;
t79 = -t30 * t81 + t31 * t80;
t72 = 0.2e1 * t85;
t71 = -0.2e1 * pkin(4) * qJD(5);
t70 = t42 * t81;
t69 = t31 * t75;
t68 = t42 * t74;
t29 = t31 ^ 2;
t63 = t29 * t68;
t62 = 0.2e1 * (t40 ^ 2 + t41 ^ 2) * qJD(3);
t59 = pkin(4) * t26 - pkin(7) * t27;
t58 = pkin(4) * t31 + pkin(7) * t30;
t56 = t3 * t44 + t4 * t42;
t54 = t15 * t27 + t7 * t30;
t52 = -t42 * t26 + t31 * t74;
t12 = t69 + t81;
t13 = t30 * t74 + t82;
t48 = -qJD(5) * t56 - t1 * t44 - t2 * t42;
t18 = t31 * t23;
t17 = t31 * t83;
t11 = t30 * t75 - t80;
t9 = t23 + t83;
t6 = t31 * t64 + t70;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t34 * t62, -0.2e1 * t84, 0.2e1 * t95, 0, t72, 0, 0, 0.2e1 * t32 * t27, -0.2e1 * t32 * t26, -0.2e1 * t15 * t26 - 0.2e1 * t16 * t27 - 0.2e1 * t45 * t30 + 0.2e1 * t7 * t31, 0.2e1 * t45 * t16 + 0.2e1 * t92, -0.2e1 * t18 - 0.2e1 * t63, 0.2e1 * t29 * t64 + 0.4e1 * t31 * t70, -0.2e1 * t30 * t69 + 0.2e1 * t79, -0.2e1 * t17 + 0.2e1 * t63, -0.2e1 * t30 * t52 - 0.2e1 * t31 * t82, t72, 0.2e1 * t15 * t52 + 0.2e1 * t2 * t30 + 0.2e1 * t3 * t27 + 0.2e1 * t31 * t89, 0.2e1 * t1 * t30 - 0.2e1 * t12 * t15 - 0.2e1 * t4 * t27 + 0.2e1 * t31 * t88, 0.2e1 * t56 * t26 + 0.2e1 * t93 * t31, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t26 + t45 * t31 + t54, 0, 0, 0, 0, 0, 0, 0, t44 * t95 + t79, 0, t26 * t55 + t31 * t48 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t84 + 0.2e1 * t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17 - 0.2e1 * t18 + 0.2e1 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t13, t9, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, -t27, 0, -t7, -t45, 0, 0, -t6, -0.4e1 * t31 * t68 - t23 + t83, t13, t6, -t11, 0, -t88 + t59 * t42 + (t15 * t42 - t44 * t58) * qJD(5), t89 + t59 * t44 + (t15 * t44 + t42 * t58) * qJD(5), t48, -t7 * pkin(4) + pkin(7) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t26, 0, 0, 0, 0, 0, 0, 0, 0, t11, t13, -t9, -t90 + (-t38 - t39) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t68, -0.2e1 * t64, 0, -0.2e1 * t68, 0, 0, t42 * t71, t44 * t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, -t52, t27, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t74, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, -t75, 0, -pkin(7) * t74, pkin(7) * t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
