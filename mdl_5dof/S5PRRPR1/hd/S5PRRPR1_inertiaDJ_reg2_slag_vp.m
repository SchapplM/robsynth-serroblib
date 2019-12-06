% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:51
% EndTime: 2019-12-05 16:15:54
% DurationCPUTime: 0.56s
% Computational Cost: add. (423->69), mult. (1016->121), div. (0->0), fcn. (854->6), ass. (0->59)
t51 = sin(pkin(9));
t52 = cos(pkin(9));
t73 = t51 ^ 2 + t52 ^ 2;
t89 = t73 * qJD(4);
t90 = 0.2e1 * t89;
t55 = cos(qJ(3));
t72 = pkin(2) * qJD(3);
t70 = t55 * t72;
t39 = qJD(4) + t70;
t88 = t73 * t39;
t87 = t55 * pkin(2);
t54 = sin(qJ(3));
t43 = t54 * pkin(2) + qJ(4);
t86 = -pkin(7) - t43;
t85 = cos(qJ(5));
t68 = t85 * t51;
t53 = sin(qJ(5));
t80 = t53 * t52;
t34 = t68 + t80;
t29 = t34 * qJD(5);
t67 = t85 * t52;
t33 = t53 * t51 - t67;
t84 = t33 * t29;
t62 = qJD(5) * t85;
t71 = qJD(5) * t53;
t28 = t51 * t71 - t52 * t62;
t83 = t34 * t28;
t44 = -t52 * pkin(4) - pkin(3);
t82 = t44 * t28;
t81 = t44 * t29;
t79 = -pkin(7) - qJ(4);
t37 = t44 - t87;
t45 = t54 * t72;
t77 = t37 * t29 + t33 * t45;
t76 = -t37 * t28 + t34 * t45;
t69 = t53 * t86;
t48 = t52 * pkin(7);
t30 = t52 * t43 + t48;
t57 = t86 * t68;
t15 = -t53 * t30 + t57;
t16 = t85 * t30 + t51 * t69;
t4 = -qJD(5) * t57 - t39 * t67 + (qJD(5) * t30 + t39 * t51) * t53;
t5 = -t30 * t62 - t39 * t80 + (-qJD(5) * t69 - t85 * t39) * t51;
t66 = t15 * t28 - t16 * t29 + t4 * t33 - t5 * t34;
t65 = t79 * t51;
t38 = t52 * qJ(4) + t48;
t56 = t85 * t65;
t61 = t85 * qJD(4);
t11 = -qJD(5) * t56 - t52 * t61 + (qJD(4) * t51 + qJD(5) * t38) * t53;
t12 = -t38 * t62 - qJD(4) * t80 + (-t79 * t71 - t61) * t51;
t21 = -t53 * t38 + t56;
t22 = t85 * t38 + t53 * t65;
t63 = t11 * t33 - t12 * t34 + t21 * t28 - t22 * t29;
t59 = t51 * t45;
t58 = t52 * t45;
t18 = -0.2e1 * t83;
t17 = 0.2e1 * t84;
t6 = 0.2e1 * t33 * t28 - 0.2e1 * t34 * t29;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t83 + 0.2e1 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29 * t15 - t28 * t16 - t33 * t5 - t34 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t45, -0.2e1 * t70, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t58, 0.2e1 * t59, 0.2e1 * t88, 0.2e1 * (-pkin(3) - t87) * t45 + 0.2e1 * t43 * t88, t18, t6, 0, t17, 0, 0, 0.2e1 * t77, 0.2e1 * t76, 0.2e1 * t66, 0.2e1 * t15 * t5 - 0.2e1 * t16 * t4 + 0.2e1 * t37 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t11 - t33 * t12 - t29 * t21 - t28 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t70, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t59, t89 + t88, -pkin(3) * t45 + qJ(4) * t88 + t43 * t89, t18, t6, 0, t17, 0, 0, t77 + t81, t76 - t82, t63 + t66, -t16 * t11 + t15 * t12 + t5 * t21 - t4 * t22 + t44 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, qJ(4) * t90, t18, t6, 0, t17, 0, 0, 0.2e1 * t81, -0.2e1 * t82, 0.2e1 * t63, -0.2e1 * t22 * t11 + 0.2e1 * t21 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, t29, -t28, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, -t29, 0, t5, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, -t29, 0, t12, t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
