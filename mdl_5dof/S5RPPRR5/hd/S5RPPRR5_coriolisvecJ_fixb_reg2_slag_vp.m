% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:38
% EndTime: 2019-12-31 17:56:40
% DurationCPUTime: 0.39s
% Computational Cost: add. (966->74), mult. (1650->112), div. (0->0), fcn. (788->6), ass. (0->62)
t39 = cos(qJ(4));
t57 = qJD(3) * qJD(1);
t23 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(3);
t17 = t23 * qJD(1) + qJD(3);
t24 = sin(pkin(8)) * pkin(1) + qJ(3);
t21 = qJD(1) * t24;
t37 = sin(qJ(4));
t9 = t39 * t17 - t37 * t21;
t5 = t9 * qJD(4) + t39 * t57;
t56 = qJD(1) - qJD(4);
t68 = t56 * t9;
t79 = t5 + t68;
t40 = qJD(5) ^ 2;
t77 = t56 ^ 2;
t78 = t37 * (t40 + t77);
t70 = t56 * pkin(4);
t7 = -t9 + t70;
t76 = t56 * t7;
t10 = t37 * t17 + t39 * t21;
t6 = qJD(4) * t10 + t37 * t57;
t75 = -t10 * t56 - t6;
t73 = qJD(5) * t56;
t62 = t37 * t23 + t39 * t24;
t12 = t37 * qJD(3) + t62 * qJD(4);
t72 = t12 * t56 + t6;
t71 = t39 * t77;
t36 = sin(qJ(5));
t38 = cos(qJ(5));
t8 = -pkin(7) * t56 + t10;
t3 = -t38 * qJD(2) - t36 * t8;
t1 = t3 * qJD(5) + t38 * t5;
t58 = t36 * qJD(2);
t26 = qJD(5) * t58;
t59 = qJD(5) * t38;
t2 = -t36 * t5 - t8 * t59 + t26;
t69 = t38 * t8;
t4 = -t58 + t69;
t41 = -(t3 * t38 + t36 * t4) * qJD(5) + t1 * t38 - t2 * t36;
t47 = t39 * t23 - t37 * t24;
t11 = t39 * qJD(3) + t47 * qJD(4);
t66 = t11 * t56;
t32 = t36 ^ 2;
t33 = t38 ^ 2;
t61 = t32 - t33;
t60 = t32 + t33;
t55 = t36 * t77 * t38;
t54 = -t5 + t76;
t51 = t36 * t56 * t59;
t48 = t3 * t36 - t38 * t4;
t46 = pkin(7) * t40 - t75;
t45 = qJD(5) * (t7 + t9 + t70);
t14 = -pkin(7) + t62;
t44 = t14 * t40 - t72;
t13 = pkin(4) - t47;
t43 = qJD(5) * (-t13 * t56 - t11 - t7);
t42 = 0.2e1 * t39 * t73;
t29 = t40 * t38;
t28 = t40 * t36;
t19 = 0.2e1 * t51;
t18 = -0.2e1 * t51;
t16 = t61 * t73;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t57, 0.2e1 * t21 * qJD(3), 0, 0, 0, 0, 0, 0, t72, t5 + t66, 0, t10 * t11 - t9 * t12 - t6 * t47 + t5 * t62, t19, -0.2e1 * t16, -t29, t18, t28, 0, t36 * t43 - t44 * t38, t44 * t36 + t38 * t43, -t60 * t66 - t41, -t48 * t11 + t7 * t12 + t6 * t13 + t41 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t29, 0, t48 * qJD(5) - t1 * t36 - t2 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1) ^ 2, -t21 * qJD(1), 0, 0, 0, 0, 0, 0, -t37 * t77, -t71, 0, t79 * t37 + t75 * t39, 0, 0, 0, 0, 0, 0, t36 * t42 - t38 * t78, t36 * t78 + t38 * t42, t60 * t71, (t48 * t56 - t6) * t39 + (t41 - t76) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t79, 0, 0, t18, 0.2e1 * t16, t29, t19, -t28, 0, t36 * t45 - t46 * t38, t46 * t36 + t38 * t45, t60 * t68 + t41, -t6 * pkin(4) + t41 * pkin(7) - t7 * t10 + t48 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t61 * t77, 0, t55, 0, 0, t26 + t54 * t36 + (t4 - t69) * qJD(5), t54 * t38, 0, 0;];
tauc_reg = t15;
