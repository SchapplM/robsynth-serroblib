% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:49
% EndTime: 2019-12-31 18:17:50
% DurationCPUTime: 0.32s
% Computational Cost: add. (739->87), mult. (1442->107), div. (0->0), fcn. (740->6), ass. (0->74)
t43 = cos(qJ(3));
t84 = sin(pkin(8)) * pkin(1);
t65 = qJD(3) * t84;
t59 = qJD(1) * t65;
t32 = cos(pkin(8)) * pkin(1) + pkin(2);
t30 = t32 * qJD(1);
t41 = sin(qJ(3));
t78 = t41 * t30;
t14 = qJD(3) * t78 + t43 * t59;
t66 = qJD(1) * t84;
t16 = t43 * t66 + t78;
t35 = qJD(1) + qJD(3);
t72 = t35 * qJ(4);
t13 = t16 + t72;
t82 = t13 * t35;
t86 = t14 - t82;
t42 = cos(qJ(5));
t12 = t42 * t14;
t40 = sin(qJ(5));
t44 = -pkin(3) - pkin(7);
t15 = -t43 * t30 + t41 * t66;
t50 = qJD(4) + t15;
t8 = t35 * t44 + t50;
t5 = qJD(2) * t42 + t40 * t8;
t71 = t5 * qJD(5);
t2 = t12 - t71;
t4 = -qJD(2) * t40 + t42 * t8;
t69 = qJD(5) * t40;
t1 = t4 * qJD(5) + t40 * t14;
t83 = t1 * t40;
t48 = -(t2 + t71) * t42 + t4 * t69 - t83;
t34 = t35 ^ 2;
t33 = t35 * qJD(4);
t70 = qJD(3) * t43;
t51 = -t30 * t70 + t41 * t59;
t10 = -t33 + t51;
t68 = qJD(5) * t42;
t85 = -t10 * t40 + t13 * t68;
t81 = t16 * t35;
t18 = t32 * t70 - t41 * t65;
t17 = -qJD(4) - t18;
t80 = t17 * t35;
t52 = t32 * t41 + t43 * t84;
t19 = t52 * qJD(3);
t79 = t19 * t35;
t45 = qJD(5) ^ 2;
t77 = t45 * t40;
t76 = t45 * t42;
t75 = -t34 - t45;
t36 = t40 ^ 2;
t37 = t42 ^ 2;
t74 = t36 - t37;
t73 = t36 + t37;
t67 = t42 * t34 * t40;
t22 = qJ(4) + t52;
t62 = t22 * t35 + t19;
t61 = t32 * t43 - t41 * t84;
t60 = t35 * t40 * t68;
t58 = -pkin(3) - t61;
t57 = -t14 + t81;
t56 = t14 + t79;
t55 = t4 * t42 + t40 * t5;
t54 = -t10 * t22 - t13 * t17;
t21 = -pkin(7) + t58;
t53 = -t21 * t45 - t80;
t49 = -t10 * qJ(4) + t13 * t50;
t47 = -t15 * t35 + t51;
t46 = t83 + t2 * t42 + (-t4 * t40 + t42 * t5) * qJD(5);
t27 = -0.2e1 * t60;
t26 = 0.2e1 * t60;
t20 = 0.2e1 * t74 * t35 * qJD(5);
t11 = -pkin(3) * t35 + t50;
t7 = t10 * t42;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t18 * t35 + t51, 0, -t14 * t61 + t15 * t19 + t16 * t18 - t51 * t52, 0, 0, 0, 0, 0, 0, 0, t56, -t10 - t80, t11 * t19 + t14 * t58 + t54, t27, t20, -t77, t26, -t76, 0, t40 * t53 + t62 * t68 + t85, -t7 + t53 * t42 + (-t13 - t62) * t69, -t73 * t79 + t48, t19 * t55 + t21 * t46 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t77, 0, -qJD(5) * t55 + t1 * t42 - t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0.2e1 * t33 - t47, -t14 * pkin(3) - t11 * t16 + t49, t27, t20, -t77, t26, -t76, 0, -t16 * t68 - t44 * t77 + (qJ(4) * t68 + t40 * t50) * t35 + t85, -t7 + (t35 * t50 - t44 * t45) * t42 + (-t13 + t16 - t72) * t69, t73 * t81 + t48, -t16 * t55 + t44 * t46 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t86, 0, 0, 0, 0, 0, 0, t75 * t40, t75 * t42, 0, -t48 - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t74 * t34, 0, -t67, 0, 0, -t42 * t82 + t12, -t86 * t40, 0, 0;];
tauc_reg = t3;
