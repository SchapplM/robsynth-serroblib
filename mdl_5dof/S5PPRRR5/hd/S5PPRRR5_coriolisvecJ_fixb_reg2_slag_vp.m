% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:45
% EndTime: 2019-12-31 17:35:47
% DurationCPUTime: 0.48s
% Computational Cost: add. (640->88), mult. (1402->134), div. (0->0), fcn. (952->6), ass. (0->73)
t41 = cos(qJ(3));
t66 = t41 * qJD(2);
t26 = qJD(3) * pkin(3) + t66;
t38 = sin(qJ(3));
t37 = sin(qJ(4));
t40 = cos(qJ(4));
t21 = t37 * t41 + t40 * t38;
t50 = t21 * qJD(3);
t69 = qJD(4) * t40;
t44 = (t38 * t69 + t50) * qJD(2);
t70 = qJD(4) * t37;
t6 = t26 * t70 + t44;
t33 = qJD(3) + qJD(4);
t39 = cos(qJ(5));
t59 = qJD(3) * t66;
t71 = qJD(2) * t38;
t61 = t37 * t71;
t74 = t33 * t61;
t5 = (qJD(4) * t26 + t59) * t40 - t74;
t16 = t37 * t26 + t40 * t71;
t14 = t33 * pkin(7) + t16;
t36 = sin(qJ(5));
t7 = -t39 * qJD(1) - t36 * t14;
t2 = t7 * qJD(5) + t39 * t5;
t67 = t36 * qJD(1);
t27 = qJD(5) * t67;
t68 = qJD(5) * t39;
t3 = -t14 * t68 - t36 * t5 + t27;
t76 = t39 * t14;
t8 = -t67 + t76;
t89 = t2 * t39 - t3 * t36 + (-t36 * t8 - t39 * t7) * qJD(5);
t20 = t37 * t38 - t40 * t41;
t88 = t20 * t33;
t19 = t20 * qJD(2);
t87 = pkin(3) * t69 + t19;
t18 = t21 * qJD(2);
t28 = t37 * pkin(3) + pkin(7);
t42 = qJD(5) ^ 2;
t85 = (pkin(3) * t70 - t18) * t33 + t28 * t42;
t84 = t33 * pkin(4);
t83 = t6 * t20;
t82 = t88 * t33;
t15 = t40 * t26 - t61;
t13 = -t15 - t84;
t81 = t13 * t68 + t6 * t36;
t10 = t21 * qJD(4) + t50;
t80 = t10 * t33;
t79 = t15 * t33;
t78 = t16 * t33;
t75 = t42 * t36;
t34 = t36 ^ 2;
t35 = t39 ^ 2;
t73 = t34 - t35;
t72 = t34 + t35;
t32 = t33 ^ 2;
t65 = t36 * t32 * t39;
t60 = -t13 * t33 - t5;
t58 = t36 * t33 * t68;
t56 = t36 * t7 - t39 * t8;
t55 = pkin(7) * t42 - t78;
t54 = t21 * t42 + t80;
t53 = qJD(5) * (t15 - t84);
t52 = 0.2e1 * qJD(5) * t88;
t51 = (-pkin(3) * t33 - t26) * qJD(4);
t29 = -t40 * pkin(3) - pkin(4);
t48 = qJD(5) * (t29 * t33 - t87);
t43 = qJD(3) ^ 2;
t31 = t42 * t39;
t23 = -0.2e1 * t58;
t22 = 0.2e1 * t58;
t17 = -0.2e1 * t73 * t33 * qJD(5);
t11 = t13 * qJD(5) * t36;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t31, 0, t56 * qJD(5) - t2 * t36 - t3 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t38, -t43 * t41, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t82, 0, -t15 * t10 - t16 * t88 + t5 * t21 + t83, 0, 0, 0, 0, 0, 0, t36 * t52 - t54 * t39, t54 * t36 + t39 * t52, -t72 * t82, t13 * t10 + t21 * t89 + t56 * t88 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t33 + t37 * t51 - t44, -t19 * t33 + (t51 - t59) * t40 + t74, 0, t15 * t18 + t16 * t19 + (t37 * t5 - t40 * t6 + (-t15 * t37 + t16 * t40) * qJD(4)) * pkin(3), t22, t17, t31, t23, -t75, 0, t11 + t36 * t48 + (-t6 - t85) * t39, t85 * t36 + t39 * t48 + t81, t87 * t33 * t72 + t89, -t13 * t18 + t6 * t29 - t56 * t19 + (t13 * t37 - t56 * t40) * qJD(4) * pkin(3) + t89 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 - t6, -t5 + t79, 0, 0, t22, t17, t31, t23, -t75, 0, t11 + t36 * t53 + (-t55 - t6) * t39, t55 * t36 + t39 * t53 + t81, -t72 * t79 + t89, -t6 * pkin(4) + pkin(7) * t89 - t13 * t16 + t56 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t73 * t32, 0, t65, 0, 0, t27 + t60 * t36 + (t8 - t76) * qJD(5), t60 * t39, 0, 0;];
tauc_reg = t1;
