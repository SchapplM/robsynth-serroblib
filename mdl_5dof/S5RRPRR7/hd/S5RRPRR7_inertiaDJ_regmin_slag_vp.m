% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:40
% EndTime: 2019-12-31 20:15:41
% DurationCPUTime: 0.32s
% Computational Cost: add. (304->77), mult. (676->106), div. (0->0), fcn. (524->6), ass. (0->64)
t78 = qJD(4) + qJD(5);
t56 = 2 * qJD(3);
t54 = cos(qJ(2));
t59 = -t54 * pkin(1) - pkin(2);
t39 = -pkin(7) + t59;
t77 = -pkin(8) + t39;
t55 = -pkin(2) - pkin(7);
t76 = -pkin(8) + t55;
t49 = sin(qJ(5));
t50 = sin(qJ(4));
t66 = t50 * qJD(4);
t67 = qJD(5) * t49;
t52 = cos(qJ(5));
t53 = cos(qJ(4));
t72 = t52 * t53;
t12 = -t49 * t66 - t50 * t67 + t78 * t72;
t68 = pkin(1) * qJD(2);
t63 = t54 * t68;
t35 = qJD(3) + t63;
t65 = t53 * qJD(4);
t44 = pkin(4) * t65;
t23 = t35 + t44;
t25 = t49 * t53 + t52 * t50;
t51 = sin(qJ(2));
t40 = t51 * pkin(1) + qJ(3);
t48 = t50 * pkin(4);
t32 = t40 + t48;
t75 = t32 * t12 + t23 * t25;
t11 = t78 * t25;
t26 = -t49 * t50 + t72;
t74 = -t32 * t11 + t23 * t26;
t36 = qJD(3) + t44;
t41 = qJ(3) + t48;
t73 = -t41 * t11 + t36 * t26;
t71 = t41 * t12 + t36 * t25;
t70 = t35 * t50 + t40 * t65;
t64 = qJ(3) * qJD(4);
t69 = qJD(3) * t50 + t53 * t64;
t45 = t51 * t68;
t62 = pkin(4) * t67;
t61 = qJD(5) * t52 * pkin(4);
t60 = t55 * t66;
t20 = t77 * t53;
t31 = t76 * t53;
t58 = t50 * t45;
t57 = -t39 * t66 + t53 * t45;
t47 = qJD(3) * t53;
t43 = pkin(8) * t66;
t34 = -0.2e1 * t50 * t65;
t30 = t76 * t50;
t28 = t35 * t53;
t24 = 0.2e1 * (t50 ^ 2 - t53 ^ 2) * qJD(4);
t22 = qJD(4) * t31;
t21 = t43 - t60;
t19 = t77 * t50;
t16 = qJD(4) * t20 + t58;
t15 = t43 + t57;
t6 = -0.2e1 * t26 * t11;
t5 = t52 * t21 - t49 * t22 + (-t30 * t52 - t31 * t49) * qJD(5);
t4 = -t49 * t21 - t52 * t22 + (t30 * t49 - t31 * t52) * qJD(5);
t3 = 0.2e1 * t11 * t25 - 0.2e1 * t26 * t12;
t2 = t52 * t15 - t49 * t16 + (-t19 * t52 - t20 * t49) * qJD(5);
t1 = -t49 * t15 - t52 * t16 + (t19 * t49 - t20 * t52) * qJD(5);
t7 = [0, 0, 0, 0, -0.2e1 * t45, -0.2e1 * t63, 0.2e1 * t45, 0.2e1 * t35, 0.2e1 * t40 * t35 + 0.2e1 * t59 * t45, t34, t24, 0, 0, 0, 0.2e1 * t70, -0.2e1 * t40 * t66 + 0.2e1 * t28, t6, t3, 0, 0, 0, 0.2e1 * t75, 0.2e1 * t74; 0, 0, 0, 0, -t45, -t63, t45, t56 + t63, -pkin(2) * t45 + t35 * qJ(3) + t40 * qJD(3), t34, t24, 0, 0, 0, t69 + t70, t28 + t47 + (-qJ(3) - t40) * t66, t6, t3, 0, 0, 0, t71 + t75, t73 + t74; 0, 0, 0, 0, 0, 0, 0, t56, qJ(3) * t56, t34, t24, 0, 0, 0, 0.2e1 * t69, -0.2e1 * t50 * t64 + 0.2e1 * t47, t6, t3, 0, 0, 0, 0.2e1 * t71, 0.2e1 * t73; 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, t57, -t39 * t65 - t58, 0, 0, -t11, -t12, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, -t60, -t55 * t65, 0, 0, -t11, -t12, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t62, -0.2e1 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
