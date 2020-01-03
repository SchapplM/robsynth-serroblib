% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:35
% EndTime: 2020-01-03 12:00:37
% DurationCPUTime: 0.60s
% Computational Cost: add. (645->75), mult. (1640->129), div. (0->0), fcn. (1198->8), ass. (0->71)
t78 = cos(pkin(9));
t82 = sin(qJ(2));
t56 = t78 * t82;
t53 = pkin(1) * t56;
t84 = cos(qJ(2));
t71 = t84 * pkin(1);
t60 = t71 + pkin(2);
t77 = sin(pkin(9));
t30 = t77 * t60 + t53;
t43 = sin(qJ(4));
t55 = t77 * t82;
t52 = pkin(1) * t55;
t47 = t78 * t60 - t52;
t46 = pkin(3) + t47;
t83 = cos(qJ(4));
t45 = t83 * t46;
t58 = t84 * t77;
t79 = pkin(1) * qJD(2);
t48 = (-t56 - t58) * t79;
t35 = qJD(2) * t52;
t57 = t78 * t84;
t54 = pkin(1) * t57;
t50 = qJD(2) * t54 - t35;
t74 = -qJD(4) * t45 - t43 * t48 - t83 * t50;
t76 = qJD(4) * t43;
t5 = t30 * t76 + t74;
t42 = sin(qJ(5));
t40 = t42 ^ 2;
t44 = cos(qJ(5));
t41 = t44 ^ 2;
t80 = -t40 - t41;
t87 = t80 * t5;
t59 = t78 * pkin(2) + pkin(3);
t51 = t83 * t59;
t33 = qJD(4) * t51;
t65 = t77 * pkin(2);
t61 = t43 * t65;
t25 = qJD(4) * t61 - t33;
t86 = t80 * t25;
t12 = -t43 * t30 + t45;
t10 = -pkin(4) - t12;
t39 = t44 * qJD(5);
t24 = t83 * t30;
t13 = t43 * t46 + t24;
t49 = -t43 * t50 + t83 * t48;
t6 = t13 * qJD(4) - t49;
t85 = t10 * t39 + t6 * t42;
t34 = t43 * t59;
t38 = t83 * t65;
t31 = t38 + t34;
t26 = t31 * qJD(4);
t29 = -t61 + t51;
t27 = -pkin(4) - t29;
t81 = t26 * t42 + t27 * t39;
t75 = t42 * qJD(5);
t73 = pkin(4) * t75;
t72 = pkin(4) * t39;
t70 = t42 * t39;
t7 = t10 * t75;
t69 = -t6 * t44 + t7;
t28 = pkin(8) + t31;
t66 = t80 * t28;
t17 = t27 * t75;
t64 = -t26 * t44 + t17;
t63 = qJD(2) * t71;
t62 = t82 * t79;
t37 = -0.2e1 * t70;
t36 = 0.2e1 * t70;
t32 = 0.2e1 * (-t40 + t41) * qJD(5);
t11 = pkin(8) + t13;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t62, -0.2e1 * t63, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t48, -0.2e1 * t50, 0, -0.2e1 * t30 * t35 + 0.2e1 * (t30 * t54 + t47 * (-pkin(1) * t58 - t53)) * qJD(2), 0, 0, 0, 0, 0, 0, -0.2e1 * t6, 0.2e1 * t5, 0, -0.2e1 * t12 * t6 - 0.2e1 * t13 * t5, t36, t32, 0, t37, 0, 0, 0.2e1 * t69, 0.2e1 * t85, 0.2e1 * t87, 0.2e1 * t10 * t6 + 0.2e1 * t11 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t63, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t50, 0, (t78 * t48 + t50 * t77) * pkin(2), 0, 0, 0, 0, 0, 0, (-t38 - 0.2e1 * t34 - t24 - t43 * (t57 - t55) * pkin(1)) * qJD(4) + t49, -t33 + (t65 + t30) * t76 + t74, 0, -t12 * t26 - t13 * t25 - t6 * t29 - t5 * t31, t36, t32, 0, t37, 0, 0, t17 + t7 + (-t26 - t6) * t44, t81 + t85, t86 + t87, t10 * t26 + t11 * t86 + t6 * t27 + t5 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t26, 0.2e1 * t25, 0, -0.2e1 * t31 * t25 - 0.2e1 * t29 * t26, t36, t32, 0, t37, 0, 0, 0.2e1 * t64, 0.2e1 * t81, 0.2e1 * t86, 0.2e1 * t25 * t66 + 0.2e1 * t27 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, t36, t32, 0, t37, 0, 0, t69 - t73, -t72 + t85, t87, -t6 * pkin(4) + pkin(8) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t25, 0, 0, t36, t32, 0, t37, 0, 0, t64 - t73, -t72 + t81, t86, -t26 * pkin(4) + pkin(8) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t32, 0, t37, 0, 0, -0.2e1 * t73, -0.2e1 * t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t75, 0, -t11 * t39 + t42 * t5, t11 * t75 + t44 * t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t75, 0, t42 * t25 - t28 * t39, t44 * t25 + t28 * t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t75, 0, -pkin(8) * t39, pkin(8) * t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
