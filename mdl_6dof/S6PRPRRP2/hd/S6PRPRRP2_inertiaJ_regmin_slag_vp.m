% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:38:15
% EndTime: 2019-05-04 23:38:16
% DurationCPUTime: 0.44s
% Computational Cost: add. (382->91), mult. (870->160), div. (0->0), fcn. (998->10), ass. (0->63)
t37 = sin(qJ(5));
t40 = cos(qJ(5));
t46 = -t40 * pkin(5) - t37 * qJ(6);
t20 = -pkin(4) + t46;
t70 = -0.2e1 * t20;
t38 = sin(qJ(4));
t69 = 0.2e1 * t38;
t68 = pkin(4) * t37;
t67 = pkin(4) * t40;
t66 = t37 * pkin(9);
t65 = t40 * pkin(9);
t33 = sin(pkin(11));
t34 = sin(pkin(6));
t35 = cos(pkin(11));
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t15 = (t33 * t42 + t35 * t39) * t34;
t36 = cos(pkin(6));
t41 = cos(qJ(4));
t9 = t15 * t38 - t36 * t41;
t64 = t9 * t37;
t63 = t9 * t40;
t23 = t33 * pkin(2) + pkin(8);
t62 = t23 * t37;
t61 = t23 * t40;
t29 = t37 ^ 2;
t60 = t29 * t38;
t59 = t34 * t39;
t58 = t34 * t42;
t57 = t37 * t38;
t56 = t37 * t40;
t55 = t37 * t41;
t27 = t40 * t38;
t28 = t40 * t41;
t54 = t41 * t23;
t24 = -t35 * pkin(2) - pkin(3);
t18 = -t41 * pkin(4) - t38 * pkin(9) + t24;
t8 = t37 * t18 + t40 * t54;
t31 = t40 ^ 2;
t53 = t29 + t31;
t52 = t41 * qJ(6);
t51 = t41 * t69;
t10 = t15 * t41 + t36 * t38;
t13 = t33 * t59 - t35 * t58;
t2 = t10 * t37 - t13 * t40;
t50 = t2 * t41 + t9 * t57;
t49 = t53 * pkin(9);
t3 = t10 * t40 + t13 * t37;
t48 = t2 * t37 + t3 * t40;
t5 = -t52 + t8;
t17 = t40 * t18;
t6 = -t17 + (pkin(5) + t62) * t41;
t47 = t6 * t37 + t5 * t40;
t45 = -pkin(5) * t37 + t40 * qJ(6);
t32 = t41 ^ 2;
t30 = t38 ^ 2;
t26 = t31 * t38;
t25 = t31 * t30;
t22 = pkin(9) * t55;
t11 = (t23 - t45) * t38;
t7 = -t37 * t54 + t17;
t1 = t9 * t27 + t3 * t41;
t4 = [1, 0, 0, 0, t13 ^ 2 + t15 ^ 2 + t36 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 ^ 2 + t3 ^ 2 + t9 ^ 2; 0, 0, t58, -t59 (-t13 * t35 + t15 * t33) * pkin(2), 0, 0, 0, 0, 0, -t13 * t41, t13 * t38, 0, 0, 0, 0, 0, t50, t1, t50 (t2 * t40 - t3 * t37) * t38, -t1, t9 * t11 + t2 * t6 + t3 * t5; 0, 1, 0, 0 (t33 ^ 2 + t35 ^ 2) * pkin(2) ^ 2, t30, t51, 0, 0, 0, -0.2e1 * t24 * t41, t24 * t69, t25, -0.2e1 * t30 * t56, -0.2e1 * t38 * t28, t37 * t51, t32, 0.2e1 * t30 * t62 - 0.2e1 * t7 * t41, 0.2e1 * t30 * t61 + 0.2e1 * t8 * t41, 0.2e1 * t11 * t57 + 0.2e1 * t6 * t41 (-t37 * t5 + t40 * t6) * t69, -0.2e1 * t11 * t27 - 0.2e1 * t5 * t41, t11 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t38 - t9 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t41 + t47 * t38; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t30 + t25 + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, 0, 0, 0, 0, -t63, t64, -t63, t48, -t64, t48 * pkin(9) + t9 * t20; 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * t23, -t54, t37 * t27, t26 - t60, -t55, -t28, 0, t22 + (-t61 - t68) * t38, pkin(9) * t28 + (t62 - t67) * t38, -t11 * t40 + t20 * t57 + t22, t47, -t11 * t37 + (-pkin(9) * t41 - t20 * t38) * t40, t47 * pkin(9) + t11 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t38, 0, 0, 0, 0, 0, t28, -t55, t28, t26 + t60, t55, -t41 * t20 + t38 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t29, 0.2e1 * t56, 0, 0, 0, 0.2e1 * t67, -0.2e1 * t68, t40 * t70, 0.2e1 * t49, t37 * t70, t53 * pkin(9) ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t3, -t2, 0, t3, -t2 * pkin(5) + t3 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t57, -t41, t7, -t8, t17 + (-0.2e1 * pkin(5) - t62) * t41, t46 * t38, -0.2e1 * t52 + t8, -t6 * pkin(5) + t5 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t27, -t57, 0, t27, t45 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t40, 0, -t66, -t65, -t66, t45, t65, t45 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t27, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
