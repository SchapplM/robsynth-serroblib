% Calculate minimal parameter regressor of joint inertia matrix for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:33:59
% EndTime: 2019-05-04 20:34:00
% DurationCPUTime: 0.53s
% Computational Cost: add. (481->107), mult. (1278->185), div. (0->0), fcn. (1587->12), ass. (0->69)
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t49 = -t42 * pkin(5) - t39 * qJ(6);
t24 = -pkin(4) + t49;
t75 = -0.2e1 * t24;
t40 = sin(qJ(4));
t74 = -0.2e1 * t40;
t43 = cos(qJ(4));
t73 = 0.2e1 * t43;
t72 = pkin(4) * t39;
t71 = pkin(4) * t42;
t70 = pkin(9) * t39;
t69 = pkin(9) * t42;
t68 = t39 * pkin(10);
t67 = t42 * pkin(10);
t33 = sin(pkin(12));
t35 = sin(pkin(6));
t38 = cos(pkin(6));
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t36 = cos(pkin(12));
t37 = cos(pkin(7));
t60 = t36 * t37;
t34 = sin(pkin(7));
t62 = t34 * t41;
t10 = t38 * t62 + (t33 * t44 + t41 * t60) * t35;
t19 = -t35 * t36 * t34 + t38 * t37;
t5 = t10 * t40 - t19 * t43;
t66 = t5 * t39;
t65 = t5 * t42;
t20 = -t43 * t37 + t40 * t62;
t64 = t20 * t39;
t63 = t20 * t42;
t61 = t34 * t44;
t59 = t39 * t40;
t58 = t39 * t42;
t57 = t39 * t43;
t28 = t42 * t40;
t56 = t42 * t43;
t25 = -t43 * pkin(4) - t40 * pkin(10) - pkin(3);
t17 = pkin(9) * t56 + t39 * t25;
t30 = t39 ^ 2;
t32 = t42 ^ 2;
t55 = t30 + t32;
t54 = t43 * qJ(6);
t53 = t40 * t73;
t6 = t10 * t43 + t19 * t40;
t9 = -t38 * t61 + (t33 * t41 - t44 * t60) * t35;
t2 = t6 * t39 - t9 * t42;
t52 = t2 * t43 + t5 * t59;
t21 = t40 * t37 + t43 * t62;
t11 = t39 * t21 + t42 * t61;
t51 = t11 * t43 + t20 * t59;
t3 = t9 * t39 + t6 * t42;
t50 = t2 * t39 + t3 * t42;
t48 = -pkin(5) * t39 + t42 * qJ(6);
t12 = t42 * t21 - t39 * t61;
t47 = t11 * t39 + t12 * t42;
t13 = -t54 + t17;
t23 = t42 * t25;
t14 = -t23 + (pkin(5) + t70) * t43;
t46 = t13 * t42 + t14 * t39;
t31 = t40 ^ 2;
t26 = pkin(10) * t57;
t18 = (pkin(9) - t48) * t40;
t16 = -pkin(9) * t57 + t23;
t7 = t12 * t43 + t20 * t28;
t1 = t5 * t28 + t3 * t43;
t4 = [1, t38 ^ 2 + (t33 ^ 2 + t36 ^ 2) * t35 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 ^ 2 + t3 ^ 2 + t5 ^ 2; 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * t11 + t3 * t12 + t5 * t20; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 ^ 2 + t12 ^ 2 + t20 ^ 2; 0, 0, 0, -t9, -t10, 0, 0, 0, 0, 0, -t9 * t43, t9 * t40, 0, 0, 0, 0, 0, t52, t1, t52 (t2 * t42 - t3 * t39) * t40, -t1, t3 * t13 + t2 * t14 + t5 * t18; 0, 0, 0, t61, -t62, 0, 0, 0, 0, 0, t43 * t61, -t40 * t61, 0, 0, 0, 0, 0, t51, t7, t51 (t11 * t42 - t12 * t39) * t40, -t7, t11 * t14 + t12 * t13 + t20 * t18; 0, 0, 1, 0, 0, t31, t53, 0, 0, 0, pkin(3) * t73, pkin(3) * t74, t32 * t31, -0.2e1 * t31 * t58, t56 * t74, t39 * t53, t43 ^ 2, -0.2e1 * t16 * t43 + 0.2e1 * t31 * t70, 0.2e1 * t17 * t43 + 0.2e1 * t31 * t69, 0.2e1 * t14 * t43 + 0.2e1 * t18 * t59, 0.2e1 * (-t13 * t39 + t14 * t42) * t40, -0.2e1 * t13 * t43 - 0.2e1 * t18 * t28, t13 ^ 2 + t14 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t65, t66, -t65, t50, -t66, t50 * pkin(10) + t5 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, 0, 0, 0, 0, -t63, t64, -t63, t47, -t64, t47 * pkin(10) + t20 * t24; 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t40 * pkin(9), -t43 * pkin(9), t39 * t28 (-t30 + t32) * t40, -t57, -t56, 0, t26 + (-t69 - t72) * t40, pkin(10) * t56 + (t70 - t71) * t40, -t18 * t42 + t24 * t59 + t26, t46, -t18 * t39 + (-pkin(10) * t43 - t24 * t40) * t42, t46 * pkin(10) + t18 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t30, 0.2e1 * t58, 0, 0, 0, 0.2e1 * t71, -0.2e1 * t72, t42 * t75, 0.2e1 * t55 * pkin(10), t39 * t75, t55 * pkin(10) ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t3, -t2, 0, t3, -t2 * pkin(5) + t3 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -t11, 0, t12, -t11 * pkin(5) + t12 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t59, -t43, t16, -t17, t23 + (-0.2e1 * pkin(5) - t70) * t43, t49 * t40, -0.2e1 * t54 + t17, -t14 * pkin(5) + t13 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t42, 0, -t68, -t67, -t68, t48, t67, t48 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t28, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
