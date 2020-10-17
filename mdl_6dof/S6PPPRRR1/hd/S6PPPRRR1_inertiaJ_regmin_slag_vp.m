% Calculate minimal parameter regressor of joint inertia matrix for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPPRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_inertiaJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:33:58
% EndTime: 2019-05-04 19:33:59
% DurationCPUTime: 0.37s
% Computational Cost: add. (313->80), mult. (949->152), div. (0->0), fcn. (1256->16), ass. (0->60)
t41 = sin(qJ(5));
t61 = -0.2e1 * t41;
t44 = cos(qJ(5));
t60 = 0.2e1 * t44;
t43 = cos(qJ(6));
t59 = pkin(5) * t43;
t40 = sin(qJ(6));
t58 = pkin(10) * t40;
t32 = sin(pkin(8));
t42 = sin(qJ(4));
t57 = t32 * t42;
t45 = cos(qJ(4));
t56 = t32 * t45;
t33 = sin(pkin(7));
t39 = cos(pkin(6));
t55 = t33 * t39;
t35 = cos(pkin(14));
t37 = cos(pkin(8));
t54 = t35 * t37;
t36 = cos(pkin(13));
t38 = cos(pkin(7));
t53 = t36 * t38;
t52 = t40 * t41;
t51 = t43 * t40;
t50 = t43 * t41;
t49 = t43 * t44;
t48 = t44 * t40;
t47 = t41 * t60;
t30 = sin(pkin(14));
t31 = sin(pkin(13));
t34 = sin(pkin(6));
t12 = t35 * t55 + (-t30 * t31 + t35 * t53) * t34;
t21 = -t34 * t36 * t33 + t39 * t38;
t46 = t12 * t37 + t21 * t32;
t29 = t43 ^ 2;
t28 = t41 ^ 2;
t27 = t40 ^ 2;
t24 = -t44 * pkin(5) - t41 * pkin(11) - pkin(4);
t23 = t41 * t37 + t44 * t57;
t22 = -t44 * t37 + t41 * t57;
t20 = -t32 * t33 * t35 + t37 * t38;
t19 = pkin(10) * t49 + t40 * t24;
t18 = -pkin(10) * t48 + t43 * t24;
t17 = t43 * t23 - t40 * t56;
t16 = -t40 * t23 - t43 * t56;
t15 = t38 * t57 + (t30 * t45 + t42 * t54) * t33;
t14 = -t38 * t56 + (t30 * t42 - t45 * t54) * t33;
t13 = t34 * t31 * t35 + (t34 * t53 + t55) * t30;
t11 = t44 * t15 + t41 * t20;
t10 = t41 * t15 - t44 * t20;
t9 = -t12 * t32 + t21 * t37;
t8 = t43 * t11 + t40 * t14;
t7 = -t40 * t11 + t43 * t14;
t6 = t13 * t45 + t46 * t42;
t5 = t13 * t42 - t46 * t45;
t4 = t9 * t41 + t6 * t44;
t3 = t6 * t41 - t9 * t44;
t2 = t4 * t43 + t5 * t40;
t1 = -t4 * t40 + t5 * t43;
t25 = [1, t39 ^ 2 + (t31 ^ 2 + t36 ^ 2) * t34 ^ 2, t12 ^ 2 + t13 ^ 2 + t21 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t39, t21 * t38 + (t12 * t35 + t13 * t30) * t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, t38 ^ 2 + (t30 ^ 2 + t35 ^ 2) * t33 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t44, t5 * t41, 0, 0, 0, 0, 0, -t1 * t44 + t3 * t52, t2 * t44 + t3 * t50; 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t14 * t44, t14 * t41, 0, 0, 0, 0, 0, t10 * t52 - t7 * t44, t10 * t50 + t8 * t44; 0, 0, 0, 0, t56, -t57, 0, 0, 0, 0, 0, t44 * t56, -t41 * t56, 0, 0, 0, 0, 0, -t16 * t44 + t22 * t52, t17 * t44 + t22 * t50; 0, 0, 0, 1, 0, 0, t28, t47, 0, 0, 0, pkin(4) * t60, pkin(4) * t61, t29 * t28, -0.2e1 * t28 * t51, t49 * t61, t40 * t47, t44 ^ 2, -0.2e1 * t18 * t44 + 0.2e1 * t28 * t58, 0.2e1 * t28 * pkin(10) * t43 + 0.2e1 * t19 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, -t3 * t43, t3 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t10 * t43, t10 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, 0, 0, 0, 0, -t22 * t43, t22 * t40; 0, 0, 0, 0, 0, 0, 0, 0, t41, t44, 0, -t41 * pkin(10), -t44 * pkin(10), t40 * t50 (-t27 + t29) * t41, -t48, -t49, 0, -pkin(10) * t50 + (-pkin(5) * t41 + pkin(11) * t44) * t40, pkin(11) * t49 + (t58 - t59) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t27, 0.2e1 * t51, 0, 0, 0, 0.2e1 * t59, -0.2e1 * pkin(5) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t52, -t44, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t40 * pkin(11), -t43 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t25;
