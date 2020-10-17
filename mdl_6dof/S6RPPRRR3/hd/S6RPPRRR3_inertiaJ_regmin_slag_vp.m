% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:31:48
% EndTime: 2019-05-05 15:31:50
% DurationCPUTime: 0.51s
% Computational Cost: add. (263->78), mult. (531->132), div. (0->0), fcn. (606->8), ass. (0->59)
t38 = sin(qJ(6));
t39 = sin(qJ(5));
t41 = cos(qJ(6));
t42 = cos(qJ(5));
t19 = t38 * t42 + t41 * t39;
t43 = cos(qJ(4));
t12 = t43 * t19;
t65 = -0.2e1 * t12;
t36 = sin(pkin(10));
t25 = t36 * pkin(1) + qJ(3);
t64 = 0.2e1 * t25;
t31 = -t42 * pkin(5) - pkin(4);
t63 = 0.2e1 * t31;
t40 = sin(qJ(4));
t62 = 0.2e1 * t40;
t61 = 0.2e1 * t42;
t60 = pkin(8) + pkin(9);
t59 = t38 * pkin(5);
t58 = t41 * pkin(5);
t17 = t40 * pkin(4) - t43 * pkin(8) + t25;
t37 = cos(pkin(10));
t27 = -t37 * pkin(1) - pkin(2);
t24 = -pkin(7) + t27;
t51 = t42 * t40;
t46 = t24 * t51;
t5 = t46 + (-pkin(9) * t43 + t17) * t39;
t57 = t41 * t5;
t16 = t19 * t40;
t56 = t24 * t39;
t28 = t39 * t40;
t55 = t39 * t42;
t54 = t39 * t43;
t18 = t38 * t39 - t41 * t42;
t53 = t40 * t18;
t52 = t40 * t24;
t30 = t42 * t43;
t50 = t43 * t24;
t49 = t43 * t40;
t33 = t40 ^ 2;
t35 = t43 ^ 2;
t48 = -t33 - t35;
t47 = -0.2e1 * t49;
t10 = t42 * t17;
t4 = -pkin(9) * t30 + t10 + (pkin(5) - t56) * t40;
t1 = -t38 * t5 + t41 * t4;
t45 = -pkin(4) * t43 - pkin(8) * t40;
t34 = t42 ^ 2;
t32 = t39 ^ 2;
t21 = t60 * t42;
t20 = t60 * t39;
t15 = (pkin(5) * t39 - t24) * t43;
t14 = t41 * t30 - t38 * t54;
t13 = -t38 * t28 + t41 * t51;
t9 = -t38 * t20 + t41 * t21;
t8 = -t41 * t20 - t38 * t21;
t7 = t39 * t17 + t46;
t6 = -t39 * t52 + t10;
t2 = t38 * t4 + t57;
t3 = [1, 0, 0 (t36 ^ 2 + t37 ^ 2) * pkin(1) ^ 2, 0.2e1 * t27, t64, t25 ^ 2 + t27 ^ 2, t35, t47, 0, 0, 0, t25 * t62, t43 * t64, t34 * t35, -0.2e1 * t35 * t55, t49 * t61, t39 * t47, t33, -0.2e1 * t35 * t56 + 0.2e1 * t6 * t40, -0.2e1 * t35 * t24 * t42 - 0.2e1 * t7 * t40, t14 ^ 2, t14 * t65, t14 * t62, t40 * t65, t33, 0.2e1 * t1 * t40 + 0.2e1 * t15 * t12, 0.2e1 * t15 * t14 - 0.2e1 * t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t39, t48 * t42, 0, 0, 0, 0, 0, -t43 * t12 - t16 * t40, -t13 * t40 - t43 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t40, 0, t50, -t52, t39 * t30 (-t32 + t34) * t43, t28, t51, 0, t45 * t39 + t42 * t50, -t39 * t50 + t45 * t42, t14 * t19, -t19 * t12 - t14 * t18, t16, -t53, 0, t31 * t12 + t15 * t18 + t8 * t40, t31 * t14 + t15 * t19 - t9 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t43, 0, 0, 0, 0, 0, -t51, t28, 0, 0, 0, 0, 0, t53, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t40, 0, 0, 0, 0, 0, t30, -t54, 0, 0, 0, 0, 0, -t43 * t18, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t32, 0.2e1 * t55, 0, 0, 0, pkin(4) * t61, -0.2e1 * pkin(4) * t39, t19 ^ 2, -0.2e1 * t19 * t18, 0, 0, 0, t18 * t63, t19 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t54, t40, t6, -t7, 0, 0, t14, -t12, t40, t40 * t58 + t1, -t57 + (-t40 * pkin(5) - t4) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t30, 0, 0, 0, 0, 0, -t12, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t51, 0, 0, 0, 0, 0, -t16, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t42, 0, -t39 * pkin(8), -t42 * pkin(8), 0, 0, t19, -t18, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t12, t40, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t58, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
