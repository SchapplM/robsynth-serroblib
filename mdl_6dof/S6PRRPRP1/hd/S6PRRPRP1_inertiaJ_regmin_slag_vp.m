% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:38:48
% EndTime: 2021-01-16 02:38:51
% DurationCPUTime: 0.63s
% Computational Cost: add. (674->105), mult. (1422->193), div. (0->0), fcn. (1726->10), ass. (0->67)
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t40 = sin(pkin(6));
t55 = t40 * sin(qJ(2));
t56 = cos(pkin(6));
t23 = t56 * t43 + t45 * t55;
t39 = sin(pkin(11));
t41 = cos(pkin(11));
t48 = -t43 * t55 + t56 * t45;
t11 = t39 * t23 - t41 * t48;
t74 = t11 ^ 2;
t36 = -t45 * pkin(3) - pkin(2);
t73 = 0.2e1 * t36;
t42 = sin(qJ(5));
t72 = 0.2e1 * t42;
t44 = cos(qJ(5));
t71 = -0.2e1 * t44;
t70 = 0.2e1 * t45;
t27 = t39 * t43 - t41 * t45;
t69 = t27 * pkin(5);
t68 = t39 * pkin(3);
t67 = t41 * pkin(3);
t66 = t44 * pkin(5);
t65 = t11 * t44;
t46 = cos(qJ(2));
t64 = t40 * t46;
t20 = t42 * t27;
t28 = t39 * t45 + t41 * t43;
t63 = t42 * t28;
t62 = t42 * t44;
t60 = qJ(4) + pkin(8);
t31 = t60 * t45;
t54 = t60 * t43;
t18 = t41 * t31 - t39 * t54;
t61 = t44 * t18;
t22 = t44 * t28;
t37 = t42 ^ 2;
t38 = t44 ^ 2;
t59 = t37 + t38;
t58 = qJ(6) * t28;
t34 = pkin(9) + t68;
t57 = qJ(6) + t34;
t35 = -pkin(4) - t67;
t15 = t27 * pkin(4) - t28 * pkin(9) + t36;
t5 = t44 * t15 - t42 * t18;
t16 = t39 * t31 + t41 * t54;
t49 = -t44 * t58 + t5;
t3 = t49 + t69;
t4 = t61 + (t15 - t58) * t42;
t53 = t3 * t44 + t4 * t42;
t13 = t41 * t23 + t39 * t48;
t7 = -t42 * t13 - t44 * t64;
t8 = t44 * t13 - t42 * t64;
t52 = t8 * t42 + t7 * t44;
t24 = t57 * t42;
t25 = t57 * t44;
t51 = -t24 * t44 + t25 * t42;
t50 = -t27 * t34 + t28 * t35;
t30 = t35 - t66;
t26 = t28 ^ 2;
t21 = t44 * t27;
t10 = t11 * t42;
t9 = pkin(5) * t63 + t16;
t6 = t42 * t15 + t61;
t2 = t11 * t22 - t8 * t27;
t1 = t11 * t63 + t7 * t27;
t12 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 ^ 2 * t46 ^ 2 + t13 ^ 2 + t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t8 ^ 2 + t74; 0, 0, t64, -t55, 0, 0, 0, 0, 0, t45 * t64, -t43 * t64, -t27 * t64, -t28 * t64, t11 * t28 - t13 * t27, t11 * t16 + t13 * t18 - t36 * t64, 0, 0, 0, 0, 0, t1, t2, t1, t2, -t52 * t28, t11 * t9 + t7 * t3 + t8 * t4; 0, 1, 0, 0, t43 ^ 2, t43 * t70, 0, 0, 0, pkin(2) * t70, -0.2e1 * pkin(2) * t43, t27 * t73, t28 * t73, 0.2e1 * t16 * t28 - 0.2e1 * t18 * t27, t16 ^ 2 + t18 ^ 2 + t36 ^ 2, t38 * t26, -0.2e1 * t26 * t62, 0.2e1 * t27 * t22, -0.2e1 * t27 * t63, t27 ^ 2, 0.2e1 * t16 * t63 + 0.2e1 * t5 * t27, 0.2e1 * t16 * t22 - 0.2e1 * t6 * t27, 0.2e1 * t3 * t27 + 0.2e1 * t9 * t63, 0.2e1 * t9 * t22 - 0.2e1 * t4 * t27, -0.2e1 * t53 * t28, t3 ^ 2 + t4 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t23, -t11, -t13, 0, (-t11 * t41 + t13 * t39) * pkin(3), 0, 0, 0, 0, 0, -t65, t10, -t65, t10, -t7 * t42 + t8 * t44, t11 * t30 - t7 * t24 + t8 * t25; 0, 0, 0, 0, 0, 0, t43, t45, 0, -t43 * pkin(8), -t45 * pkin(8), -t16, -t18, (-t27 * t39 - t28 * t41) * pkin(3), (-t16 * t41 + t18 * t39) * pkin(3), t42 * t22, (-t37 + t38) * t28, t20, t21, 0, -t16 * t44 + t50 * t42, t16 * t42 + t50 * t44, -t24 * t27 + t30 * t63 - t9 * t44, t30 * t22 - t25 * t27 + t9 * t42, -t51 * t28 - t3 * t42 + t4 * t44, -t3 * t24 + t4 * t25 + t9 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t67, -0.2e1 * t68, 0, (t39 ^ 2 + t41 ^ 2) * pkin(3) ^ 2, t37, 0.2e1 * t62, 0, 0, 0, t35 * t71, t35 * t72, t30 * t71, t30 * t72, 0.2e1 * t24 * t42 + 0.2e1 * t25 * t44, t24 ^ 2 + t25 ^ 2 + t30 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t28, 0, t36, 0, 0, 0, 0, 0, t21, -t20, t21, -t20, -t59 * t28, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, t7, -t8, 0, t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t63, t27, t5, -t6, t49 + 0.2e1 * t69, -t4, -pkin(5) * t22, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t44, 0, -t42 * t34, -t44 * t34, -t24, -t25, -t42 * pkin(5), -t24 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t42, t44, -t42, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t22, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t42, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t12;
