% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:45:58
% DurationCPUTime: 0.75s
% Computational Cost: add. (677->91), mult. (1869->171), div. (0->0), fcn. (2038->8), ass. (0->56)
t38 = sin(pkin(5));
t41 = sin(qJ(4));
t42 = sin(qJ(3));
t44 = cos(qJ(4));
t45 = cos(qJ(3));
t20 = (-t41 * t42 + t44 * t45) * t38;
t22 = (t41 * t45 + t42 * t44) * t38;
t40 = sin(qJ(5));
t43 = cos(qJ(5));
t10 = -t43 * t20 + t40 * t22;
t63 = -0.2e1 * t10;
t39 = cos(pkin(5));
t62 = 0.2e1 * t39;
t61 = pkin(7) + pkin(8);
t60 = pkin(2) * t39;
t59 = pkin(2) * t42;
t58 = t39 * pkin(3);
t57 = t39 * pkin(4);
t56 = t40 * pkin(4);
t55 = t41 * pkin(3);
t34 = t43 * pkin(4);
t29 = t45 * t60;
t52 = t38 * t42;
t15 = -t61 * t52 + t29 + t58;
t32 = t38 * t45;
t49 = t39 * t59;
t16 = t61 * t32 + t49;
t51 = t44 * t16;
t7 = t41 * t15 + t51;
t5 = t20 * pkin(9) + t7;
t54 = t43 * t5;
t35 = t44 * pkin(3);
t36 = t38 ^ 2;
t53 = t36 * t45;
t50 = t38 * t62;
t48 = t43 * t55;
t6 = t44 * t15 - t41 * t16;
t4 = -t22 * pkin(9) + t57 + t6;
t1 = t43 * t4 - t40 * t5;
t33 = t35 + pkin(4);
t23 = t43 * t33 - t40 * t55;
t27 = (-pkin(3) * t45 - pkin(2)) * t38;
t2 = t40 * t4 + t54;
t37 = t39 ^ 2;
t31 = t36 * t45 ^ 2;
t30 = t36 * t42 ^ 2;
t26 = pkin(7) * t32 + t49;
t25 = -pkin(7) * t52 + t29;
t24 = t40 * t33 + t48;
t19 = t22 ^ 2;
t18 = t20 ^ 2;
t13 = -t20 * pkin(4) + t27;
t12 = t40 * t20 + t43 * t22;
t9 = t12 ^ 2;
t8 = t10 ^ 2;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 + t31 + t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 + t18 + t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 + t8 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t25 * t45 + t26 * t42 - t60) * t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t6 + t22 * t7 + t39 * t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t1 + t12 * t2 + t39 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t30, 0.2e1 * t42 * t53, t42 * t50, t31, t45 * t50, t37, 0.2e1 * pkin(2) * t53 + 0.2e1 * t25 * t39, -0.2e1 * t26 * t39 - 0.2e1 * t36 * t59, 0.2e1 * (-t25 * t42 + t26 * t45) * t38, t36 * pkin(2) ^ 2 + t25 ^ 2 + t26 ^ 2, t19, 0.2e1 * t22 * t20, t22 * t62, t18, t20 * t62, t37, -0.2e1 * t27 * t20 + 0.2e1 * t6 * t39, 0.2e1 * t27 * t22 - 0.2e1 * t7 * t39, 0.2e1 * t7 * t20 - 0.2e1 * t6 * t22, t27 ^ 2 + t6 ^ 2 + t7 ^ 2, t9, t12 * t63, t12 * t62, t8, t39 * t63, t37, 0.2e1 * t1 * t39 + 0.2e1 * t13 * t10, 0.2e1 * t13 * t12 - 0.2e1 * t2 * t39, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t13 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t52, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t22, 0, (t20 * t44 + t22 * t41) * pkin(3), 0, 0, 0, 0, 0, 0, -t10, -t12, 0, -t10 * t23 + t12 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t32, t39, t25, -t26, 0, 0, 0, 0, t22, 0, t20, t39, t39 * t35 + t6, -t51 + (-t15 - t58) * t41, (t20 * t41 - t22 * t44) * pkin(3), (t41 * t7 + t44 * t6) * pkin(3), 0, 0, t12, 0, -t10, t39, t23 * t39 + t1, -t24 * t39 - t2, -t24 * t10 - t23 * t12, t1 * t23 + t2 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t55, 0, (t41 ^ 2 + t44 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t23, -0.2e1 * t24, 0, t23 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t22, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t12, 0, (-t10 * t43 + t12 * t40) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t20, t39, t6, -t7, 0, 0, 0, 0, t12, 0, -t10, t39, t39 * t34 + t1, -t54 + (-t4 - t57) * t40, (-t10 * t40 - t12 * t43) * pkin(4), (t1 * t43 + t2 * t40) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t55, 0, 0, 0, 0, 0, 0, 0, 1, t23 + t34, -t48 + (-pkin(4) - t33) * t40, 0, (t23 * t43 + t24 * t40) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t56, 0, (t40 ^ 2 + t43 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, t39, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t23, -t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t56, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
