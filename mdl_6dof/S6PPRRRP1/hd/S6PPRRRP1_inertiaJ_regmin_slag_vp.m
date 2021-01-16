% Calculate minimal parameter regressor of joint inertia matrix for
% S6PPRRRP1
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
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:11
% EndTime: 2021-01-16 00:50:14
% DurationCPUTime: 0.56s
% Computational Cost: add. (466->105), mult. (1243->191), div. (0->0), fcn. (1549->12), ass. (0->63)
t41 = sin(qJ(4));
t64 = -0.2e1 * t41;
t44 = cos(qJ(4));
t63 = 0.2e1 * t44;
t43 = cos(qJ(5));
t62 = pkin(4) * t43;
t40 = sin(qJ(5));
t61 = pkin(9) * t40;
t60 = t40 * pkin(5);
t34 = sin(pkin(12));
t36 = sin(pkin(6));
t39 = cos(pkin(6));
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t37 = cos(pkin(12));
t38 = cos(pkin(7));
t55 = t37 * t38;
t35 = sin(pkin(7));
t57 = t35 * t42;
t12 = t39 * t57 + (t34 * t45 + t42 * t55) * t36;
t19 = -t35 * t36 * t37 + t38 * t39;
t6 = t12 * t41 - t19 * t44;
t59 = t6 * t43;
t20 = -t38 * t44 + t41 * t57;
t58 = t20 * t43;
t56 = t35 * t45;
t54 = t40 * t41;
t53 = t40 * t43;
t52 = t40 * t44;
t28 = t43 * t41;
t51 = t43 * t44;
t50 = -qJ(6) - pkin(10);
t49 = qJ(6) * t41;
t48 = t41 * t63;
t47 = pkin(9) * t51;
t24 = -pkin(4) * t44 - pkin(10) * t41 - pkin(3);
t22 = t43 * t24;
t46 = -t43 * t49 + t22;
t33 = t43 ^ 2;
t32 = t41 ^ 2;
t31 = t40 ^ 2;
t29 = -pkin(5) * t43 - pkin(4);
t26 = t50 * t43;
t25 = t50 * t40;
t23 = (pkin(9) + t60) * t41;
t21 = t38 * t41 + t44 * t57;
t18 = t20 * t40;
t17 = t24 * t40 + t47;
t16 = -pkin(9) * t52 + t22;
t15 = t47 + (t24 - t49) * t40;
t14 = t21 * t43 - t40 * t56;
t13 = -t21 * t40 - t43 * t56;
t11 = -t39 * t56 + (t34 * t42 - t45 * t55) * t36;
t10 = (-pkin(5) - t61) * t44 + t46;
t9 = t14 * t44 + t20 * t28;
t8 = -t13 * t44 + t20 * t54;
t7 = t12 * t44 + t19 * t41;
t5 = t6 * t40;
t4 = t11 * t40 + t43 * t7;
t3 = t11 * t43 - t40 * t7;
t2 = t28 * t6 + t4 * t44;
t1 = -t3 * t44 + t54 * t6;
t27 = [1, t39 ^ 2 + (t34 ^ 2 + t37 ^ 2) * t36 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t6 ^ 2; 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t3 + t14 * t4 + t20 * t6; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 ^ 2 + t14 ^ 2 + t20 ^ 2; 0, 0, 0, -t11, -t12, 0, 0, 0, 0, 0, -t11 * t44, t11 * t41, 0, 0, 0, 0, 0, t1, t2, t1, t2, (-t3 * t43 - t4 * t40) * t41, t10 * t3 + t15 * t4 + t23 * t6; 0, 0, 0, t56, -t57, 0, 0, 0, 0, 0, t44 * t56, -t41 * t56, 0, 0, 0, 0, 0, t8, t9, t8, t9, (-t13 * t43 - t14 * t40) * t41, t10 * t13 + t14 * t15 + t20 * t23; 0, 0, 1, 0, 0, t32, t48, 0, 0, 0, pkin(3) * t63, pkin(3) * t64, t33 * t32, -0.2e1 * t32 * t53, t51 * t64, t40 * t48, t44 ^ 2, -0.2e1 * t16 * t44 + 0.2e1 * t32 * t61, 0.2e1 * pkin(9) * t32 * t43 + 0.2e1 * t17 * t44, -0.2e1 * t10 * t44 + 0.2e1 * t23 * t54, 0.2e1 * t15 * t44 + 0.2e1 * t23 * t28, 0.2e1 * (-t10 * t43 - t15 * t40) * t41, t10 ^ 2 + t15 ^ 2 + t23 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t59, t5, -t59, t5, -t3 * t40 + t4 * t43, t25 * t3 - t26 * t4 + t29 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, 0, 0, 0, 0, -t58, t18, -t58, t18, -t13 * t40 + t14 * t43, t13 * t25 - t14 * t26 + t20 * t29; 0, 0, 0, 0, 0, 0, 0, t41, t44, 0, -t41 * pkin(9), -t44 * pkin(9), t40 * t28, (-t31 + t33) * t41, -t52, -t51, 0, -pkin(9) * t28 + (-pkin(4) * t41 + pkin(10) * t44) * t40, pkin(10) * t51 + (t61 - t62) * t41, -t23 * t43 - t25 * t44 + t29 * t54, t23 * t40 - t26 * t44 + t28 * t29, (-t25 * t41 + t15) * t43 + (t26 * t41 - t10) * t40, t10 * t25 - t15 * t26 + t23 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t31, 0.2e1 * t53, 0, 0, 0, 0.2e1 * t62, -0.2e1 * pkin(4) * t40, -0.2e1 * t29 * t43, 0.2e1 * t29 * t40, -0.2e1 * t25 * t40 - 0.2e1 * t26 * t43, t25 ^ 2 + t26 ^ 2 + t29 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, t3, -t4, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, t13, -t14, 0, t13 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t54, -t44, t16, -t17, (-0.2e1 * pkin(5) - t61) * t44 + t46, -t15, -pkin(5) * t28, t10 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t40 * pkin(10), -t43 * pkin(10), t25, t26, -t60, t25 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t28, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t40, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t27;
