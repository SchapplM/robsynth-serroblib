% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:24
% EndTime: 2019-12-31 21:09:25
% DurationCPUTime: 0.42s
% Computational Cost: add. (306->85), mult. (587->156), div. (0->0), fcn. (532->4), ass. (0->54)
t33 = cos(qJ(3));
t21 = qJ(4) * t33;
t31 = sin(qJ(3));
t63 = -pkin(3) * t31 + t21;
t29 = pkin(3) + qJ(5);
t45 = t31 * qJ(4);
t62 = -t29 * t33 - t45;
t61 = -0.2e1 * t31;
t32 = sin(qJ(2));
t60 = -0.2e1 * t32;
t59 = 0.2e1 * t32;
t34 = cos(qJ(2));
t58 = 0.2e1 * t34;
t57 = pkin(2) * t33;
t55 = pkin(6) * t31;
t54 = pkin(7) * t34;
t52 = t31 * t32;
t51 = t31 * t33;
t50 = t31 * t34;
t19 = t33 * t32;
t49 = t33 * t34;
t23 = t32 * pkin(6);
t48 = pkin(3) * t52 + t23;
t13 = -t34 * pkin(2) - t32 * pkin(7) - pkin(1);
t47 = pkin(6) * t50 - t33 * t13;
t7 = pkin(6) * t49 + t31 * t13;
t26 = t31 ^ 2;
t28 = t33 ^ 2;
t46 = t26 + t28;
t44 = t34 * qJ(4);
t43 = t32 * t58;
t25 = t34 * pkin(3);
t5 = t25 + t47;
t4 = t44 - t7;
t42 = t5 * t31 - t4 * t33;
t40 = -t33 * pkin(3) - t45;
t12 = -pkin(2) + t40;
t41 = -t12 * t32 - t54;
t39 = -pkin(4) * t52 + t7;
t38 = -pkin(4) * t19 - t5;
t36 = qJ(4) ^ 2;
t35 = 0.2e1 * qJ(4);
t27 = t32 ^ 2;
t24 = t33 * pkin(7);
t22 = t31 * pkin(7);
t20 = -0.2e1 * t44;
t15 = t33 * pkin(4) + t24;
t14 = t31 * pkin(4) + t22;
t9 = -pkin(2) + t62;
t8 = -t32 * t21 + t48;
t3 = (qJ(5) * t31 - t21) * t32 + t48;
t2 = t39 - t44;
t1 = t34 * qJ(5) - t38;
t6 = [1, 0, 0, t27, t43, 0, 0, 0, pkin(1) * t58, pkin(1) * t60, t28 * t27, -0.2e1 * t27 * t51, t49 * t60, t31 * t43, t34 ^ 2, 0.2e1 * t27 * t55 + 0.2e1 * t34 * t47, 0.2e1 * t27 * pkin(6) * t33 + 0.2e1 * t7 * t34, (t31 * t4 + t33 * t5) * t59, -0.2e1 * t5 * t34 - 0.2e1 * t8 * t52, -0.2e1 * t8 * t19 + 0.2e1 * t4 * t34, t4 ^ 2 + t5 ^ 2 + t8 ^ 2, (t1 * t33 - t2 * t31) * t59, -0.2e1 * t3 * t19 - 0.2e1 * t2 * t34, 0.2e1 * t1 * t34 + 0.2e1 * t3 * t52, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t32, t34, 0, -t23, -t34 * pkin(6), t31 * t19, (-t26 + t28) * t32, -t50, -t49, 0, -pkin(6) * t19 + (-pkin(2) * t32 + t54) * t31, pkin(7) * t49 + (t55 - t57) * t32, t42, t41 * t31 + t8 * t33, -t8 * t31 + t41 * t33, t42 * pkin(7) + t8 * t12, (t14 * t32 + t2) * t33 + (-t15 * t32 + t1) * t31, -t15 * t34 - t9 * t19 - t3 * t31, t14 * t34 - t3 * t33 + t9 * t52, t1 * t14 + t2 * t15 + t3 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t26, 0.2e1 * t51, 0, 0, 0, 0.2e1 * t57, pkin(2) * t61, 0.2e1 * t46 * pkin(7), 0.2e1 * t12 * t33, t12 * t61, t46 * pkin(7) ^ 2 + t12 ^ 2, 0.2e1 * t14 * t31 + 0.2e1 * t15 * t33, t9 * t61, -0.2e1 * t9 * t33, t14 ^ 2 + t15 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t52, -t34, -t47, -t7, t40 * t32, 0.2e1 * t25 + t47, t20 + t7, -t5 * pkin(3) - t4 * qJ(4), t62 * t32, t20 + t39, (-qJ(5) - t29) * t34 + t38, t2 * qJ(4) - t1 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t33, 0, -t22, -t24, t63, t22, t24, t63 * pkin(7), -t29 * t31 + t21, t15, -t14, t15 * qJ(4) - t14 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t35, pkin(3) ^ 2 + t36, 0, t35, 0.2e1 * t29, t29 ^ 2 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t34, 0, t5, t19, 0, t34, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, t22, t31, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, -1, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t34, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
