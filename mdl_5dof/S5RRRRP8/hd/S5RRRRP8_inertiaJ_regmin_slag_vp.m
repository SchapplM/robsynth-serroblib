% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:08
% EndTime: 2019-12-31 22:02:10
% DurationCPUTime: 0.39s
% Computational Cost: add. (406->82), mult. (854->158), div. (0->0), fcn. (909->6), ass. (0->55)
t37 = sin(qJ(4));
t40 = cos(qJ(4));
t39 = sin(qJ(2));
t41 = cos(qJ(3));
t46 = t41 * t39;
t38 = sin(qJ(3));
t50 = t38 * t39;
t17 = -t37 * t50 + t40 * t46;
t59 = -0.2e1 * t17;
t30 = -t41 * pkin(3) - pkin(2);
t58 = 0.2e1 * t30;
t57 = -0.2e1 * t39;
t42 = cos(qJ(2));
t56 = 0.2e1 * t42;
t55 = pkin(7) + pkin(8);
t54 = pkin(2) * t41;
t53 = pkin(6) * t38;
t52 = t37 * pkin(3);
t51 = t42 * pkin(3);
t49 = t38 * t41;
t48 = t38 * t42;
t23 = -t42 * pkin(2) - t39 * pkin(7) - pkin(1);
t45 = t41 * t42;
t43 = pkin(6) * t45;
t10 = t43 + (-pkin(8) * t39 + t23) * t38;
t47 = t40 * t10;
t31 = t39 * pkin(6);
t22 = pkin(3) * t50 + t31;
t44 = t39 * t56;
t18 = t41 * t23;
t8 = -pkin(8) * t46 + t18 + (-pkin(3) - t53) * t42;
t3 = -t37 * t10 + t40 * t8;
t24 = t55 * t38;
t25 = t55 * t41;
t11 = -t40 * t24 - t37 * t25;
t4 = t37 * t8 + t47;
t12 = -t37 * t24 + t40 * t25;
t20 = t37 * t41 + t40 * t38;
t36 = t42 ^ 2;
t35 = t41 ^ 2;
t34 = t39 ^ 2;
t33 = t38 ^ 2;
t32 = t40 * pkin(3);
t29 = t32 + pkin(4);
t19 = t37 * t38 - t40 * t41;
t16 = t20 * t39;
t15 = t19 * pkin(4) + t30;
t14 = t38 * t23 + t43;
t13 = -pkin(6) * t48 + t18;
t9 = t16 * pkin(4) + t22;
t6 = -t19 * qJ(5) + t12;
t5 = -t20 * qJ(5) + t11;
t2 = -t16 * qJ(5) + t4;
t1 = -t42 * pkin(4) - t17 * qJ(5) + t3;
t7 = [1, 0, 0, t34, t44, 0, 0, 0, pkin(1) * t56, pkin(1) * t57, t35 * t34, -0.2e1 * t34 * t49, t45 * t57, t38 * t44, t36, -0.2e1 * t13 * t42 + 0.2e1 * t34 * t53, 0.2e1 * t34 * pkin(6) * t41 + 0.2e1 * t14 * t42, t17 ^ 2, t16 * t59, t42 * t59, t16 * t56, t36, 0.2e1 * t22 * t16 - 0.2e1 * t3 * t42, 0.2e1 * t22 * t17 + 0.2e1 * t4 * t42, -0.2e1 * t1 * t17 - 0.2e1 * t2 * t16, t1 ^ 2 + t2 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, t39, t42, 0, -t31, -t42 * pkin(6), t38 * t46, (-t33 + t35) * t39, -t48, -t45, 0, -pkin(6) * t46 + (-pkin(2) * t39 + pkin(7) * t42) * t38, pkin(7) * t45 + (t53 - t54) * t39, t17 * t20, -t20 * t16 - t17 * t19, -t20 * t42, t19 * t42, 0, -t11 * t42 + t30 * t16 + t22 * t19, t12 * t42 + t30 * t17 + t22 * t20, -t1 * t20 - t6 * t16 - t5 * t17 - t2 * t19, t1 * t5 + t9 * t15 + t2 * t6; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t33, 0.2e1 * t49, 0, 0, 0, 0.2e1 * t54, -0.2e1 * pkin(2) * t38, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t58, t20 * t58, -0.2e1 * t6 * t19 - 0.2e1 * t5 * t20, t15 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t50, -t42, t13, -t14, 0, 0, t17, -t16, -t42, -t40 * t51 + t3, -t47 + (-t8 + t51) * t37, -t16 * t52 - t29 * t17, t1 * t29 + t2 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * pkin(7), -t41 * pkin(7), 0, 0, t20, -t19, 0, t11, -t12, -t19 * t52 - t29 * t20, t5 * t29 + t6 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t32, -0.2e1 * t52, 0, t37 ^ 2 * pkin(3) ^ 2 + t29 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t42, t3, -t4, -pkin(4) * t17, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t11, -t12, -pkin(4) * t20, t5 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t32, -t52, 0, t29 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
