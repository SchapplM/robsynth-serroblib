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
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-16 00:21:20
% EndTime: 2021-01-16 00:21:23
% DurationCPUTime: 0.48s
% Computational Cost: add. (559->96), mult. (1159->179), div. (0->0), fcn. (1238->6), ass. (0->60)
t39 = sin(qJ(4));
t40 = sin(qJ(3));
t42 = cos(qJ(4));
t43 = cos(qJ(3));
t19 = t39 * t40 - t42 * t43;
t31 = -pkin(3) * t43 - pkin(2);
t15 = pkin(4) * t19 + t31;
t64 = 0.2e1 * t15;
t41 = sin(qJ(2));
t50 = t43 * t41;
t53 = t40 * t41;
t17 = -t39 * t53 + t42 * t50;
t63 = -0.2e1 * t17;
t62 = 0.2e1 * t31;
t61 = -0.2e1 * t41;
t44 = cos(qJ(2));
t60 = 0.2e1 * t44;
t59 = pkin(7) + pkin(8);
t58 = pkin(2) * t43;
t57 = pkin(6) * t40;
t56 = t39 * pkin(3);
t55 = t44 * pkin(3);
t54 = t44 * pkin(4);
t52 = t40 * t43;
t51 = t40 * t44;
t49 = t43 * t44;
t33 = t41 * pkin(6);
t22 = pkin(3) * t53 + t33;
t48 = t41 * t60;
t47 = pkin(6) * t49;
t23 = -pkin(2) * t44 - pkin(7) * t41 - pkin(1);
t10 = t47 + (-pkin(8) * t41 + t23) * t40;
t18 = t43 * t23;
t8 = -pkin(8) * t50 + t18 + (-pkin(3) - t57) * t44;
t3 = -t39 * t10 + t42 * t8;
t24 = t59 * t40;
t25 = t59 * t43;
t11 = -t42 * t24 - t25 * t39;
t4 = t10 * t42 + t39 * t8;
t12 = -t24 * t39 + t25 * t42;
t20 = t39 * t43 + t40 * t42;
t46 = -t17 * qJ(5) + t3;
t16 = t20 * t41;
t2 = -qJ(5) * t16 + t4;
t45 = 0.2e1 * pkin(4);
t38 = t44 ^ 2;
t37 = t43 ^ 2;
t36 = t41 ^ 2;
t35 = t40 ^ 2;
t34 = t42 * pkin(3);
t32 = -0.2e1 * t56;
t30 = t34 + pkin(4);
t27 = t39 * t55;
t14 = t23 * t40 + t47;
t13 = -pkin(6) * t51 + t18;
t9 = pkin(4) * t16 + t22;
t6 = -qJ(5) * t19 + t12;
t5 = -qJ(5) * t20 + t11;
t1 = t46 - t54;
t7 = [1, 0, 0, t36, t48, 0, 0, 0, pkin(1) * t60, pkin(1) * t61, t37 * t36, -0.2e1 * t36 * t52, t49 * t61, t40 * t48, t38, -0.2e1 * t13 * t44 + 0.2e1 * t36 * t57, 0.2e1 * pkin(6) * t36 * t43 + 0.2e1 * t14 * t44, t17 ^ 2, t16 * t63, t44 * t63, t16 * t60, t38, 0.2e1 * t16 * t22 - 0.2e1 * t3 * t44, 0.2e1 * t17 * t22 + 0.2e1 * t4 * t44, -0.2e1 * t1 * t44 + 0.2e1 * t16 * t9, 0.2e1 * t17 * t9 + 0.2e1 * t2 * t44, -0.2e1 * t1 * t17 - 0.2e1 * t16 * t2, t1 ^ 2 + t2 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, t41, t44, 0, -t33, -t44 * pkin(6), t40 * t50, (-t35 + t37) * t41, -t51, -t49, 0, -pkin(6) * t50 + (-pkin(2) * t41 + pkin(7) * t44) * t40, pkin(7) * t49 + (t57 - t58) * t41, t17 * t20, -t16 * t20 - t17 * t19, -t20 * t44, t19 * t44, 0, -t11 * t44 + t16 * t31 + t19 * t22, t12 * t44 + t17 * t31 + t20 * t22, t15 * t16 + t19 * t9 - t44 * t5, t15 * t17 + t20 * t9 + t44 * t6, -t1 * t20 - t16 * t6 - t17 * t5 - t19 * t2, t1 * t5 + t15 * t9 + t2 * t6; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t35, 0.2e1 * t52, 0, 0, 0, 0.2e1 * t58, -0.2e1 * pkin(2) * t40, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t62, t20 * t62, t19 * t64, t20 * t64, -0.2e1 * t19 * t6 - 0.2e1 * t20 * t5, t15 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t53, -t44, t13, -t14, 0, 0, t17, -t16, -t44, -t42 * t55 + t3, t27 - t4, (-pkin(4) - t30) * t44 + t46, -t2 + t27, -t16 * t56 - t17 * t30, t1 * t30 + t2 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t40 * pkin(7), -t43 * pkin(7), 0, 0, t20, -t19, 0, t11, -t12, t5, -t6, -t19 * t56 - t20 * t30, t30 * t5 + t56 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, t32, 0.2e1 * t30, t32, 0, pkin(3) ^ 2 * t39 ^ 2 + t30 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t44, t3, -t4, t46 - 0.2e1 * t54, -t2, -t17 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t11, -t12, t5, -t6, -t20 * pkin(4), t5 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t56, t45 + t34, -t56, 0, t30 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t45, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
