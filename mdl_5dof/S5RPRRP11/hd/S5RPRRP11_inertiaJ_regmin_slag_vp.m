% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:32
% EndTime: 2019-12-31 18:54:33
% DurationCPUTime: 0.35s
% Computational Cost: add. (395->60), mult. (772->120), div. (0->0), fcn. (869->6), ass. (0->47)
t27 = sin(pkin(8));
t28 = cos(pkin(8));
t30 = sin(qJ(3));
t48 = cos(qJ(3));
t17 = t48 * t27 + t30 * t28;
t55 = -0.2e1 * t17;
t22 = -t28 * pkin(2) - pkin(1);
t54 = 0.2e1 * t22;
t29 = sin(qJ(4));
t53 = -0.2e1 * t29;
t45 = pkin(6) + qJ(2);
t19 = t45 * t27;
t20 = t45 * t28;
t11 = -t30 * t19 + t48 * t20;
t31 = cos(qJ(4));
t16 = t30 * t27 - t48 * t28;
t8 = t16 * pkin(3) - t17 * pkin(7) + t22;
t4 = t31 * t11 + t29 * t8;
t52 = pkin(7) * t16;
t51 = t16 * pkin(4);
t50 = t29 * pkin(7);
t49 = t31 * pkin(7);
t12 = t29 * t16;
t47 = t29 * t17;
t46 = t29 * t31;
t13 = t31 * t16;
t14 = t31 * t17;
t44 = t27 ^ 2 + t28 ^ 2;
t25 = t29 ^ 2;
t26 = t31 ^ 2;
t43 = t25 + t26;
t42 = t16 * qJ(5);
t41 = t16 * t55;
t40 = t29 * t11 - t31 * t8;
t39 = -pkin(3) * t17 - t52;
t1 = t42 + t4;
t2 = t40 - t51;
t38 = t1 * t31 + t2 * t29;
t37 = t1 * t29 - t2 * t31;
t35 = t31 * pkin(4) + t29 * qJ(5);
t18 = -pkin(3) - t35;
t36 = -t17 * t18 + t52;
t34 = pkin(4) * t29 - t31 * qJ(5);
t10 = t48 * t19 + t30 * t20;
t15 = t17 ^ 2;
t5 = t34 * t17 + t10;
t3 = [1, 0, 0, 0.2e1 * pkin(1) * t28, -0.2e1 * pkin(1) * t27, 0.2e1 * t44 * qJ(2), t44 * qJ(2) ^ 2 + pkin(1) ^ 2, t15, t41, 0, 0, 0, t16 * t54, t17 * t54, t26 * t15, -0.2e1 * t15 * t46, 0.2e1 * t16 * t14, t29 * t41, t16 ^ 2, 0.2e1 * t10 * t47 - 0.2e1 * t16 * t40, 0.2e1 * t10 * t14 - 0.2e1 * t4 * t16, -0.2e1 * t2 * t16 + 0.2e1 * t5 * t47, t37 * t55, 0.2e1 * t1 * t16 - 0.2e1 * t5 * t14, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -t28, t27, 0, -pkin(1), 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, t13, -t12, t13, -t43 * t17, t12, t37; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t10, -t11, t29 * t14, (-t25 + t26) * t17, t12, t13, 0, -t10 * t31 + t39 * t29, t10 * t29 + t39 * t31, -t36 * t29 - t5 * t31, t38, -t5 * t29 + t36 * t31, t38 * pkin(7) + t5 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t25, 0.2e1 * t46, 0, 0, 0, 0.2e1 * pkin(3) * t31, pkin(3) * t53, -0.2e1 * t18 * t31, 0.2e1 * t43 * pkin(7), t18 * t53, t43 * pkin(7) ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t47, t16, -t40, -t4, -t40 + 0.2e1 * t51, -t35 * t17, 0.2e1 * t42 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, t31, 0, t29, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t31, 0, -t50, -t49, -t50, -t34, t49, -t34 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t14, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
