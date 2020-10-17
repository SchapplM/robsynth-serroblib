% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:53
% EndTime: 2019-12-31 17:59:55
% DurationCPUTime: 0.49s
% Computational Cost: add. (194->53), mult. (351->93), div. (0->0), fcn. (312->6), ass. (0->49)
t21 = sin(pkin(8));
t47 = t21 * pkin(1);
t9 = qJ(3) + t47;
t52 = t9 ^ 2;
t51 = 0.2e1 * t9;
t25 = cos(qJ(5));
t50 = 0.2e1 * t25;
t26 = cos(qJ(4));
t49 = 0.2e1 * t26;
t20 = t26 ^ 2;
t22 = cos(pkin(8));
t46 = t22 * pkin(1);
t11 = -pkin(2) - t46;
t8 = -pkin(6) + t11;
t48 = t20 * t8;
t24 = sin(qJ(4));
t45 = t24 * pkin(4);
t44 = t24 * t8;
t43 = t26 * pkin(4);
t42 = t26 * t8;
t23 = sin(qJ(5));
t17 = t23 ^ 2;
t41 = t17 * t26;
t13 = t23 * t24;
t40 = t23 * t25;
t39 = t23 * t26;
t38 = t25 * t24;
t16 = t25 * t26;
t37 = t26 * t24;
t19 = t25 ^ 2;
t36 = t17 + t19;
t18 = t24 ^ 2;
t35 = t18 + t20;
t34 = -0.2e1 * t37;
t33 = t23 * t16;
t32 = t36 * pkin(7);
t31 = t36 * t24;
t30 = -pkin(7) * t24 - t43;
t4 = -t26 * pkin(7) + t45 + t9;
t1 = -t8 * t13 + t25 * t4;
t2 = t23 * t4 + t8 * t38;
t29 = -t1 * t23 + t2 * t25;
t15 = t19 * t26;
t14 = t19 * t20;
t12 = t17 * t20;
t6 = t8 ^ 2;
t5 = t20 * t6;
t3 = t35 * t8;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t46, -0.2e1 * t47, 0, (t21 ^ 2 + t22 ^ 2) * pkin(1) ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t11, t51, t11 ^ 2 + t52, t20, t34, 0, t18, 0, 0, t24 * t51, t9 * t49, -0.2e1 * t3, t18 * t6 + t5 + t52, t14, -0.2e1 * t20 * t40, t37 * t50, t12, t23 * t34, t18, 0.2e1 * t1 * t24 - 0.2e1 * t23 * t48, -0.2e1 * t2 * t24 - 0.2e1 * t25 * t48, (-t1 * t25 - t2 * t23) * t49, t1 ^ 2 + t2 ^ 2 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t29 - t44) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 + t12 + t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t3, 0, 0, 0, 0, 0, 0, -t35 * t23, -t35 * t25, 0, t29 * t24 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t36) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t18 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t24, 0, t42, -t44, 0, 0, t33, t15 - t41, t13, -t33, t38, 0, t8 * t16 + t30 * t23, t30 * t25 - t8 * t39, t29, pkin(4) * t42 + t29 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t13, t15 + t41, t26 * t32 - t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t24, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t39, t31, pkin(7) * t31 + t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t17, 0.2e1 * t40, 0, t19, 0, 0, pkin(4) * t50, -0.2e1 * pkin(4) * t23, 0.2e1 * t32, t36 * pkin(7) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t39, t24, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t25, 0, -t23 * pkin(7), -t25 * pkin(7), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t7;
