% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:04
% EndTime: 2019-12-31 17:10:06
% DurationCPUTime: 0.40s
% Computational Cost: add. (291->68), mult. (661->150), div. (0->0), fcn. (671->6), ass. (0->47)
t31 = sin(pkin(7));
t32 = cos(pkin(7));
t33 = sin(qJ(4));
t52 = cos(qJ(4));
t56 = -t33 * t31 + t52 * t32;
t24 = -t32 * pkin(3) - pkin(2);
t55 = 0.2e1 * t24;
t35 = cos(qJ(2));
t54 = -0.2e1 * t35;
t34 = sin(qJ(2));
t29 = t34 ^ 2;
t53 = t29 * pkin(5);
t26 = t34 * pkin(5);
t51 = t31 * t32;
t50 = t31 * t34;
t49 = t31 * t35;
t48 = t32 * t34;
t47 = t32 * t35;
t45 = t34 * t35;
t44 = pkin(6) + qJ(3);
t18 = -t35 * pkin(2) - t34 * qJ(3) - pkin(1);
t8 = pkin(5) * t47 + t31 * t18;
t27 = t31 ^ 2;
t28 = t32 ^ 2;
t43 = t27 + t28;
t42 = 0.2e1 * t45;
t41 = t31 * t48;
t13 = t32 * t18;
t7 = -pkin(5) * t49 + t13;
t39 = -t7 * t31 + t8 * t32;
t38 = -pkin(2) * t34 + qJ(3) * t35;
t16 = t52 * t31 + t33 * t32;
t37 = pkin(5) ^ 2;
t30 = t35 ^ 2;
t25 = t29 * t37;
t20 = t44 * t32;
t19 = t44 * t31;
t17 = pkin(3) * t50 + t26;
t11 = t56 * t34;
t9 = t16 * t34;
t6 = -t33 * t19 + t52 * t20;
t5 = -t52 * t19 - t33 * t20;
t4 = -pkin(6) * t50 + t8;
t3 = -pkin(6) * t48 + t13 + (-pkin(5) * t31 - pkin(3)) * t35;
t2 = t33 * t3 + t52 * t4;
t1 = t52 * t3 - t33 * t4;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t29, t42, 0, t30, 0, 0, 0.2e1 * pkin(1) * t35, -0.2e1 * pkin(1) * t34, 0.2e1 * (t29 + t30) * pkin(5), pkin(1) ^ 2 + t30 * t37 + t25, t28 * t29, -0.2e1 * t29 * t51, -0.2e1 * t32 * t45, t27 * t29, t31 * t42, t30, 0.2e1 * t31 * t53 - 0.2e1 * t7 * t35, 0.2e1 * t32 * t53 + 0.2e1 * t8 * t35, 0.2e1 * (-t31 * t8 - t32 * t7) * t34, t7 ^ 2 + t8 ^ 2 + t25, t11 ^ 2, -0.2e1 * t11 * t9, t11 * t54, t9 ^ 2, -t9 * t54, t30, -0.2e1 * t1 * t35 + 0.2e1 * t17 * t9, 0.2e1 * t17 * t11 + 0.2e1 * t2 * t35, -0.2e1 * t1 * t11 - 0.2e1 * t2 * t9, t1 ^ 2 + t17 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t35, 0, -t26, -t35 * pkin(5), 0, 0, t41, (-t27 + t28) * t34, -t49, -t41, -t47, 0, -pkin(5) * t48 + t38 * t31, pkin(5) * t50 + t38 * t32, t39, -pkin(2) * t26 + t39 * qJ(3), t11 * t16, t11 * t56 - t16 * t9, -t16 * t35, -t9 * t56, -t56 * t35, 0, -t17 * t56 + t24 * t9 - t5 * t35, t24 * t11 + t17 * t16 + t6 * t35, -t1 * t16 - t5 * t11 + t2 * t56 - t6 * t9, t1 * t5 + t17 * t24 + t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t27, 0.2e1 * t51, 0, t28, 0, 0, 0.2e1 * pkin(2) * t32, -0.2e1 * pkin(2) * t31, 0.2e1 * t43 * qJ(3), t43 * qJ(3) ^ 2 + pkin(2) ^ 2, t16 ^ 2, 0.2e1 * t16 * t56, 0, t56 ^ 2, 0, 0, -t56 * t55, t16 * t55, -0.2e1 * t5 * t16 + 0.2e1 * t56 * t6, t24 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t48, 0, t26, 0, 0, 0, 0, 0, 0, t9, t11, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t31, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t56, t16, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, -t35, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t56, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t10;
