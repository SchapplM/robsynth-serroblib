% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:41
% EndTime: 2019-12-05 16:19:43
% DurationCPUTime: 0.44s
% Computational Cost: add. (345->52), mult. (701->99), div. (0->0), fcn. (828->6), ass. (0->37)
t29 = sin(pkin(9));
t30 = cos(pkin(9));
t32 = sin(qJ(3));
t33 = cos(qJ(3));
t18 = t29 * t32 - t30 * t33;
t26 = -t33 * pkin(3) - pkin(2);
t12 = t18 * pkin(4) + t26;
t43 = 0.2e1 * t12;
t42 = 0.2e1 * t26;
t41 = 0.2e1 * t33;
t40 = t29 * pkin(3);
t39 = t30 * pkin(3);
t38 = cos(qJ(5));
t37 = -qJ(4) - pkin(6);
t27 = t32 ^ 2;
t28 = t33 ^ 2;
t36 = t27 + t28;
t22 = t37 * t32;
t23 = t37 * t33;
t10 = t30 * t22 + t29 * t23;
t11 = t29 * t22 - t30 * t23;
t31 = sin(qJ(5));
t25 = pkin(4) + t39;
t20 = t29 * t33 + t30 * t32;
t17 = t20 ^ 2;
t16 = t18 ^ 2;
t14 = t31 * t25 + t38 * t40;
t13 = t38 * t25 - t31 * t40;
t9 = -t31 * t18 + t38 * t20;
t7 = t38 * t18 + t31 * t20;
t6 = t9 ^ 2;
t5 = t7 ^ 2;
t4 = -t18 * pkin(7) + t11;
t3 = -t20 * pkin(7) + t10;
t2 = t31 * t3 + t38 * t4;
t1 = t38 * t3 - t31 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 + t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18 * t10 + t20 * t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t1 + t9 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t27, t32 * t41, 0, t28, 0, 0, pkin(2) * t41, -0.2e1 * pkin(2) * t32, 0.2e1 * t36 * pkin(6), t36 * pkin(6) ^ 2 + pkin(2) ^ 2, t17, -0.2e1 * t20 * t18, 0, t16, 0, 0, t18 * t42, t20 * t42, -0.2e1 * t10 * t20 - 0.2e1 * t11 * t18, t10 ^ 2 + t11 ^ 2 + t26 ^ 2, t6, -0.2e1 * t9 * t7, 0, t5, 0, 0, t7 * t43, t9 * t43, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t20, 0, (-t18 * t30 + t20 * t29) * pkin(3), 0, 0, 0, 0, 0, 0, -t7, -t9, 0, -t7 * t13 + t9 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, t33, 0, -t32 * pkin(6), -t33 * pkin(6), 0, 0, 0, 0, t20, 0, -t18, 0, t10, -t11, (-t18 * t29 - t20 * t30) * pkin(3), (t10 * t30 + t11 * t29) * pkin(3), 0, 0, t9, 0, -t7, 0, t1, -t2, -t13 * t9 - t14 * t7, t1 * t13 + t2 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t39, -0.2e1 * t40, 0, (t29 ^ 2 + t30 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t13, -0.2e1 * t14, 0, t13 ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t20, 0, t26, 0, 0, 0, 0, 0, 0, t7, t9, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t13, -t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
