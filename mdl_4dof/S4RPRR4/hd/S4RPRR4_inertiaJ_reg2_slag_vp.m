% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:35
% EndTime: 2019-12-31 16:50:36
% DurationCPUTime: 0.25s
% Computational Cost: add. (127->43), mult. (272->85), div. (0->0), fcn. (239->6), ass. (0->41)
t18 = sin(qJ(3));
t42 = 0.2e1 * t18;
t19 = cos(qJ(4));
t41 = pkin(3) * t19;
t15 = sin(pkin(7));
t40 = t15 * pkin(1);
t16 = cos(pkin(7));
t39 = t16 * pkin(1);
t17 = sin(qJ(4));
t6 = pkin(5) + t40;
t38 = t17 * t6;
t37 = t18 * t6;
t20 = cos(qJ(3));
t36 = t20 * pkin(3);
t35 = t20 * t6;
t11 = t17 ^ 2;
t34 = t11 * t18;
t33 = t17 * t18;
t32 = t17 * t19;
t31 = t17 * t20;
t30 = t19 * t18;
t29 = t19 * t20;
t13 = t19 ^ 2;
t28 = t11 + t13;
t12 = t18 ^ 2;
t14 = t20 ^ 2;
t27 = t12 + t14;
t26 = t20 * t42;
t25 = t17 * t30;
t7 = -pkin(2) - t39;
t24 = t28 * pkin(6);
t3 = -t18 * pkin(6) - t36 + t7;
t1 = t19 * t3 - t6 * t31;
t2 = t17 * t3 + t6 * t29;
t23 = -t1 * t17 + t2 * t19;
t10 = t13 * t18;
t9 = t13 * t12;
t8 = t11 * t12;
t5 = t6 ^ 2;
t4 = t12 * t5;
t21 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t39, -0.2e1 * t40, 0, (t15 ^ 2 + t16 ^ 2) * pkin(1) ^ 2, t12, t26, 0, t14, 0, 0, -0.2e1 * t7 * t20, t7 * t42, 0.2e1 * t27 * t6, t14 * t5 + t7 ^ 2 + t4, t9, -0.2e1 * t12 * t32, -0.2e1 * t18 * t29, t8, t17 * t26, t14, -0.2e1 * t1 * t20 + 0.2e1 * t12 * t38, 0.2e1 * t12 * t6 * t19 + 0.2e1 * t2 * t20, (-t1 * t19 - t17 * t2) * t42, t1 ^ 2 + t2 ^ 2 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t23 - t35) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 + t8 + t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t20, 0, -t37, -t35, 0, 0, t25, t10 - t34, -t31, -t25, -t29, 0, -t6 * t30 + (-pkin(3) * t18 + pkin(6) * t20) * t17, pkin(6) * t29 + (t38 - t41) * t18, t23, -pkin(3) * t37 + t23 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t31, t10 + t34, t18 * t24 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t11, 0.2e1 * t32, 0, t13, 0, 0, 0.2e1 * t41, -0.2e1 * pkin(3) * t17, 0.2e1 * t24, t28 * pkin(6) ^ 2 + pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, -t33, -t20, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, t19, 0, -t17 * pkin(6), -t19 * pkin(6), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t21;
