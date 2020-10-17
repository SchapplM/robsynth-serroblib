% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:39:17
% EndTime: 2020-01-03 11:39:20
% DurationCPUTime: 0.48s
% Computational Cost: add. (431->55), mult. (796->107), div. (0->0), fcn. (917->8), ass. (0->43)
t32 = sin(pkin(9));
t34 = cos(pkin(9));
t37 = sin(qJ(3));
t38 = cos(qJ(3));
t21 = t32 * t37 - t34 * t38;
t35 = cos(pkin(8));
t44 = t35 * pkin(1);
t29 = -pkin(2) - t44;
t24 = -t38 * pkin(3) + t29;
t12 = t21 * pkin(4) + t24;
t50 = 0.2e1 * t12;
t49 = 0.2e1 * t24;
t48 = 0.2e1 * t37;
t47 = t32 * pkin(3);
t33 = sin(pkin(8));
t46 = t33 * pkin(1);
t45 = t34 * pkin(3);
t43 = cos(qJ(5));
t30 = t37 ^ 2;
t31 = t38 ^ 2;
t42 = t30 + t31;
t27 = pkin(6) + t46;
t41 = qJ(4) + t27;
t17 = t41 * t37;
t18 = t41 * t38;
t5 = -t34 * t17 - t32 * t18;
t6 = -t32 * t17 + t34 * t18;
t36 = sin(qJ(5));
t28 = pkin(4) + t45;
t23 = t32 * t38 + t34 * t37;
t20 = t23 ^ 2;
t19 = t21 ^ 2;
t15 = t36 * t28 + t43 * t47;
t14 = t43 * t28 - t36 * t47;
t11 = -t36 * t21 + t43 * t23;
t9 = t43 * t21 + t36 * t23;
t8 = t11 ^ 2;
t7 = t9 ^ 2;
t4 = -t21 * pkin(7) + t6;
t3 = -t23 * pkin(7) + t5;
t2 = t36 * t3 + t43 * t4;
t1 = t43 * t3 - t36 * t4;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t44, -0.2e1 * t46, 0, (t33 ^ 2 + t35 ^ 2) * pkin(1) ^ 2, t30, t38 * t48, 0, t31, 0, 0, -0.2e1 * t29 * t38, t29 * t48, 0.2e1 * t42 * t27, t27 ^ 2 * t42 + t29 ^ 2, t20, -0.2e1 * t23 * t21, 0, t19, 0, 0, t21 * t49, t23 * t49, -0.2e1 * t6 * t21 - 0.2e1 * t5 * t23, t24 ^ 2 + t5 ^ 2 + t6 ^ 2, t8, -0.2e1 * t11 * t9, 0, t7, 0, 0, t9 * t50, t11 * t50, -0.2e1 * t1 * t11 - 0.2e1 * t2 * t9, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t21 + t6 * t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t9 + t2 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 + t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t38, 0, -t37 * t27, -t38 * t27, 0, 0, 0, 0, t23, 0, -t21, 0, t5, -t6, (-t21 * t32 - t23 * t34) * pkin(3), (t32 * t6 + t34 * t5) * pkin(3), 0, 0, t11, 0, -t9, 0, t1, -t2, -t14 * t11 - t15 * t9, t1 * t14 + t2 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t23, 0, (-t21 * t34 + t23 * t32) * pkin(3), 0, 0, 0, 0, 0, 0, -t9, -t11, 0, t11 * t15 - t9 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t45, -0.2e1 * t47, 0, (t32 ^ 2 + t34 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t14, -0.2e1 * t15, 0, t14 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t23, 0, t24, 0, 0, 0, 0, 0, 0, t9, t11, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t14, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t10;
