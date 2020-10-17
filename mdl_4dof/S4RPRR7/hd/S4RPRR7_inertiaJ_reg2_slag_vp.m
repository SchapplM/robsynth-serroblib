% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPRR7
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:05
% EndTime: 2019-12-31 16:54:06
% DurationCPUTime: 0.28s
% Computational Cost: add. (241->41), mult. (506->95), div. (0->0), fcn. (552->6), ass. (0->40)
t20 = sin(pkin(7));
t21 = cos(pkin(7));
t23 = sin(qJ(3));
t40 = cos(qJ(3));
t12 = t40 * t20 + t23 * t21;
t45 = -0.2e1 * t12;
t35 = pkin(5) + qJ(2);
t13 = t35 * t21;
t30 = t35 * t20;
t4 = t23 * t13 + t40 * t30;
t44 = t4 ^ 2;
t10 = t23 * t20 - t40 * t21;
t43 = t10 ^ 2;
t15 = -t21 * pkin(2) - pkin(1);
t42 = 0.2e1 * t15;
t41 = 0.2e1 * t21;
t22 = sin(qJ(4));
t39 = t22 * t10;
t38 = t22 * t12;
t24 = cos(qJ(4));
t37 = t22 * t24;
t36 = t24 * t12;
t16 = t20 ^ 2;
t17 = t21 ^ 2;
t34 = t16 + t17;
t18 = t22 ^ 2;
t19 = t24 ^ 2;
t33 = t18 + t19;
t32 = t10 * t45;
t31 = t22 * t36;
t29 = -pkin(3) * t12 - pkin(6) * t10;
t3 = t10 * pkin(3) - t12 * pkin(6) + t15;
t6 = t40 * t13 - t23 * t30;
t1 = -t22 * t6 + t24 * t3;
t2 = t22 * t3 + t24 * t6;
t28 = t1 * t24 + t2 * t22;
t27 = -t1 * t22 + t2 * t24;
t8 = t12 ^ 2;
t7 = t24 * t10;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t16, t20 * t41, 0, t17, 0, 0, pkin(1) * t41, -0.2e1 * pkin(1) * t20, 0.2e1 * t34 * qJ(2), t34 * qJ(2) ^ 2 + pkin(1) ^ 2, t8, t32, 0, t43, 0, 0, t10 * t42, t12 * t42, -0.2e1 * t6 * t10 + 0.2e1 * t4 * t12, t15 ^ 2 + t6 ^ 2 + t44, t19 * t8, -0.2e1 * t8 * t37, 0.2e1 * t10 * t36, t18 * t8, t22 * t32, t43, 0.2e1 * t1 * t10 + 0.2e1 * t4 * t38, -0.2e1 * t2 * t10 + 0.2e1 * t4 * t36, t28 * t45, t1 ^ 2 + t2 ^ 2 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t20, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t10, t12, 0, t15, 0, 0, 0, 0, 0, 0, t7, -t39, -t33 * t12, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, 0, -t4, -t6, 0, 0, t31, (-t18 + t19) * t12, t39, -t31, t7, 0, t29 * t22 - t4 * t24, t4 * t22 + t29 * t24, t27, -t4 * pkin(3) + t27 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t18, 0.2e1 * t37, 0, t19, 0, 0, 0.2e1 * pkin(3) * t24, -0.2e1 * pkin(3) * t22, 0.2e1 * t33 * pkin(6), t33 * pkin(6) ^ 2 + pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, -t38, t10, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t24, 0, -t22 * pkin(6), -t24 * pkin(6), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
