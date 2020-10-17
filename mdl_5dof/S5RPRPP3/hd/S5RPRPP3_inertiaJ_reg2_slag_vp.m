% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:54
% EndTime: 2019-12-31 18:12:55
% DurationCPUTime: 0.40s
% Computational Cost: add. (216->40), mult. (421->71), div. (0->0), fcn. (466->4), ass. (0->34)
t21 = sin(pkin(7));
t22 = cos(pkin(7));
t25 = sin(qJ(3));
t38 = cos(qJ(3));
t12 = t25 * t21 - t38 * t22;
t9 = t12 ^ 2;
t14 = t38 * t21 + t25 * t22;
t10 = t14 ^ 2;
t41 = -0.2e1 * t14;
t18 = -t22 * pkin(2) - pkin(1);
t40 = 0.2e1 * t18;
t39 = 0.2e1 * t22;
t37 = t12 * t14;
t23 = pkin(3) + qJ(5);
t36 = pkin(6) + qJ(2);
t19 = t21 ^ 2;
t20 = t22 ^ 2;
t35 = t19 + t20;
t34 = qJ(4) * t12;
t33 = -0.2e1 * t37;
t16 = t36 * t22;
t31 = t36 * t21;
t5 = t25 * t16 + t38 * t31;
t7 = t38 * t16 - t25 * t31;
t32 = t5 ^ 2 + t7 ^ 2;
t30 = -t14 * qJ(4) + t18;
t29 = -0.2e1 * t7 * t12 + 0.2e1 * t5 * t14;
t27 = qJ(4) ^ 2;
t26 = 0.2e1 * qJ(4);
t4 = t12 * pkin(3) + t30;
t3 = -t12 * pkin(4) + t7;
t2 = t14 * pkin(4) + t5;
t1 = t23 * t12 + t30;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t19, t21 * t39, 0, t20, 0, 0, pkin(1) * t39, -0.2e1 * pkin(1) * t21, 0.2e1 * t35 * qJ(2), t35 * qJ(2) ^ 2 + pkin(1) ^ 2, t10, t33, 0, t9, 0, 0, t12 * t40, t14 * t40, t29, t18 ^ 2 + t32, 0, 0, 0, t10, t33, t9, t29, -0.2e1 * t4 * t12, t4 * t41, t4 ^ 2 + t32, 0, 0, 0, t9, 0.2e1 * t37, t10, -0.2e1 * t3 * t12 + 0.2e1 * t2 * t14, t1 * t41, 0.2e1 * t1 * t12, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t12, t14, 0, t18, 0, 0, 0, 0, 0, 0, 0, -t12, -t14, t4, 0, 0, 0, 0, 0, 0, 0, -t14, t12, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, -t12, 0, -t5, -t7, 0, 0, 0, -t14, t12, 0, 0, 0, -pkin(3) * t14 - t34, t5, t7, -t5 * pkin(3) + t7 * qJ(4), 0, t12, t14, 0, 0, 0, -t23 * t14 - t34, t3, -t2, t3 * qJ(4) - t2 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(3), t26, pkin(3) ^ 2 + t27, 1, 0, 0, 0, 0, 0, 0, t26, 0.2e1 * t23, t23 ^ 2 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, t5, 0, 0, 0, 0, 0, 0, t14, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -1, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
