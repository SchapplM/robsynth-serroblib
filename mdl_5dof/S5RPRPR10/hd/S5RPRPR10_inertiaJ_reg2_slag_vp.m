% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:05
% EndTime: 2019-12-31 18:26:06
% DurationCPUTime: 0.42s
% Computational Cost: add. (307->55), mult. (417->88), div. (0->0), fcn. (405->6), ass. (0->33)
t28 = sin(pkin(8));
t29 = cos(pkin(8));
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t12 = t28 * t31 - t29 * t33;
t43 = t12 ^ 2;
t42 = -pkin(1) - pkin(2);
t41 = t29 * pkin(3);
t16 = t31 * qJ(2) - t33 * t42;
t15 = -pkin(3) - t16;
t18 = t33 * qJ(2) + t31 * t42;
t7 = t28 * t15 + t29 * t18;
t23 = -pkin(4) - t41;
t5 = -t29 * t15 + t28 * t18;
t3 = pkin(4) + t5;
t40 = t23 - t3;
t30 = sin(qJ(5));
t39 = t12 * t30;
t32 = cos(qJ(5));
t38 = t12 * t32;
t37 = t30 * t32;
t26 = t30 ^ 2;
t27 = t32 ^ 2;
t36 = t26 + t27;
t14 = t28 * t33 + t29 * t31;
t1 = t36 * t14;
t25 = t28 * pkin(3);
t22 = t25 + pkin(7);
t35 = t36 * t22;
t20 = 0.2e1 * t37;
t11 = t14 ^ 2;
t4 = -pkin(7) + t7;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2 * pkin(1), 0, 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t16, 0.2e1 * t18, 0, t16 ^ 2 + t18 ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t5, 0.2e1 * t7, 0, t5 ^ 2 + t7 ^ 2, t26, t20, 0, t27, 0, 0, 0.2e1 * t3 * t32, -0.2e1 * t3 * t30, -0.2e1 * t36 * t4, t36 * t4 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t33, t31, 0, -t16 * t33 + t18 * t31, 0, 0, 0, 0, 0, 0, t12, t14, 0, t5 * t12 + t7 * t14, 0, 0, 0, 0, 0, 0, t38, -t39, -t1, t4 * t1 + t3 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 ^ 2 + t33 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t11 + t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t16, -t18, 0, 0, 0, 0, 0, 0, 0, -1, -t5 - t41, t25 - t7, 0, (t28 * t7 - t29 * t5) * pkin(3), -t26, -0.2e1 * t37, 0, -t27, 0, 0, t40 * t32, -t40 * t30, (-t22 + t4) * t36, t3 * t23 + t4 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t14, 0, (-t12 * t29 + t14 * t28) * pkin(3), 0, 0, 0, 0, 0, 0, -t38, t39, t1, t12 * t23 + t14 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t25, 0, (t28 ^ 2 + t29 ^ 2) * pkin(3) ^ 2, t26, t20, 0, t27, 0, 0, -0.2e1 * t23 * t32, 0.2e1 * t23 * t30, 0.2e1 * t35, t36 * t22 ^ 2 + t23 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, -t32, 0, -t30 * t4, -t32 * t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t14, -t32 * t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t32, 0, -t30 * t22, -t32 * t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t2;
