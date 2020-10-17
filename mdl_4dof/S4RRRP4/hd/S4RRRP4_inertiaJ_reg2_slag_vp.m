% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:41
% EndTime: 2019-12-31 17:15:42
% DurationCPUTime: 0.25s
% Computational Cost: add. (156->38), mult. (340->74), div. (0->0), fcn. (352->4), ass. (0->33)
t24 = sin(qJ(3));
t25 = sin(qJ(2));
t26 = cos(qJ(3));
t27 = cos(qJ(2));
t10 = t24 * t25 - t26 * t27;
t18 = -t27 * pkin(2) - pkin(1);
t6 = t10 * pkin(3) + t18;
t36 = 0.2e1 * t6;
t35 = 0.2e1 * t18;
t34 = 0.2e1 * t27;
t33 = -pkin(6) - pkin(5);
t32 = t24 * pkin(2);
t21 = t26 * pkin(2);
t22 = t25 ^ 2;
t23 = t27 ^ 2;
t31 = t22 + t23;
t14 = t33 * t25;
t15 = t33 * t27;
t4 = t26 * t14 + t24 * t15;
t5 = t24 * t14 - t26 * t15;
t30 = pkin(2) ^ 2;
t28 = 0.2e1 * pkin(3);
t20 = t24 ^ 2 * t30;
t19 = -0.2e1 * t32;
t17 = t21 + pkin(3);
t12 = t24 * t27 + t26 * t25;
t9 = t12 ^ 2;
t8 = t10 ^ 2;
t7 = t10 * t32;
t3 = -0.2e1 * t12 * t10;
t2 = -t10 * qJ(4) + t5;
t1 = -t12 * qJ(4) + t4;
t11 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t22, t25 * t34, 0, t23, 0, 0, pkin(1) * t34, -0.2e1 * pkin(1) * t25, 0.2e1 * t31 * pkin(5), t31 * pkin(5) ^ 2 + pkin(1) ^ 2, t9, t3, 0, t8, 0, 0, t10 * t35, t12 * t35, -0.2e1 * t5 * t10 - 0.2e1 * t4 * t12, t18 ^ 2 + t4 ^ 2 + t5 ^ 2, t9, t3, 0, t8, 0, 0, t10 * t36, t12 * t36, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t27, 0, -t25 * pkin(5), -t27 * pkin(5), 0, 0, 0, 0, t12, 0, -t10, 0, t4, -t5, -t12 * t21 - t7, (t24 * t5 + t26 * t4) * pkin(2), 0, 0, t12, 0, -t10, 0, t1, -t2, -t17 * t12 - t7, t1 * t17 + t2 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t21, t19, 0, t26 ^ 2 * t30 + t20, 0, 0, 0, 0, 0, 1, 0.2e1 * t17, t19, 0, t17 ^ 2 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, 0, t4, -t5, 0, 0, 0, 0, t12, 0, -t10, 0, t1, -t2, -t12 * pkin(3), t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t21, -t32, 0, 0, 0, 0, 0, 0, 0, 1, t28 + t21, -t32, 0, t17 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t28, 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t11;
