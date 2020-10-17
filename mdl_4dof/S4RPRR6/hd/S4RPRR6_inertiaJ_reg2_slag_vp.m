% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPRR6
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:41
% EndTime: 2019-12-31 16:52:42
% DurationCPUTime: 0.23s
% Computational Cost: add. (219->33), mult. (447->76), div. (0->0), fcn. (510->6), ass. (0->30)
t22 = sin(pkin(7));
t23 = cos(pkin(7));
t25 = sin(qJ(3));
t27 = cos(qJ(3));
t13 = t25 * t22 - t27 * t23;
t19 = -t23 * pkin(2) - pkin(1);
t10 = t13 * pkin(3) + t19;
t36 = 0.2e1 * t10;
t35 = 0.2e1 * t19;
t34 = 0.2e1 * t23;
t24 = sin(qJ(4));
t33 = t24 * pkin(3);
t26 = cos(qJ(4));
t32 = t26 * pkin(3);
t31 = pkin(5) + qJ(2);
t20 = t22 ^ 2;
t21 = t23 ^ 2;
t30 = t20 + t21;
t16 = t31 * t22;
t17 = t31 * t23;
t8 = -t27 * t16 - t25 * t17;
t9 = -t25 * t16 + t27 * t17;
t15 = t27 * t22 + t25 * t23;
t7 = -t24 * t13 + t26 * t15;
t5 = t26 * t13 + t24 * t15;
t4 = -t13 * pkin(6) + t9;
t3 = -t15 * pkin(6) + t8;
t2 = t24 * t3 + t26 * t4;
t1 = -t24 * t4 + t26 * t3;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t20, t22 * t34, 0, t21, 0, 0, pkin(1) * t34, -0.2e1 * pkin(1) * t22, 0.2e1 * t30 * qJ(2), t30 * qJ(2) ^ 2 + pkin(1) ^ 2, t15 ^ 2, -0.2e1 * t15 * t13, 0, t13 ^ 2, 0, 0, t13 * t35, t15 * t35, -0.2e1 * t9 * t13 - 0.2e1 * t8 * t15, t19 ^ 2 + t8 ^ 2 + t9 ^ 2, t7 ^ 2, -0.2e1 * t7 * t5, 0, t5 ^ 2, 0, 0, t5 * t36, t7 * t36, -0.2e1 * t1 * t7 - 0.2e1 * t2 * t5, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t22, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t13, t15, 0, t19, 0, 0, 0, 0, 0, 0, t5, t7, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t13, 0, t8, -t9, 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, (-t24 * t5 - t26 * t7) * pkin(3), (t1 * t26 + t2 * t24) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t32, -0.2e1 * t33, 0, (t24 ^ 2 + t26 ^ 2) * pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t32, -t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
