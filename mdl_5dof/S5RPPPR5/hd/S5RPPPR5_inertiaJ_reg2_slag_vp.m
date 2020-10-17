% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPPR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:30
% EndTime: 2019-12-31 17:46:31
% DurationCPUTime: 0.36s
% Computational Cost: add. (226->46), mult. (335->74), div. (0->0), fcn. (365->6), ass. (0->32)
t27 = sin(pkin(7));
t26 = sin(pkin(8));
t22 = t26 ^ 2;
t28 = cos(pkin(8));
t24 = t28 ^ 2;
t33 = t22 + t24;
t40 = t33 * t27;
t30 = sin(qJ(5));
t34 = cos(qJ(5));
t32 = -t30 * t26 + t34 * t28;
t39 = t32 ^ 2;
t11 = t34 * t26 + t30 * t28;
t38 = t11 ^ 2;
t37 = -0.2e1 * t11;
t36 = 0.2e1 * t28;
t29 = cos(pkin(7));
t31 = -pkin(1) - pkin(2);
t17 = t29 * qJ(2) + t27 * t31;
t13 = -qJ(4) + t17;
t35 = pkin(6) - t13;
t15 = t27 * qJ(2) - t29 * t31;
t14 = pkin(3) + t15;
t25 = t29 ^ 2;
t23 = t27 ^ 2;
t7 = t28 * pkin(4) + t14;
t6 = t32 * t27;
t5 = t11 * t27;
t4 = t35 * t28;
t3 = t35 * t26;
t2 = t30 * t3 - t34 * t4;
t1 = t34 * t3 + t30 * t4;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2 * pkin(1), 0, 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t15, 0.2e1 * t17, 0, t15 ^ 2 + t17 ^ 2, t22, t26 * t36, 0, t24, 0, 0, t14 * t36, -0.2e1 * t14 * t26, -0.2e1 * t33 * t13, t33 * t13 ^ 2 + t14 ^ 2, t38, -t32 * t37, 0, t39, 0, 0, 0.2e1 * t7 * t32, t7 * t37, 0.2e1 * t1 * t11 - 0.2e1 * t2 * t32, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t29, t27, 0, -t15 * t29 + t17 * t27, 0, 0, 0, 0, 0, 0, -t29 * t28, t29 * t26, -t40, t13 * t40 - t14 * t29, 0, 0, 0, 0, 0, 0, -t29 * t32, t29 * t11, -t5 * t11 - t32 * t6, -t1 * t5 + t2 * t6 - t7 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 + t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 * t23 + t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t32 + t2 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t11 - t32 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26, 0, t14, 0, 0, 0, 0, 0, 0, t32, -t11, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t32, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
