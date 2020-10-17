% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPR6
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:54
% EndTime: 2019-12-31 17:04:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (249->39), mult. (513->90), div. (0->0), fcn. (570->6), ass. (0->33)
t25 = sin(pkin(7));
t26 = cos(pkin(7));
t28 = sin(qJ(2));
t29 = cos(qJ(2));
t14 = t25 * t28 - t26 * t29;
t22 = -t29 * pkin(2) - pkin(1);
t10 = t14 * pkin(3) + t22;
t39 = 0.2e1 * t10;
t38 = 0.2e1 * t22;
t37 = 0.2e1 * t29;
t36 = t25 * pkin(2);
t35 = t26 * pkin(2);
t34 = cos(qJ(4));
t33 = -qJ(3) - pkin(5);
t23 = t28 ^ 2;
t24 = t29 ^ 2;
t32 = t23 + t24;
t18 = t33 * t28;
t19 = t33 * t29;
t8 = t26 * t18 + t25 * t19;
t9 = t25 * t18 - t26 * t19;
t27 = sin(qJ(4));
t21 = pkin(3) + t35;
t16 = t25 * t29 + t26 * t28;
t12 = t27 * t21 + t34 * t36;
t11 = t34 * t21 - t27 * t36;
t7 = -t27 * t14 + t34 * t16;
t5 = t34 * t14 + t27 * t16;
t4 = -t14 * pkin(6) + t9;
t3 = -t16 * pkin(6) + t8;
t2 = t27 * t3 + t34 * t4;
t1 = -t27 * t4 + t34 * t3;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t23, t28 * t37, 0, t24, 0, 0, pkin(1) * t37, -0.2e1 * pkin(1) * t28, 0.2e1 * t32 * pkin(5), t32 * pkin(5) ^ 2 + pkin(1) ^ 2, t16 ^ 2, -0.2e1 * t16 * t14, 0, t14 ^ 2, 0, 0, t14 * t38, t16 * t38, -0.2e1 * t9 * t14 - 0.2e1 * t8 * t16, t22 ^ 2 + t8 ^ 2 + t9 ^ 2, t7 ^ 2, -0.2e1 * t7 * t5, 0, t5 ^ 2, 0, 0, t5 * t39, t7 * t39, -0.2e1 * t1 * t7 - 0.2e1 * t2 * t5, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t29, 0, -t28 * pkin(5), -t29 * pkin(5), 0, 0, 0, 0, t16, 0, -t14, 0, t8, -t9, (-t14 * t25 - t16 * t26) * pkin(2), (t25 * t9 + t26 * t8) * pkin(2), 0, 0, t7, 0, -t5, 0, t1, -t2, -t11 * t7 - t12 * t5, t1 * t11 + t2 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t36, 0, (t25 ^ 2 + t26 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t11, -0.2e1 * t12, 0, t11 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, t22, 0, 0, 0, 0, 0, 0, t5, t7, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t11, -t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
