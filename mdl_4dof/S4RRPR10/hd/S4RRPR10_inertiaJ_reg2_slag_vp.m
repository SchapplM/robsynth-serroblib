% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:55
% EndTime: 2019-12-31 17:11:56
% DurationCPUTime: 0.28s
% Computational Cost: add. (136->42), mult. (273->80), div. (0->0), fcn. (235->4), ass. (0->41)
t23 = sin(qJ(2));
t18 = t23 ^ 2;
t25 = cos(qJ(2));
t20 = t25 ^ 2;
t44 = t18 + t20;
t43 = -0.2e1 * t23;
t42 = 0.2e1 * t25;
t41 = 2 * qJ(3);
t26 = -pkin(2) - pkin(6);
t22 = sin(qJ(4));
t40 = t22 * t23;
t39 = t22 * t25;
t38 = t23 * t25;
t24 = cos(qJ(4));
t37 = t24 * t22;
t36 = t24 * t25;
t35 = t44 * pkin(5) ^ 2;
t17 = t22 ^ 2;
t19 = t24 ^ 2;
t10 = t17 + t19;
t34 = t25 * qJ(3);
t33 = -0.2e1 * t38;
t32 = t22 * t36;
t31 = -t23 * qJ(3) - pkin(1);
t4 = t26 * t25 + t31;
t14 = t23 * pkin(5);
t8 = t23 * pkin(3) + t14;
t2 = -t22 * t4 + t24 * t8;
t3 = t22 * t8 + t24 * t4;
t1 = t2 * t24 + t3 * t22;
t30 = -t23 * pkin(2) + t34;
t29 = t23 * t26 + t34;
t27 = qJ(3) ^ 2;
t16 = t25 * pkin(5);
t12 = t24 * t23;
t11 = 0.2e1 * t38;
t9 = t25 * pkin(3) + t16;
t7 = -t25 * pkin(2) + t31;
t6 = 0.2e1 * t44 * pkin(5);
t5 = t10 * t26;
t13 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t18, t11, 0, t20, 0, 0, pkin(1) * t42, pkin(1) * t43, t6, pkin(1) ^ 2 + t35, 0, 0, 0, t18, t11, t20, t6, t7 * t42, t7 * t43, t7 ^ 2 + t35, t17 * t20, 0.2e1 * t20 * t37, t22 * t33, t19 * t20, t24 * t33, t18, 0.2e1 * t2 * t23 + 0.2e1 * t9 * t36, -0.2e1 * t3 * t23 - 0.2e1 * t9 * t39, (t2 * t22 - t24 * t3) * t42, t2 ^ 2 + t3 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t25, 0, -t14, -t16, 0, 0, 0, -t23, -t25, 0, 0, 0, t30, t14, t16, t30 * pkin(5), -t32, (t17 - t19) * t25, t12, t32, -t40, 0, t9 * t22 + t29 * t24, -t29 * t22 + t9 * t24, -t1, t9 * qJ(3) + t1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(2), t41, pkin(2) ^ 2 + t27, t19, -0.2e1 * t37, 0, t17, 0, 0, t22 * t41, t24 * t41, -0.2e1 * t5, t10 * t26 ^ 2 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, t14, 0, 0, 0, 0, 0, 0, t12, -t40, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t10, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, -t36, t23, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, -t22, 0, t24 * t26, -t22 * t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t13;
