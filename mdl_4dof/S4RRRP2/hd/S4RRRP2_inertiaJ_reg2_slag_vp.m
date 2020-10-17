% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRP2
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:09
% EndTime: 2019-12-31 17:13:10
% DurationCPUTime: 0.28s
% Computational Cost: add. (106->39), mult. (224->71), div. (0->0), fcn. (178->4), ass. (0->31)
t20 = sin(qJ(2));
t32 = t20 * pkin(1);
t12 = pkin(6) + t32;
t19 = sin(qJ(3));
t17 = t19 ^ 2;
t21 = cos(qJ(3));
t18 = t21 ^ 2;
t26 = t17 + t18;
t34 = t26 * t12;
t37 = 0.2e1 * t19;
t36 = -0.2e1 * t21;
t35 = 0.2e1 * t21;
t33 = t19 * pkin(3);
t22 = cos(qJ(2));
t31 = t22 * pkin(1);
t13 = -pkin(2) - t31;
t30 = pkin(2) - t13;
t14 = -t21 * pkin(3) - pkin(2);
t5 = t14 - t31;
t29 = t14 + t5;
t28 = -qJ(4) - pkin(6);
t27 = t26 * pkin(6);
t25 = qJ(4) + t12;
t10 = t19 * t35;
t7 = t28 * t21;
t6 = t28 * t19;
t4 = t7 * t21;
t3 = t25 * t21;
t2 = t25 * t19;
t1 = t3 * t21;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t31, -0.2e1 * t32, 0, (t20 ^ 2 + t22 ^ 2) * pkin(1) ^ 2, t17, t10, 0, t18, 0, 0, t13 * t36, t13 * t37, 0.2e1 * t34, t26 * t12 ^ 2 + t13 ^ 2, t17, t10, 0, t18, 0, 0, t5 * t36, t5 * t37, 0.2e1 * t2 * t19 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t31, -t32, 0, 0, t17, t10, 0, t18, 0, 0, t30 * t21, -t30 * t19, t27 + t34, -t13 * pkin(2) + pkin(6) * t34, t17, t10, 0, t18, 0, 0, -t29 * t21, t29 * t19, t1 - t4 + (t2 - t6) * t19, t5 * t14 - t2 * t6 - t3 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t17, t10, 0, t18, 0, 0, pkin(2) * t35, -0.2e1 * pkin(2) * t19, 0.2e1 * t27, t26 * pkin(6) ^ 2 + pkin(2) ^ 2, t17, t10, 0, t18, 0, 0, t14 * t36, t14 * t37, -0.2e1 * t6 * t19 - 0.2e1 * t4, t14 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t21, 0, -t19 * t12, -t21 * t12, 0, 0, 0, 0, t19, 0, t21, 0, -t2, -t3, -t33, -t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t21, 0, -t19 * pkin(6), -t21 * pkin(6), 0, 0, 0, 0, t19, 0, t21, 0, t6, t7, -t33, t6 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(3), 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t19, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t19, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
