% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPRR3
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:18
% EndTime: 2019-12-31 16:49:19
% DurationCPUTime: 0.19s
% Computational Cost: add. (122->32), mult. (247->66), div. (0->0), fcn. (257->6), ass. (0->28)
t18 = cos(pkin(7));
t29 = t18 * pkin(1);
t13 = -pkin(2) - t29;
t22 = cos(qJ(3));
t10 = -t22 * pkin(3) + t13;
t32 = 0.2e1 * t10;
t20 = sin(qJ(3));
t31 = 0.2e1 * t20;
t17 = sin(pkin(7));
t30 = t17 * pkin(1);
t19 = sin(qJ(4));
t28 = t19 * pkin(3);
t21 = cos(qJ(4));
t27 = t21 * pkin(3);
t12 = pkin(5) + t30;
t26 = pkin(6) + t12;
t15 = t20 ^ 2;
t16 = t22 ^ 2;
t25 = t15 + t16;
t9 = t19 * t22 + t21 * t20;
t7 = t19 * t20 - t21 * t22;
t6 = t9 ^ 2;
t5 = t7 ^ 2;
t4 = t26 * t22;
t3 = t26 * t20;
t2 = -t19 * t3 + t21 * t4;
t1 = -t19 * t4 - t21 * t3;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t30, 0, (t17 ^ 2 + t18 ^ 2) * pkin(1) ^ 2, t15, t22 * t31, 0, t16, 0, 0, -0.2e1 * t13 * t22, t13 * t31, 0.2e1 * t25 * t12, t25 * t12 ^ 2 + t13 ^ 2, t6, -0.2e1 * t9 * t7, 0, t5, 0, 0, t7 * t32, t9 * t32, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t7 + t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, t22, 0, -t20 * t12, -t22 * t12, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, (-t19 * t7 - t21 * t9) * pkin(3), (t1 * t21 + t19 * t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, 0, (t19 * t9 - t21 * t7) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t27, -0.2e1 * t28, 0, (t19 ^ 2 + t21 ^ 2) * pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
