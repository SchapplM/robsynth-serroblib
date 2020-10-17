% Calculate inertial parameters regressor of joint inertia matrix for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRPR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:41
% EndTime: 2019-12-31 16:24:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (78->30), mult. (185->56), div. (0->0), fcn. (205->6), ass. (0->25)
t17 = cos(pkin(7));
t11 = -t17 * pkin(3) - pkin(2);
t27 = 0.2e1 * t11;
t26 = 0.2e1 * t17;
t25 = pkin(5) + qJ(3);
t16 = sin(pkin(7));
t12 = t16 ^ 2;
t13 = t17 ^ 2;
t24 = t12 + t13;
t23 = t24 * qJ(3);
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t7 = t20 * t16 + t18 * t17;
t5 = t18 * t16 - t20 * t17;
t21 = cos(qJ(2));
t19 = sin(qJ(2));
t15 = t21 ^ 2;
t14 = t19 ^ 2;
t9 = t25 * t17;
t8 = t25 * t16;
t4 = t5 * t19;
t3 = t7 * t19;
t2 = -t18 * t8 + t20 * t9;
t1 = -t18 * t9 - t20 * t8;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 + t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t14 + t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t17, -t21 * t16, t24 * t19, t21 * pkin(2) + t19 * t23, 0, 0, 0, 0, 0, 0, -t21 * t5, -t21 * t7, t3 * t7 + t4 * t5, -t3 * t1 - t21 * t11 - t4 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t12, t16 * t26, 0, t13, 0, 0, pkin(2) * t26, -0.2e1 * pkin(2) * t16, 0.2e1 * t23, t24 * qJ(3) ^ 2 + pkin(2) ^ 2, t7 ^ 2, -0.2e1 * t7 * t5, 0, t5 ^ 2, 0, 0, t5 * t27, t7 * t27, -0.2e1 * t1 * t7 - 0.2e1 * t2 * t5, t1 ^ 2 + t11 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, -pkin(2), 0, 0, 0, 0, 0, 0, t5, t7, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
