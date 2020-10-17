% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPR3
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:35
% DurationCPUTime: 0.22s
% Computational Cost: add. (106->27), mult. (217->55), div. (0->0), fcn. (172->6), ass. (0->28)
t19 = sin(pkin(7));
t34 = t19 * pkin(2);
t13 = pkin(6) + t34;
t21 = sin(qJ(4));
t17 = t21 ^ 2;
t23 = cos(qJ(4));
t18 = t23 ^ 2;
t29 = t17 + t18;
t31 = t29 * t13;
t37 = 0.2e1 * t21;
t36 = -0.2e1 * t23;
t24 = cos(qJ(2));
t16 = t24 * pkin(1);
t15 = t16 + pkin(2);
t20 = cos(pkin(7));
t22 = sin(qJ(2));
t32 = t22 * pkin(1);
t28 = t20 * t32;
t7 = t19 * t15 + t28;
t5 = pkin(6) + t7;
t35 = t29 * t5;
t33 = t20 * pkin(2);
t14 = -pkin(3) - t33;
t27 = -t20 * t15 + t19 * t32;
t4 = -pkin(3) + t27;
t30 = t14 + t4;
t12 = t23 * t37;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t16, -0.2e1 * t32, 0, (t22 ^ 2 + t24 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t27, -0.2e1 * t7, 0, t27 ^ 2 + t7 ^ 2, t17, t12, 0, t18, 0, 0, t4 * t36, t4 * t37, 0.2e1 * t35, t29 * t5 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t16, -t32, 0, 0, 0, 0, 0, 0, 0, 1, -t27 + t33, -t28 + (-pkin(2) - t15) * t19, 0, (t19 * t7 - t20 * t27) * pkin(2), t17, t12, 0, t18, 0, 0, -t30 * t23, t30 * t21, t31 + t35, t4 * t14 + t5 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t33, -0.2e1 * t34, 0, (t19 ^ 2 + t20 ^ 2) * pkin(2) ^ 2, t17, t12, 0, t18, 0, 0, t14 * t36, t14 * t37, 0.2e1 * t31, t29 * t13 ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, t23, 0, -t21 * t5, -t23 * t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, t23, 0, -t21 * t13, -t23 * t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t21, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
