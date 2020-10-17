% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPPP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:07:31
% EndTime: 2019-05-04 19:07:32
% DurationCPUTime: 0.22s
% Computational Cost: add. (89->34), mult. (289->79), div. (0->0), fcn. (262->4), ass. (0->33)
t22 = sin(pkin(4));
t34 = 0.2e1 * t22;
t21 = sin(pkin(6));
t33 = pkin(1) * t21;
t19 = t22 ^ 2;
t23 = cos(pkin(6));
t32 = t19 * t23;
t17 = t22 * t21;
t18 = t22 * t23;
t24 = cos(pkin(4));
t31 = t23 * t24;
t30 = qJ(2) * t22;
t8 = t23 * t30 + t24 * t33;
t29 = t21 * t32;
t28 = t24 * t17;
t27 = t22 * t31;
t26 = -pkin(1) * t23 - pkin(2);
t25 = -qJ(3) * t21 - pkin(1);
t4 = -t24 * qJ(3) - t8;
t20 = t24 ^ 2;
t16 = t19 * t23 ^ 2;
t15 = t19 * t21 ^ 2;
t12 = t21 * t30;
t11 = -0.2e1 * t27;
t10 = 0.2e1 * t28;
t9 = 0.2e1 * t29;
t7 = pkin(1) * t31 - t12;
t6 = (-pkin(2) * t23 + t25) * t22;
t5 = t26 * t24 + t12;
t3 = ((-pkin(2) - qJ(4)) * t23 + t25) * t22;
t2 = pkin(3) * t18 - t4;
t1 = pkin(3) * t17 + t12 + (-qJ(4) + t26) * t24;
t13 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t15, t9, t10, t16, 0.2e1 * t27, t20, 0.2e1 * pkin(1) * t32 + 0.2e1 * t7 * t24, -0.2e1 * t19 * t33 - 0.2e1 * t8 * t24 (-t21 * t7 + t23 * t8) * t34, t19 * pkin(1) ^ 2 + t7 ^ 2 + t8 ^ 2, t20, -0.2e1 * t28, t11, t15, t9, t16 (t21 * t5 - t23 * t4) * t34, 0.2e1 * t6 * t18 + 0.2e1 * t5 * t24, -0.2e1 * t6 * t17 - 0.2e1 * t4 * t24, t4 ^ 2 + t5 ^ 2 + t6 ^ 2, t20, t11, t10, t16, -0.2e1 * t29, t15 (t1 * t21 + t2 * t23) * t34, -0.2e1 * t3 * t17 + 0.2e1 * t2 * t24, -0.2e1 * t1 * t24 - 0.2e1 * t3 * t18, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0, -t22 * pkin(1), 0, 0, 0, 0, 0, 0, 0, t18, -t17, t6, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t24, 0, t5, 0, 0, 0, 0, 0, 0, t17, 0, -t24, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t24, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t13;
