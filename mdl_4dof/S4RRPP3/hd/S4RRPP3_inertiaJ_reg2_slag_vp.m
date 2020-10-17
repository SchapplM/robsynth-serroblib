% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:51
% EndTime: 2019-12-31 16:57:52
% DurationCPUTime: 0.24s
% Computational Cost: add. (130->30), mult. (282->69), div. (0->0), fcn. (289->4), ass. (0->28)
t20 = sin(pkin(6));
t21 = cos(pkin(6));
t22 = sin(qJ(2));
t23 = cos(qJ(2));
t8 = t20 * t22 - t21 * t23;
t36 = t8 ^ 2;
t17 = -t23 * pkin(2) - pkin(1);
t35 = 0.2e1 * t17;
t34 = 0.2e1 * t23;
t10 = t20 * t23 + t21 * t22;
t33 = t10 * t8;
t32 = t20 * pkin(2);
t31 = t21 * pkin(2);
t30 = -qJ(3) - pkin(5);
t18 = t22 ^ 2;
t19 = t23 ^ 2;
t29 = t18 + t19;
t12 = t30 * t23;
t27 = t30 * t22;
t4 = -t20 * t12 - t21 * t27;
t6 = -t21 * t12 + t20 * t27;
t28 = t4 ^ 2 + t6 ^ 2;
t26 = 0.2e1 * t4 * t10 - 0.2e1 * t6 * t8;
t15 = pkin(3) + t31;
t13 = qJ(4) + t32;
t7 = t10 ^ 2;
t2 = t8 * pkin(3) - t10 * qJ(4) + t17;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t18, t22 * t34, 0, t19, 0, 0, pkin(1) * t34, -0.2e1 * pkin(1) * t22, 0.2e1 * t29 * pkin(5), t29 * pkin(5) ^ 2 + pkin(1) ^ 2, t7, -0.2e1 * t33, 0, t36, 0, 0, t8 * t35, t10 * t35, t26, t17 ^ 2 + t28, t7, 0, 0.2e1 * t33, 0, 0, t36, 0.2e1 * t2 * t8, t26, -0.2e1 * t2 * t10, t2 ^ 2 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t23, 0, -t22 * pkin(5), -t23 * pkin(5), 0, 0, 0, 0, t10, 0, -t8, 0, -t4, -t6, (-t10 * t21 - t20 * t8) * pkin(2), (t20 * t6 - t21 * t4) * pkin(2), 0, t10, 0, 0, t8, 0, -t4, -t15 * t10 - t13 * t8, t6, t6 * t13 - t4 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t31, -0.2e1 * t32, 0, (t20 ^ 2 + t21 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t15, 0, 0.2e1 * t13, t13 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10, 0, t17, 0, 0, 0, 0, 0, 0, t8, 0, -t10, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
