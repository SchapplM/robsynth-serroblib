% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPR5
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:33
% EndTime: 2019-12-31 17:03:34
% DurationCPUTime: 0.20s
% Computational Cost: add. (80->25), mult. (136->41), div. (0->0), fcn. (87->4), ass. (0->24)
t13 = sin(qJ(4));
t10 = t13 ^ 2;
t15 = cos(qJ(4));
t11 = t15 ^ 2;
t3 = t10 + t11;
t14 = sin(qJ(2));
t26 = t14 * pkin(1);
t7 = qJ(3) + t26;
t28 = t7 ^ 2;
t27 = 0.2e1 * t7;
t18 = 0.2e1 * qJ(3);
t16 = cos(qJ(2));
t25 = t16 * pkin(1);
t23 = t7 * qJ(3);
t22 = qJ(3) + t7;
t9 = -pkin(2) - t25;
t17 = -pkin(2) - pkin(6);
t2 = t3 * t17;
t20 = -0.2e1 * pkin(2);
t19 = qJ(3) ^ 2;
t6 = -0.2e1 * t15 * t13;
t5 = -pkin(6) + t9;
t1 = t3 * t5;
t4 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t26, 0, (t14 ^ 2 + t16 ^ 2) * pkin(1) ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t9, t27, t9 ^ 2 + t28, t11, t6, 0, t10, 0, 0, t13 * t27, t15 * t27, -0.2e1 * t1, t3 * t5 ^ 2 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, -t26, 0, 0, 1, 0, 0, 0, 0, 0, 0, t20 - t25, t18 + t26, -t9 * pkin(2) + t23, t11, t6, 0, t10, 0, 0, t22 * t13, t22 * t15, t3 * (-t17 - t5), t5 * t2 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, t20, t18, pkin(2) ^ 2 + t19, t11, t6, 0, t10, 0, 0, t13 * t18, t15 * t18, -0.2e1 * t2, t3 * t17 ^ 2 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t13, 0, t15 * t5, -t13 * t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t13, 0, t15 * t17, -t13 * t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t4;
