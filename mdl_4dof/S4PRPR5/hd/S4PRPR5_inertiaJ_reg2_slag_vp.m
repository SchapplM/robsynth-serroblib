% Calculate inertial parameters regressor of joint inertia matrix for
% S4PRPR5
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRPR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:17
% EndTime: 2019-12-31 16:23:18
% DurationCPUTime: 0.15s
% Computational Cost: add. (50->21), mult. (113->41), div. (0->0), fcn. (122->6), ass. (0->20)
t11 = sin(pkin(7));
t12 = cos(pkin(7));
t14 = sin(qJ(2));
t16 = cos(qJ(2));
t2 = t11 * t14 - t12 * t16;
t23 = t2 ^ 2;
t13 = sin(qJ(4));
t22 = 0.2e1 * t13;
t21 = t11 * pkin(2);
t20 = t12 * pkin(2);
t15 = cos(qJ(4));
t10 = t15 ^ 2;
t9 = t13 ^ 2;
t19 = t9 + t10;
t7 = pkin(5) + t21;
t18 = t19 * t7;
t8 = -pkin(3) - t20;
t4 = t11 * t16 + t12 * t14;
t1 = t4 ^ 2;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t16 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 + t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t1 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t4, 0, (t11 * t4 - t12 * t2) * pkin(2), 0, 0, 0, 0, 0, 0, -t2 * t15, t2 * t13, t19 * t4, t4 * t18 + t2 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t20, -0.2e1 * t21, 0, (t11 ^ 2 + t12 ^ 2) * pkin(2) ^ 2, t9, t15 * t22, 0, t10, 0, 0, -0.2e1 * t8 * t15, t8 * t22, 0.2e1 * t18, t19 * t7 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t4, -t15 * t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t15, 0, -t13 * t7, -t15 * t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
