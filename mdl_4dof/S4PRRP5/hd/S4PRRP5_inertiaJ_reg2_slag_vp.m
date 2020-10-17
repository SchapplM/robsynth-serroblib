% Calculate inertial parameters regressor of joint inertia matrix for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:22
% EndTime: 2019-12-31 16:29:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (44->24), mult. (111->44), div. (0->0), fcn. (101->4), ass. (0->23)
t14 = cos(qJ(3));
t23 = 0.2e1 * t14;
t10 = t14 ^ 2;
t12 = sin(qJ(3));
t8 = t12 ^ 2;
t22 = t10 + t8;
t13 = sin(qJ(2));
t21 = t12 * t13;
t20 = t14 * t13;
t15 = cos(qJ(2));
t19 = t15 * t12;
t18 = -qJ(4) - pkin(5);
t2 = t22 * t13;
t3 = t18 * t12;
t4 = t18 * t14;
t17 = -t3 * t12 - t4 * t14;
t11 = t15 ^ 2;
t9 = t13 ^ 2;
t7 = -t14 * pkin(3) - pkin(2);
t6 = t15 * t14;
t5 = t12 * t23;
t1 = t22 * t9 + t11;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 + t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t19, t2, t15 * pkin(2) + pkin(5) * t2, 0, 0, 0, 0, 0, 0, t6, -t19, t2, t17 * t13 - t15 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t8, t5, 0, t10, 0, 0, pkin(2) * t23, -0.2e1 * pkin(2) * t12, 0.2e1 * t22 * pkin(5), t22 * pkin(5) ^ 2 + pkin(2) ^ 2, t8, t5, 0, t10, 0, 0, -0.2e1 * t7 * t14, 0.2e1 * t7 * t12, 0.2e1 * t17, t3 ^ 2 + t4 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t20, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t20, 0, -pkin(3) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t14, 0, -t12 * pkin(5), -t14 * pkin(5), 0, 0, 0, 0, t12, 0, t14, 0, t3, t4, -t12 * pkin(3), t3 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(3), 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t12, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t16;
