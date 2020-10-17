% Calculate inertial parameters regressor of joint inertia matrix for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:39
% EndTime: 2019-12-31 16:32:40
% DurationCPUTime: 0.19s
% Computational Cost: add. (85->29), mult. (201->58), div. (0->0), fcn. (217->4), ass. (0->22)
t17 = cos(qJ(3));
t11 = -t17 * pkin(3) - pkin(2);
t25 = 0.2e1 * t11;
t24 = 0.2e1 * t17;
t23 = -pkin(6) - pkin(5);
t14 = sin(qJ(4));
t22 = t14 * pkin(3);
t16 = cos(qJ(4));
t21 = t16 * pkin(3);
t15 = sin(qJ(3));
t12 = t15 ^ 2;
t13 = t17 ^ 2;
t20 = t12 + t13;
t9 = t23 * t17;
t8 = t23 * t15;
t7 = t14 * t17 + t16 * t15;
t5 = t14 * t15 - t16 * t17;
t4 = t7 ^ 2;
t3 = t5 ^ 2;
t2 = t14 * t8 - t16 * t9;
t1 = t14 * t9 + t16 * t8;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t1 + t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t12, t15 * t24, 0, t13, 0, 0, pkin(2) * t24, -0.2e1 * pkin(2) * t15, 0.2e1 * t20 * pkin(5), t20 * pkin(5) ^ 2 + pkin(2) ^ 2, t4, -0.2e1 * t7 * t5, 0, t3, 0, 0, t5 * t25, t7 * t25, -0.2e1 * t1 * t7 - 0.2e1 * t2 * t5, t1 ^ 2 + t11 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t7, 0, (t14 * t7 - t16 * t5) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, t17, 0, -t15 * pkin(5), -t17 * pkin(5), 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, (-t14 * t5 - t16 * t7) * pkin(3), (t1 * t16 + t14 * t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t22, 0, (t14 ^ 2 + t16 ^ 2) * pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t21, -t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
