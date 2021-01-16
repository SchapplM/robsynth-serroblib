% Calculate minimal parameter regressor of joint inertia matrix for
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
% MM_reg [((4+1)*4/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:29
% EndTime: 2021-01-15 11:04:30
% DurationCPUTime: 0.13s
% Computational Cost: add. (84->33), mult. (152->54), div. (0->0), fcn. (133->4), ass. (0->25)
t13 = sin(qJ(3));
t26 = 0.2e1 * t13;
t15 = cos(qJ(3));
t25 = -0.2e1 * t15;
t24 = 0.2e1 * t15;
t23 = t13 * pkin(3);
t22 = sin(qJ(2)) * pkin(1);
t21 = cos(qJ(2)) * pkin(1);
t10 = -pkin(2) - t21;
t20 = pkin(2) - t10;
t11 = -t15 * pkin(3) - pkin(2);
t5 = t11 - t21;
t19 = t11 + t5;
t18 = -qJ(4) - pkin(6);
t9 = pkin(6) + t22;
t17 = qJ(4) + t9;
t12 = t13 ^ 2;
t8 = t13 * t24;
t7 = t18 * t15;
t6 = t18 * t13;
t4 = t7 * t15;
t3 = t17 * t15;
t2 = t17 * t13;
t1 = t3 * t15;
t14 = [1, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t22, t12, t8, 0, 0, 0, t10 * t25, t10 * t26, t5 * t25, t5 * t26, 0.2e1 * t2 * t13 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t5 ^ 2; 0, 0, 0, 1, t21, -t22, t12, t8, 0, 0, 0, t20 * t15, -t20 * t13, -t19 * t15, t19 * t13, t1 - t4 + (t2 - t6) * t13, t5 * t11 - t2 * t6 - t3 * t7; 0, 0, 0, 1, 0, 0, t12, t8, 0, 0, 0, pkin(2) * t24, -0.2e1 * pkin(2) * t13, t11 * t25, t11 * t26, -0.2e1 * t6 * t13 - 0.2e1 * t4, t11 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, -t13 * t9, -t15 * t9, -t2, -t3, -t23, -t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, -t13 * pkin(6), -t15 * pkin(6), t6, t7, -t23, t6 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t13, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t13, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t14;
