% Calculate minimal parameter regressor of joint inertia matrix for
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
% MM_reg [((4+1)*4/2)x14]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:34
% DurationCPUTime: 0.10s
% Computational Cost: add. (43->18), mult. (89->33), div. (0->0), fcn. (81->6), ass. (0->19)
t15 = sin(qJ(4));
t23 = 0.2e1 * t15;
t17 = cos(qJ(4));
t22 = -0.2e1 * t17;
t11 = cos(qJ(2)) * pkin(1);
t10 = t11 + pkin(2);
t13 = sin(pkin(7));
t14 = cos(pkin(7));
t20 = sin(qJ(2)) * pkin(1);
t3 = t14 * t10 - t13 * t20;
t1 = -pkin(3) - t3;
t9 = -t14 * pkin(2) - pkin(3);
t21 = t1 + t9;
t4 = t13 * t10 + t14 * t20;
t12 = t15 ^ 2;
t8 = t13 * pkin(2) + pkin(6);
t7 = t17 * t23;
t2 = pkin(6) + t4;
t5 = [1, 0, 0, 1, 0.2e1 * t11, -0.2e1 * t20, t3 ^ 2 + t4 ^ 2, t12, t7, 0, 0, 0, t1 * t22, t1 * t23; 0, 0, 0, 1, t11, -t20, (t13 * t4 + t14 * t3) * pkin(2), t12, t7, 0, 0, 0, -t21 * t17, t21 * t15; 0, 0, 0, 1, 0, 0, (t13 ^ 2 + t14 ^ 2) * pkin(2) ^ 2, t12, t7, 0, 0, 0, t9 * t22, t9 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t17, 0, -t15 * t2, -t17 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t17, 0, -t15 * t8, -t17 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
