% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:48
% EndTime: 2019-12-31 18:02:49
% DurationCPUTime: 0.23s
% Computational Cost: add. (107->44), mult. (176->83), div. (0->0), fcn. (186->6), ass. (0->33)
t21 = sin(qJ(4));
t36 = -0.2e1 * t21;
t23 = cos(qJ(4));
t35 = 0.2e1 * t23;
t22 = cos(qJ(5));
t34 = pkin(4) * t22;
t16 = t21 ^ 2;
t20 = sin(qJ(5));
t33 = t16 * t20;
t32 = t16 * t22;
t31 = t20 * t21;
t30 = t20 * t22;
t29 = t20 * t23;
t18 = sin(pkin(8));
t28 = t21 * t18;
t27 = t22 * t21;
t13 = t22 * t23;
t26 = t23 * t18;
t19 = cos(pkin(8));
t24 = -pkin(1) - pkin(2);
t10 = t19 * qJ(2) + t18 * t24;
t25 = t21 * t35;
t8 = t18 * qJ(2) - t19 * t24;
t6 = pkin(3) + t8;
t17 = t22 ^ 2;
t15 = t20 ^ 2;
t7 = -pkin(6) + t10;
t5 = -t20 * t19 + t22 * t26;
t4 = -t22 * t19 - t20 * t26;
t3 = t23 * pkin(4) + t21 * pkin(7) + t6;
t2 = t7 * t13 + t20 * t3;
t1 = t22 * t3 - t7 * t29;
t9 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 0.2e1 * t8, 0.2e1 * t10, t10 ^ 2 + t8 ^ 2, t16, t25, 0, 0, 0, t6 * t35, t6 * t36, t17 * t16, -0.2e1 * t16 * t30, t13 * t36, t20 * t25, t23 ^ 2, 0.2e1 * t1 * t23 - 0.2e1 * t7 * t33, -0.2e1 * t2 * t23 - 0.2e1 * t7 * t32; 0, 0, 0, -1, 0, -pkin(1), -t19, t18, t10 * t18 - t8 * t19, 0, 0, 0, 0, 0, -t19 * t23, t21 * t19, 0, 0, 0, 0, 0, -t18 * t33 + t4 * t23, -t18 * t32 - t5 * t23; 0, 0, 0, 0, 0, 1, 0, 0, t18 ^ 2 + t19 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t23, 0, -t21 * t7, -t23 * t7, -t20 * t27, (t15 - t17) * t21, t29, t13, 0, -t7 * t27 + (pkin(4) * t21 - pkin(7) * t23) * t20, -pkin(7) * t13 + (t20 * t7 + t34) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t26, 0, 0, 0, 0, 0, -t18 * t27, t20 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t21, 0, 0, 0, 0, 0, t13, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t15, 0.2e1 * t30, 0, 0, 0, 0.2e1 * t34, -0.2e1 * pkin(4) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t31, t23, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t22, 0, -t20 * pkin(7), -t22 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
