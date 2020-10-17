% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR16_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:32
% EndTime: 2019-12-31 18:39:33
% DurationCPUTime: 0.24s
% Computational Cost: add. (90->39), mult. (162->71), div. (0->0), fcn. (153->4), ass. (0->34)
t21 = sin(qJ(3));
t23 = cos(qJ(3));
t5 = t21 * pkin(3) - t23 * qJ(4) + qJ(2);
t35 = -0.2e1 * t5;
t34 = 0.2e1 * qJ(2);
t33 = 0.2e1 * qJ(4);
t20 = sin(qJ(5));
t10 = t20 * t21;
t32 = t20 * t23;
t22 = cos(qJ(5));
t31 = t22 * t20;
t11 = t22 * t21;
t30 = t23 * t21;
t25 = -pkin(1) - pkin(6);
t29 = t23 * t25;
t16 = t21 ^ 2;
t18 = t23 ^ 2;
t9 = t16 + t18;
t28 = t21 * qJ(4);
t27 = 0.2e1 * t30;
t8 = t23 * pkin(3) + t28;
t24 = -pkin(3) - pkin(7);
t26 = -t23 * t24 + t28;
t17 = t22 ^ 2;
t15 = t20 ^ 2;
t13 = t21 * t25;
t12 = t22 * t23;
t7 = (pkin(4) - t25) * t23;
t6 = -t21 * pkin(4) + t13;
t4 = t9 * t25;
t3 = t21 * pkin(7) + t5;
t2 = t20 * t7 + t22 * t3;
t1 = -t20 * t3 + t22 * t7;
t14 = [1, 0, 0, -2 * pkin(1), t34, (pkin(1) ^ 2) + qJ(2) ^ 2, t18, -0.2e1 * t30, 0, 0, 0, t21 * t34, t23 * t34, -0.2e1 * t4, t21 * t35, t23 * t35, t9 * t25 ^ 2 + t5 ^ 2, t15 * t16, 0.2e1 * t16 * t31, t20 * t27, t22 * t27, t18, 0.2e1 * t1 * t23 - 0.2e1 * t6 * t11, 0.2e1 * t6 * t10 - 0.2e1 * t2 * t23; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, t4, 0, 0, 0, 0, 0, -t9 * t22, t9 * t20; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t23, -t21, 0, t29, -t13, -t8, -t29, t13, t8 * t25, t20 * t11, (-t15 + t17) * t21, t12, -t32, 0, t6 * t20 - t26 * t22, t26 * t20 + t6 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t21, 0, -t23, t21, t8, 0, 0, 0, 0, 0, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t33, pkin(3) ^ 2 + qJ(4) ^ 2, t17, -0.2e1 * t31, 0, 0, 0, t20 * t33, t22 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, -t29, 0, 0, 0, 0, 0, t12, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, t23, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, t22 * t24, -t20 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t14;
