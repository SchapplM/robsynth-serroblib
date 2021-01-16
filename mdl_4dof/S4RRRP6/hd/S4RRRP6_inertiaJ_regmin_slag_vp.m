% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRP6
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
% MM_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:14
% EndTime: 2021-01-15 14:39:16
% DurationCPUTime: 0.21s
% Computational Cost: add. (133->54), mult. (295->107), div. (0->0), fcn. (273->4), ass. (0->33)
t16 = sin(qJ(2));
t32 = -0.2e1 * t16;
t18 = cos(qJ(2));
t31 = 0.2e1 * t18;
t17 = cos(qJ(3));
t30 = pkin(2) * t17;
t15 = sin(qJ(3));
t29 = pkin(5) * t15;
t28 = t15 * pkin(3);
t27 = t15 * t16;
t26 = t15 * t17;
t25 = t15 * t18;
t10 = t17 * t16;
t24 = t17 * t18;
t23 = qJ(4) + pkin(6);
t22 = qJ(4) * t16;
t21 = t16 * t31;
t20 = pkin(5) * t24;
t7 = -t18 * pkin(2) - t16 * pkin(6) - pkin(1);
t5 = t17 * t7;
t19 = -t17 * t22 + t5;
t14 = t17 ^ 2;
t13 = t16 ^ 2;
t12 = t15 ^ 2;
t11 = -t17 * pkin(3) - pkin(2);
t9 = t23 * t17;
t8 = t23 * t15;
t6 = (pkin(5) + t28) * t16;
t4 = t15 * t7 + t20;
t3 = -pkin(5) * t25 + t5;
t2 = t20 + (t7 - t22) * t15;
t1 = (-pkin(3) - t29) * t18 + t19;
t33 = [1, 0, 0, t13, t21, 0, 0, 0, pkin(1) * t31, pkin(1) * t32, t14 * t13, -0.2e1 * t13 * t26, t24 * t32, t15 * t21, t18 ^ 2, 0.2e1 * t13 * t29 - 0.2e1 * t3 * t18, 0.2e1 * t13 * pkin(5) * t17 + 0.2e1 * t4 * t18, -0.2e1 * t1 * t18 + 0.2e1 * t6 * t27, 0.2e1 * t6 * t10 + 0.2e1 * t2 * t18, 0.2e1 * (-t1 * t17 - t15 * t2) * t16, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t16, t18, 0, -t16 * pkin(5), -t18 * pkin(5), t15 * t10, (-t12 + t14) * t16, -t25, -t24, 0, -pkin(5) * t10 + (-pkin(2) * t16 + pkin(6) * t18) * t15, pkin(6) * t24 + (t29 - t30) * t16, t11 * t27 - t6 * t17 + t8 * t18, t11 * t10 + t6 * t15 + t9 * t18, (t16 * t8 + t2) * t17 + (-t16 * t9 - t1) * t15, -t1 * t8 + t6 * t11 + t2 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t12, 0.2e1 * t26, 0, 0, 0, 0.2e1 * t30, -0.2e1 * pkin(2) * t15, -0.2e1 * t11 * t17, 0.2e1 * t11 * t15, 0.2e1 * t8 * t15 + 0.2e1 * t9 * t17, t11 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t27, -t18, t3, -t4, (-0.2e1 * pkin(3) - t29) * t18 + t19, -t2, -pkin(3) * t10, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t17, 0, -t15 * pkin(6), -t17 * pkin(6), -t8, -t9, -t28, -t8 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t10, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t15, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t33;
