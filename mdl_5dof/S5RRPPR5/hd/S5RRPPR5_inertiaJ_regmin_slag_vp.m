% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:30
% EndTime: 2021-01-15 19:36:32
% DurationCPUTime: 0.26s
% Computational Cost: add. (224->48), mult. (429->82), div. (0->0), fcn. (480->6), ass. (0->35)
t27 = sin(pkin(8));
t28 = cos(pkin(8));
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t17 = t27 * t30 - t28 * t32;
t18 = t27 * t32 + t28 * t30;
t26 = -t32 * pkin(2) - pkin(1);
t34 = t18 * qJ(4) - t26;
t43 = 0.2e1 * (-pkin(3) - pkin(4)) * t17 + 0.2e1 * t34;
t42 = 0.2e1 * t26;
t41 = 0.2e1 * t32;
t40 = t27 * pkin(2);
t39 = t28 * pkin(2);
t38 = -qJ(3) - pkin(6);
t20 = t38 * t32;
t36 = t38 * t30;
t11 = -t27 * t20 - t28 * t36;
t13 = -t28 * t20 + t27 * t36;
t37 = t11 ^ 2 + t13 ^ 2;
t24 = pkin(3) + t39;
t35 = 0.2e1 * t11 * t18 - 0.2e1 * t13 * t17;
t31 = cos(qJ(5));
t29 = sin(qJ(5));
t22 = qJ(4) + t40;
t21 = -pkin(4) - t24;
t15 = t29 * t21 + t31 * t22;
t14 = -t31 * t21 + t29 * t22;
t9 = t29 * t17 + t31 * t18;
t8 = -t31 * t17 + t29 * t18;
t7 = t17 * pkin(3) - t34;
t5 = t17 * pkin(7) + t13;
t4 = -t18 * pkin(7) + t11;
t2 = t29 * t4 + t31 * t5;
t1 = t29 * t5 - t31 * t4;
t3 = [1, 0, 0, t30 ^ 2, t30 * t41, 0, 0, 0, pkin(1) * t41, -0.2e1 * pkin(1) * t30, t17 * t42, t18 * t42, t35, t26 ^ 2 + t37, 0.2e1 * t7 * t17, t35, -0.2e1 * t7 * t18, t7 ^ 2 + t37, t9 ^ 2, -0.2e1 * t9 * t8, 0, 0, 0, t8 * t43, t9 * t43; 0, 0, 0, 0, 0, t30, t32, 0, -t30 * pkin(6), -t32 * pkin(6), -t11, -t13, (-t17 * t27 - t18 * t28) * pkin(2), (-t11 * t28 + t13 * t27) * pkin(2), -t11, -t22 * t17 - t24 * t18, t13, -t11 * t24 + t13 * t22, 0, 0, -t9, t8, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t39, -0.2e1 * t40, 0, (t27 ^ 2 + t28 ^ 2) * pkin(2) ^ 2, 0.2e1 * t24, 0, 0.2e1 * t22, t22 ^ 2 + t24 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t14, 0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, t26, t17, 0, -t18, t7, 0, 0, 0, 0, 0, -t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t24, 0, 0, 0, 0, 0, -t31, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
