% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:17
% EndTime: 2019-12-31 18:04:18
% DurationCPUTime: 0.21s
% Computational Cost: add. (132->34), mult. (265->59), div. (0->0), fcn. (314->6), ass. (0->32)
t27 = sin(pkin(8));
t28 = cos(pkin(8));
t42 = t27 ^ 2 + t28 ^ 2;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t11 = t27 * t30 + t28 * t32;
t16 = -t28 * pkin(2) - t27 * qJ(3) - pkin(1);
t9 = t28 * pkin(3) - t16;
t41 = 0.2e1 * t11 * pkin(4) + 0.2e1 * t9;
t40 = 0.2e1 * t9;
t39 = -0.2e1 * t27;
t29 = sin(qJ(5));
t38 = t29 * pkin(4);
t31 = cos(qJ(5));
t37 = t31 * pkin(4);
t36 = t42 * qJ(2) ^ 2;
t21 = t27 * qJ(2);
t17 = -t27 * pkin(6) + t21;
t18 = (-pkin(6) + qJ(2)) * t28;
t35 = t32 * t17 - t30 * t18;
t34 = -t30 * t17 - t32 * t18;
t15 = -t29 * t32 - t31 * t30;
t14 = -t29 * t30 + t31 * t32;
t13 = 0.2e1 * t42 * qJ(2);
t12 = t27 * t32 - t28 * t30;
t6 = -t29 * t11 + t31 * t12;
t5 = t31 * t11 + t29 * t12;
t4 = -t11 * pkin(7) - t34;
t3 = -t12 * pkin(7) + t35;
t2 = -t29 * t3 - t31 * t4;
t1 = -t29 * t4 + t31 * t3;
t7 = [1, 0, 0, 0.2e1 * pkin(1) * t28, pkin(1) * t39, t13, pkin(1) ^ 2 + t36, -0.2e1 * t16 * t28, t13, t16 * t39, t16 ^ 2 + t36, t12 ^ 2, -0.2e1 * t12 * t11, 0, 0, 0, t11 * t40, t12 * t40, t6 ^ 2, -0.2e1 * t6 * t5, 0, 0, 0, t5 * t41, t6 * t41; 0, 0, 0, -t28, t27, 0, -pkin(1), -t28, 0, -t27, t16, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t35, t34, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30, 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t7;
