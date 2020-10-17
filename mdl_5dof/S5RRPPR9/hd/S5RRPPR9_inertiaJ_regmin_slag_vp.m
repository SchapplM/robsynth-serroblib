% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:47
% EndTime: 2019-12-31 19:41:48
% DurationCPUTime: 0.26s
% Computational Cost: add. (114->48), mult. (196->85), div. (0->0), fcn. (177->4), ass. (0->36)
t21 = sin(qJ(2));
t41 = -0.2e1 * t21;
t23 = cos(qJ(2));
t40 = -0.2e1 * t23;
t39 = 0.2e1 * t23;
t24 = -pkin(2) - pkin(3);
t38 = t23 * pkin(6);
t20 = sin(qJ(5));
t37 = t20 * t21;
t22 = cos(qJ(5));
t36 = t20 * t22;
t35 = t20 * t23;
t34 = t22 * t21;
t33 = t22 * t23;
t16 = t21 ^ 2;
t18 = t23 ^ 2;
t32 = t16 + t18;
t31 = t23 * qJ(3);
t30 = t21 * t39;
t5 = -t23 * pkin(2) - t21 * qJ(3) - pkin(1);
t4 = t23 * pkin(3) - t5;
t29 = -t21 * pkin(2) + t31;
t14 = -pkin(7) + t24;
t19 = qJ(3) + pkin(4);
t28 = -t14 * t21 - t19 * t23;
t26 = qJ(3) ^ 2;
t25 = 0.2e1 * qJ(3);
t17 = t22 ^ 2;
t15 = t20 ^ 2;
t11 = t21 * pkin(6);
t7 = -t23 * qJ(4) + t38;
t6 = -t21 * qJ(4) + t11;
t3 = t21 * pkin(4) + t23 * pkin(7) + t4;
t2 = t20 * t3 + t22 * t6;
t1 = -t20 * t6 + t22 * t3;
t8 = [1, 0, 0, t16, t30, 0, 0, 0, pkin(1) * t39, pkin(1) * t41, t5 * t40, 0.2e1 * t32 * pkin(6), t5 * t41, t32 * pkin(6) ^ 2 + t5 ^ 2, 0.2e1 * t4 * t21, t4 * t40, -0.2e1 * t6 * t21 - 0.2e1 * t7 * t23, t4 ^ 2 + t6 ^ 2 + t7 ^ 2, t17 * t18, -0.2e1 * t18 * t36, t33 * t41, t20 * t30, t16, 0.2e1 * t1 * t21 - 0.2e1 * t7 * t35, -0.2e1 * t2 * t21 - 0.2e1 * t7 * t33; 0, 0, 0, 0, 0, t21, t23, 0, -t11, -t38, -t11, t29, t38, t29 * pkin(6), t7, t6, -t24 * t21 - t31, t7 * qJ(3) + t6 * t24, t20 * t33, (-t15 + t17) * t23, -t37, -t34, 0, t28 * t20 + t7 * t22, -t7 * t20 + t28 * t22; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, t25, pkin(2) ^ 2 + t26, t25, 0.2e1 * t24, 0, t24 ^ 2 + t26, t15, 0.2e1 * t36, 0, 0, 0, 0.2e1 * t19 * t22, -0.2e1 * t19 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, t11, 0, 0, -t21, t6, 0, 0, 0, 0, 0, -t37, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 1, 0, t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t23, 0, t4, 0, 0, 0, 0, 0, t34, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t35, t21, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t22, 0, -t20 * t14, -t22 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
