% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:15
% EndTime: 2019-12-31 17:26:15
% DurationCPUTime: 0.20s
% Computational Cost: add. (141->38), mult. (323->82), div. (0->0), fcn. (375->6), ass. (0->42)
t27 = cos(qJ(2));
t20 = -t27 * pkin(2) - pkin(1);
t43 = 0.2e1 * t20;
t42 = 0.2e1 * t27;
t41 = pkin(5) + pkin(6);
t24 = sin(qJ(3));
t40 = t24 * pkin(2);
t26 = cos(qJ(4));
t15 = t41 * t27;
t25 = sin(qJ(2));
t37 = cos(qJ(3));
t30 = t37 * t25;
t6 = t24 * t15 + t41 * t30;
t39 = t6 * t26;
t31 = t37 * pkin(2);
t19 = -t31 - pkin(3);
t38 = pkin(3) - t19;
t13 = t24 * t27 + t30;
t23 = sin(qJ(4));
t36 = t23 * t13;
t35 = t23 * t26;
t34 = t24 * t25;
t33 = t26 * t13;
t12 = -t37 * t27 + t34;
t32 = -0.2e1 * t13 * t12;
t29 = -pkin(3) * t13 - pkin(7) * t12;
t18 = pkin(7) + t40;
t28 = -t12 * t18 + t13 * t19;
t22 = t26 ^ 2;
t21 = t23 ^ 2;
t16 = 0.2e1 * t35;
t11 = t13 ^ 2;
t10 = t26 * t12;
t9 = t23 * t12;
t8 = t23 * t33;
t7 = t37 * t15 - t41 * t34;
t5 = t6 * t23;
t4 = (-t21 + t22) * t13;
t3 = t12 * pkin(3) - t13 * pkin(7) + t20;
t2 = t23 * t3 + t26 * t7;
t1 = -t23 * t7 + t26 * t3;
t14 = [1, 0, 0, t25 ^ 2, t25 * t42, 0, 0, 0, pkin(1) * t42, -0.2e1 * pkin(1) * t25, t11, t32, 0, 0, 0, t12 * t43, t13 * t43, t22 * t11, -0.2e1 * t11 * t35, 0.2e1 * t12 * t33, t23 * t32, t12 ^ 2, 0.2e1 * t1 * t12 + 0.2e1 * t6 * t36, -0.2e1 * t2 * t12 + 0.2e1 * t6 * t33; 0, 0, 0, 0, 0, t25, t27, 0, -t25 * pkin(5), -t27 * pkin(5), 0, 0, t13, -t12, 0, -t6, -t7, t8, t4, t9, t10, 0, t28 * t23 - t39, t28 * t26 + t5; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t31, -0.2e1 * t40, t21, t16, 0, 0, 0, -0.2e1 * t19 * t26, 0.2e1 * t19 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, -t6, -t7, t8, t4, t9, t10, 0, t29 * t23 - t39, t29 * t26 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t31, -t40, t21, t16, 0, 0, 0, t38 * t26, -t38 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t21, t16, 0, 0, 0, 0.2e1 * pkin(3) * t26, -0.2e1 * pkin(3) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t36, t12, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t26, 0, -t23 * t18, -t26 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t26, 0, -t23 * pkin(7), -t26 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t14;
