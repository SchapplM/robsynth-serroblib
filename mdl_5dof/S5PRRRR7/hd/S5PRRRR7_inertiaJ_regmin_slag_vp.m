% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:59
% EndTime: 2019-12-05 17:13:00
% DurationCPUTime: 0.27s
% Computational Cost: add. (184->43), mult. (403->73), div. (0->0), fcn. (507->8), ass. (0->40)
t30 = sin(qJ(4));
t31 = sin(qJ(3));
t34 = cos(qJ(4));
t35 = cos(qJ(3));
t17 = t30 * t31 - t34 * t35;
t26 = -t35 * pkin(3) - pkin(2);
t45 = 0.2e1 * t17 * pkin(4) + 0.2e1 * t26;
t44 = 0.2e1 * t26;
t43 = 0.2e1 * t35;
t42 = pkin(6) + pkin(7);
t29 = sin(qJ(5));
t41 = t29 * pkin(4);
t40 = t30 * pkin(3);
t32 = sin(qJ(2));
t39 = t31 * t32;
t38 = t35 * t32;
t33 = cos(qJ(5));
t37 = t33 * t40;
t20 = t42 * t31;
t21 = t42 * t35;
t9 = -t34 * t20 - t30 * t21;
t28 = t34 * pkin(3);
t25 = t28 + pkin(4);
t14 = t33 * t25 - t29 * t40;
t10 = t30 * t20 - t34 * t21;
t18 = t30 * t35 + t34 * t31;
t36 = cos(qJ(2));
t27 = t33 * pkin(4);
t15 = -t29 * t25 - t37;
t13 = -t30 * t39 + t34 * t38;
t12 = t18 * t32;
t8 = -t29 * t17 + t33 * t18;
t7 = t33 * t17 + t29 * t18;
t6 = -t17 * pkin(8) - t10;
t5 = -t18 * pkin(8) + t9;
t4 = t29 * t12 - t33 * t13;
t3 = -t33 * t12 - t29 * t13;
t2 = -t29 * t5 - t33 * t6;
t1 = -t29 * t6 + t33 * t5;
t11 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t36, -t32, 0, 0, 0, 0, 0, t36 * t35, -t36 * t31, 0, 0, 0, 0, 0, -t36 * t17, -t36 * t18, 0, 0, 0, 0, 0, -t36 * t7, -t36 * t8; 0, 1, 0, 0, t31 ^ 2, t31 * t43, 0, 0, 0, pkin(2) * t43, -0.2e1 * pkin(2) * t31, t18 ^ 2, -0.2e1 * t18 * t17, 0, 0, 0, t17 * t44, t18 * t44, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t45, t8 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t38, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, t31, t35, 0, -t31 * pkin(6), -t35 * pkin(6), 0, 0, t18, -t17, 0, t9, t10, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t28, -0.2e1 * t40, 0, 0, 0, 0, 1, 0.2e1 * t14, 0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, t9, t10, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t28, -t40, 0, 0, 0, 0, 1, t14 + t27, -t37 + (-pkin(4) - t25) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t27, -0.2e1 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t11;
