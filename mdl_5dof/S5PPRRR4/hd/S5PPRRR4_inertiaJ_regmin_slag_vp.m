% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:43
% EndTime: 2019-12-05 15:19:44
% DurationCPUTime: 0.22s
% Computational Cost: add. (112->52), mult. (343->103), div. (0->0), fcn. (442->12), ass. (0->43)
t26 = sin(qJ(4));
t43 = -0.2e1 * t26;
t29 = cos(qJ(4));
t42 = 0.2e1 * t29;
t28 = cos(qJ(5));
t41 = pkin(4) * t28;
t25 = sin(qJ(5));
t40 = pkin(8) * t25;
t20 = sin(pkin(6));
t27 = sin(qJ(3));
t39 = t20 * t27;
t30 = cos(qJ(3));
t38 = t20 * t30;
t22 = cos(pkin(11));
t23 = cos(pkin(6));
t37 = t22 * t23;
t36 = t25 * t26;
t35 = t25 * t28;
t34 = t25 * t29;
t33 = t28 * t26;
t32 = t28 * t29;
t31 = t26 * t42;
t24 = cos(pkin(5));
t21 = sin(pkin(5));
t19 = sin(pkin(11));
t18 = t28 ^ 2;
t17 = t26 ^ 2;
t16 = t25 ^ 2;
t14 = -t29 * pkin(4) - t26 * pkin(9) - pkin(3);
t13 = t26 * t23 + t29 * t39;
t12 = -t29 * t23 + t26 * t39;
t11 = -t21 * t22 * t20 + t24 * t23;
t10 = pkin(8) * t32 + t25 * t14;
t9 = -pkin(8) * t34 + t28 * t14;
t8 = t28 * t13 - t25 * t38;
t7 = -t25 * t13 - t28 * t38;
t6 = t24 * t39 + (t19 * t30 + t27 * t37) * t21;
t5 = -t24 * t38 + (t19 * t27 - t30 * t37) * t21;
t4 = t11 * t26 + t6 * t29;
t3 = -t11 * t29 + t6 * t26;
t2 = t5 * t25 + t4 * t28;
t1 = -t4 * t25 + t5 * t28;
t15 = [1, t24 ^ 2 + (t19 ^ 2 + t22 ^ 2) * t21 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t29, t5 * t26, 0, 0, 0, 0, 0, -t1 * t29 + t3 * t36, t2 * t29 + t3 * t33; 0, 0, 0, t38, -t39, 0, 0, 0, 0, 0, t29 * t38, -t26 * t38, 0, 0, 0, 0, 0, t12 * t36 - t7 * t29, t12 * t33 + t8 * t29; 0, 0, 1, 0, 0, t17, t31, 0, 0, 0, pkin(3) * t42, pkin(3) * t43, t18 * t17, -0.2e1 * t17 * t35, t32 * t43, t25 * t31, t29 ^ 2, 0.2e1 * t17 * t40 - 0.2e1 * t9 * t29, 0.2e1 * t17 * pkin(8) * t28 + 0.2e1 * t10 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, -t3 * t28, t3 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, -t12 * t28, t12 * t25; 0, 0, 0, 0, 0, 0, 0, t26, t29, 0, -t26 * pkin(8), -t29 * pkin(8), t25 * t33, (-t16 + t18) * t26, -t34, -t32, 0, -pkin(8) * t33 + (-pkin(4) * t26 + pkin(9) * t29) * t25, pkin(9) * t32 + (t40 - t41) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t16, 0.2e1 * t35, 0, 0, 0, 0.2e1 * t41, -0.2e1 * pkin(4) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t36, -t29, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t28, 0, -t25 * pkin(9), -t28 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t15;
