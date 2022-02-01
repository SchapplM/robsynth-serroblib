% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:49:02
% EndTime: 2022-01-20 09:49:03
% DurationCPUTime: 0.24s
% Computational Cost: add. (134->35), mult. (264->58), div. (0->0), fcn. (293->8), ass. (0->42)
t31 = cos(qJ(4));
t37 = t31 * pkin(4);
t26 = cos(pkin(9));
t20 = t26 * pkin(1) + pkin(2);
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t25 = sin(pkin(9));
t40 = pkin(1) * t25;
t11 = t32 * t20 - t29 * t40;
t9 = -pkin(3) - t11;
t8 = t9 - t37;
t44 = 0.2e1 * t8;
t22 = -pkin(3) - t37;
t43 = 0.2e1 * t22;
t42 = 0.2e1 * t31;
t41 = pkin(3) - t9;
t27 = sin(qJ(5));
t39 = t27 * pkin(4);
t30 = cos(qJ(5));
t38 = t30 * pkin(4);
t36 = t31 * pkin(7);
t35 = t22 + t8;
t12 = -t29 * t20 - t32 * t40;
t10 = pkin(7) - t12;
t34 = t31 * t10;
t28 = sin(qJ(4));
t24 = t28 ^ 2;
t23 = t31 * pkin(8);
t19 = t28 * t42;
t17 = t23 + t36;
t16 = (-pkin(7) - pkin(8)) * t28;
t15 = t27 * t31 + t30 * t28;
t14 = t27 * t28 - t30 * t31;
t13 = t15 ^ 2;
t7 = -t27 * t16 - t30 * t17;
t6 = t30 * t16 - t27 * t17;
t5 = t23 + t34;
t4 = (-pkin(8) - t10) * t28;
t3 = -0.2e1 * t15 * t14;
t2 = -t27 * t4 - t30 * t5;
t1 = -t27 * t5 + t30 * t4;
t18 = [1, 0, 0, (t25 ^ 2 + t26 ^ 2) * pkin(1) ^ 2, 1, 0.2e1 * t11, 0.2e1 * t12, t24, t19, 0, 0, 0, -0.2e1 * t9 * t31, 0.2e1 * t9 * t28, t13, t3, 0, 0, 0, t14 * t44, t15 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t11, t12, t24, t19, 0, 0, 0, t41 * t31, -t41 * t28, t13, t3, 0, 0, 0, t35 * t14, t35 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, t24, t19, 0, 0, 0, pkin(3) * t42, -0.2e1 * pkin(3) * t28, t13, t3, 0, 0, 0, t14 * t43, t15 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t31, 0, -t28 * t10, -t34, 0, 0, t15, -t14, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t28, 0, 0, 0, 0, 0, -t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t31, 0, -t28 * pkin(7), -t36, 0, 0, t15, -t14, 0, t6, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, t6, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t18;
