% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR4
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
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t25 = cos(qJ(4));
t18 = -t25 * pkin(4) - pkin(3);
t30 = cos(qJ(3)) * pkin(2);
t11 = t18 - t30;
t37 = 0.2e1 * t11;
t36 = 0.2e1 * t18;
t35 = 0.2e1 * t25;
t21 = sin(qJ(5));
t34 = t21 * pkin(4);
t33 = sin(qJ(3)) * pkin(2);
t24 = cos(qJ(5));
t32 = t24 * pkin(4);
t31 = t25 * pkin(7);
t17 = -pkin(3) - t30;
t29 = pkin(3) - t17;
t16 = pkin(7) + t33;
t28 = t25 * t16;
t27 = t11 + t18;
t22 = sin(qJ(4));
t20 = t22 ^ 2;
t19 = t25 * pkin(8);
t14 = t22 * t35;
t13 = t19 + t31;
t12 = (-pkin(7) - pkin(8)) * t22;
t10 = t21 * t25 + t24 * t22;
t9 = t21 * t22 - t24 * t25;
t8 = t10 ^ 2;
t7 = t19 + t28;
t6 = (-pkin(8) - t16) * t22;
t5 = -t21 * t12 - t24 * t13;
t4 = t24 * t12 - t21 * t13;
t3 = -0.2e1 * t10 * t9;
t2 = -t21 * t6 - t24 * t7;
t1 = -t21 * t7 + t24 * t6;
t15 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 1, 0.2e1 * t30, -0.2e1 * t33, t20, t14, 0, 0, 0, -0.2e1 * t17 * t25, 0.2e1 * t17 * t22, t8, t3, 0, 0, 0, t9 * t37, t10 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t30, -t33, t20, t14, 0, 0, 0, t29 * t25, -t29 * t22, t8, t3, 0, 0, 0, t27 * t9, t27 * t10; 0, 0, 0, 0, 1, 0, 0, t20, t14, 0, 0, 0, pkin(3) * t35, -0.2e1 * pkin(3) * t22, t8, t3, 0, 0, 0, t9 * t36, t10 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t22, 0, 0, 0, 0, 0, -t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t25, 0, -t22 * t16, -t28, 0, 0, t10, -t9, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t25, 0, -t22 * pkin(7), -t31, 0, 0, t10, -t9, 0, t4, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t32, -0.2e1 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t4, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t32, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t15;
