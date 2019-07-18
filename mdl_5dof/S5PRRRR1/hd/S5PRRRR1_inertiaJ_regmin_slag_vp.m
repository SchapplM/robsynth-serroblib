% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_inertiaJ_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t19 = sin(qJ(4));
t20 = sin(qJ(3));
t23 = cos(qJ(4));
t24 = cos(qJ(3));
t12 = t19 * t24 + t23 * t20;
t41 = -0.2e1 * t12;
t40 = pkin(2) * t24;
t39 = t19 * pkin(2);
t38 = t23 * pkin(2);
t22 = cos(qJ(5));
t21 = sin(qJ(2));
t6 = t12 * t21;
t37 = t6 * t22;
t18 = sin(qJ(5));
t36 = t18 * t12;
t35 = t18 * t22;
t34 = t20 * t21;
t33 = t22 * t12;
t32 = t24 * t21;
t11 = t19 * t20 - t23 * t24;
t31 = t11 * t41;
t30 = t18 * t38;
t29 = t22 * t38;
t28 = t11 * t40;
t27 = -0.2e1 * t28;
t26 = pkin(2) * (-t11 * t19 - t12 * t23);
t25 = cos(qJ(2));
t17 = t22 ^ 2;
t16 = t18 ^ 2;
t14 = 0.2e1 * t35;
t10 = t12 ^ 2;
t9 = t22 * t11;
t8 = t18 * t11;
t7 = -t19 * t34 + t23 * t32;
t5 = t18 * t33;
t4 = t6 * t18;
t3 = -t25 * t18 + t22 * t7;
t2 = -t18 * t7 - t25 * t22;
t1 = (-t16 + t17) * t12;
t13 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t25, -t21, 0, 0, 0, 0, 0, t25 * t24, -t25 * t20, 0, 0, 0, 0, 0, -t25 * t11, -t25 * t12, 0, 0, 0, 0, 0, t2 * t11 + t6 * t36, -t3 * t11 + t6 * t33; 0, 1, 0, 0, t20 ^ 2, 0.2e1 * t20 * t24, 0, 0, 0, 0, 0, t10, t31, 0, 0, 0, t27, t40 * t41, t17 * t10, -0.2e1 * t10 * t35, 0.2e1 * t11 * t33, t18 * t31, t11 ^ 2, t22 * t27, 0.2e1 * t18 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t32, 0, 0, 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t37, t4; 0, 0, 0, 0, 0, 0, t20, t24, 0, 0, 0, 0, 0, t12, -t11, 0, 0, 0, t5, t1, t8, t9, 0, t18 * t26, t22 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t39, t16, t14, 0, 0, 0, 0.2e1 * t29, -0.2e1 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t37, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, 0, 0, t5, t1, t8, t9, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, -t39, t16, t14, 0, 0, 0, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t16, t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t36, t11, -t22 * t40, t18 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t22, 0, -t18 * t39, -t22 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t22, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
