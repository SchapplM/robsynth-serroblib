% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:52:40
% EndTime: 2019-05-04 21:52:41
% DurationCPUTime: 0.32s
% Computational Cost: add. (121->51), mult. (297->93), div. (0->0), fcn. (357->10), ass. (0->44)
t24 = sin(pkin(11));
t14 = t24 * pkin(2) + qJ(4);
t49 = 0.2e1 * t14;
t31 = cos(qJ(6));
t48 = 0.2e1 * t31;
t26 = cos(pkin(11));
t16 = -t26 * pkin(2) - pkin(3);
t13 = -pkin(8) + t16;
t32 = cos(qJ(5));
t23 = t32 ^ 2;
t47 = t13 * t23;
t25 = sin(pkin(6));
t30 = sin(qJ(2));
t46 = t25 * t30;
t33 = cos(qJ(2));
t45 = t25 * t33;
t28 = sin(qJ(6));
t44 = t28 * t31;
t43 = t28 * t32;
t29 = sin(qJ(5));
t42 = t29 * t13;
t41 = t31 * t29;
t18 = t31 * t32;
t40 = t32 * t13;
t39 = t32 * t29;
t21 = t29 ^ 2;
t38 = -t21 - t23;
t37 = -0.2e1 * t39;
t27 = cos(pkin(6));
t7 = t24 * t46 - t26 * t45;
t9 = (-t24 * t33 - t26 * t30) * t25;
t36 = t27 ^ 2 + t7 ^ 2 + t9 ^ 2;
t35 = -pkin(5) * t32 - pkin(9) * t29;
t22 = t31 ^ 2;
t20 = t28 ^ 2;
t17 = t28 * t29;
t11 = t29 * pkin(5) - t32 * pkin(9) + t14;
t6 = t27 * t32 + t7 * t29;
t5 = t27 * t29 - t7 * t32;
t4 = t28 * t11 + t13 * t41;
t3 = t31 * t11 - t28 * t42;
t2 = -t9 * t28 + t6 * t31;
t1 = -t6 * t28 - t9 * t31;
t8 = [1, 0, 0, 0, t36, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t45, -t46 (-t24 * t9 - t26 * t7) * pkin(2), t7, -t9, -t9 * t14 + t7 * t16, 0, 0, 0, 0, 0, -t9 * t29, -t9 * t32, 0, 0, 0, 0, 0, t1 * t29 + t5 * t43, t5 * t18 - t2 * t29; 0, 1, 0, 0 (t24 ^ 2 + t26 ^ 2) * pkin(2) ^ 2, 0.2e1 * t16, t49, t14 ^ 2 + t16 ^ 2, t23, t37, 0, 0, 0, t29 * t49, t32 * t49, t22 * t23, -0.2e1 * t23 * t44, t39 * t48, t28 * t37, t21, -0.2e1 * t28 * t47 + 0.2e1 * t3 * t29, -0.2e1 * t4 * t29 - 0.2e1 * t31 * t47; 0, 0, 0, 0, t27, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t28, t38 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t31, t5 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t29, 0, t40, -t42, t28 * t18 (-t20 + t22) * t32, t17, t41, 0, t35 * t28 + t31 * t40, -t28 * t40 + t35 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t32, 0, 0, 0, 0, 0, -t41, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t29, 0, 0, 0, 0, 0, t18, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t20, 0.2e1 * t44, 0, 0, 0, pkin(5) * t48, -0.2e1 * pkin(5) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t43, t29, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t31, 0, -t28 * pkin(9), -t31 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
