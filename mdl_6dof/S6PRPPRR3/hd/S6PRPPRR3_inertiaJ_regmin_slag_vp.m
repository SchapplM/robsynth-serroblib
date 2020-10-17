% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:01:35
% EndTime: 2019-05-04 22:01:37
% DurationCPUTime: 0.35s
% Computational Cost: add. (163->68), mult. (337->121), div. (0->0), fcn. (405->10), ass. (0->47)
t33 = sin(qJ(5));
t52 = -0.2e1 * t33;
t36 = cos(qJ(5));
t51 = 0.2e1 * t36;
t35 = cos(qJ(6));
t50 = pkin(5) * t35;
t26 = t33 ^ 2;
t32 = sin(qJ(6));
t49 = t26 * t32;
t48 = t26 * t35;
t29 = sin(pkin(6));
t34 = sin(qJ(2));
t47 = t29 * t34;
t46 = t32 * t33;
t45 = t32 * t35;
t44 = t32 * t36;
t28 = sin(pkin(11));
t43 = t33 * t28;
t42 = t35 * t33;
t21 = t35 * t36;
t30 = cos(pkin(11));
t38 = -pkin(2) - pkin(3);
t17 = t30 * qJ(3) + t28 * t38;
t14 = -pkin(8) + t17;
t41 = t36 * t14;
t40 = t36 * t28;
t39 = t33 * t51;
t15 = t28 * qJ(3) - t30 * t38;
t13 = pkin(4) + t15;
t37 = cos(qJ(2));
t31 = cos(pkin(6));
t27 = t35 ^ 2;
t25 = t32 ^ 2;
t24 = t31 ^ 2;
t18 = t29 * t37;
t12 = -t32 * t30 + t35 * t40;
t11 = -t35 * t30 - t32 * t40;
t10 = (-t28 * t37 + t30 * t34) * t29;
t8 = (-t28 * t34 - t30 * t37) * t29;
t7 = t36 * pkin(5) + t33 * pkin(9) + t13;
t6 = t10 * t36 - t31 * t33;
t5 = t10 * t33 + t31 * t36;
t4 = t32 * t7 + t35 * t41;
t3 = -t32 * t41 + t35 * t7;
t2 = -t8 * t32 + t6 * t35;
t1 = -t6 * t32 - t8 * t35;
t9 = [1, 0, 0, 0, 0, 0, t24 + (t34 ^ 2 + t37 ^ 2) * t29 ^ 2, 0, 0, t10 ^ 2 + t8 ^ 2 + t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t18, -t47, t18, t47 (pkin(2) * t37 + qJ(3) * t34) * t29, -t8, t10, t10 * t17 - t8 * t15, 0, 0, 0, 0, 0, -t8 * t36, t8 * t33, 0, 0, 0, 0, 0, t1 * t36 - t5 * t46, -t2 * t36 - t5 * t42; 0, 1, 0, 0, 0.2e1 * pkin(2), 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0.2e1 * t15, 0.2e1 * t17, t15 ^ 2 + t17 ^ 2, t26, t39, 0, 0, 0, t13 * t51, t13 * t52, t27 * t26, -0.2e1 * t26 * t45, t21 * t52, t32 * t39, t36 ^ 2, -0.2e1 * t14 * t49 + 0.2e1 * t3 * t36, -0.2e1 * t14 * t48 - 0.2e1 * t4 * t36; 0, 0, 0, 0, 0, 0, -t18, 0, 0, t10 * t28 + t8 * t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -1, 0, -pkin(2), -t30, t28, -t15 * t30 + t17 * t28, 0, 0, 0, 0, 0, -t30 * t36, t33 * t30, 0, 0, 0, 0, 0, t11 * t36 - t28 * t49, -t12 * t36 - t28 * t48; 0, 0, 0, 0, 0, 0, 1, 0, 0, t28 ^ 2 + t30 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t35, t5 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t36, 0, -t33 * t14, -t41, -t32 * t42 (t25 - t27) * t33, t44, t21, 0, -t14 * t42 + (pkin(5) * t33 - pkin(9) * t36) * t32, -pkin(9) * t21 + (t14 * t32 + t50) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t40, 0, 0, 0, 0, 0, -t28 * t42, t32 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t33, 0, 0, 0, 0, 0, t21, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t25, 0.2e1 * t45, 0, 0, 0, 0.2e1 * t50, -0.2e1 * pkin(5) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t46, t36, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t35, 0, -t32 * pkin(9), -t35 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
