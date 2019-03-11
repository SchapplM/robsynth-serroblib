% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t32 = sin(qJ(2));
t58 = -0.2e1 * t32;
t34 = cos(qJ(2));
t57 = -0.2e1 * t34;
t56 = 0.2e1 * t34;
t35 = -pkin(2) - pkin(3);
t31 = sin(qJ(5));
t55 = t31 * pkin(5);
t54 = t34 * pkin(7);
t26 = t31 ^ 2;
t53 = t26 * t34;
t52 = t31 * t32;
t33 = cos(qJ(5));
t51 = t31 * t33;
t50 = t31 * t34;
t21 = t32 * pkin(7);
t13 = -t32 * qJ(4) + t21;
t49 = t33 * t13;
t48 = t33 * t32;
t47 = t33 * t34;
t30 = qJ(3) + pkin(4);
t27 = t32 ^ 2;
t29 = t34 ^ 2;
t46 = t27 + t29;
t45 = qJ(6) * t34;
t44 = t34 * qJ(3);
t25 = -pkin(8) + t35;
t43 = qJ(6) - t25;
t42 = t32 * t56;
t12 = -t34 * pkin(2) - t32 * qJ(3) - pkin(1);
t9 = t34 * pkin(3) - t12;
t7 = t32 * pkin(4) + t34 * pkin(8) + t9;
t3 = -t31 * t13 + t33 * t7;
t1 = t32 * pkin(5) + t33 * t45 + t3;
t2 = t49 + (t7 + t45) * t31;
t41 = t1 * t33 + t2 * t31;
t40 = -t32 * pkin(2) + t44;
t39 = -t25 * t32 - t30 * t34;
t37 = qJ(3) ^ 2;
t36 = 0.2e1 * qJ(3);
t28 = t33 ^ 2;
t22 = t33 * pkin(5);
t20 = t34 * qJ(4);
t18 = t28 * t34;
t17 = t22 + t30;
t16 = t28 + t26;
t14 = -t20 + t54;
t11 = t43 * t33;
t10 = t43 * t31;
t8 = -t20 + (pkin(7) - t55) * t34;
t5 = -t10 * t31 - t11 * t33;
t4 = t31 * t7 + t49;
t6 = [1, 0, 0, t27, t42, 0, 0, 0, pkin(1) * t56, pkin(1) * t58, t12 * t57, 0.2e1 * t46 * pkin(7), t12 * t58, t46 * pkin(7) ^ 2 + t12 ^ 2, 0.2e1 * t9 * t32, t9 * t57, -0.2e1 * t13 * t32 - 0.2e1 * t14 * t34, t13 ^ 2 + t14 ^ 2 + t9 ^ 2, t28 * t29, -0.2e1 * t29 * t51, t47 * t58, t31 * t42, t27, -0.2e1 * t14 * t50 + 0.2e1 * t3 * t32, -0.2e1 * t14 * t47 - 0.2e1 * t4 * t32, t41 * t56, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, t32, t34, 0, -t21, -t54, -t21, t40, t54, t40 * pkin(7), t14, t13, -t35 * t32 - t44, t14 * qJ(3) + t13 * t35, t31 * t47, t18 - t53, -t52, -t48, 0, t14 * t33 + t39 * t31, -t14 * t31 + t39 * t33 (t10 * t34 - t2) * t33 + (-t11 * t34 + t1) * t31, t1 * t10 - t2 * t11 + t8 * t17; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, t36, pkin(2) ^ 2 + t37, t36, 0.2e1 * t35, 0, t35 ^ 2 + t37, t26, 0.2e1 * t51, 0, 0, 0, 0.2e1 * t30 * t33, -0.2e1 * t30 * t31, -0.2e1 * t5, t10 ^ 2 + t11 ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, t21, 0, 0, -t32, t13, 0, 0, 0, 0, 0, -t52, -t48, 0, -t1 * t31 + t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 1, 0, t35, 0, 0, 0, 0, 0, 0, 0, -t16, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t34, 0, t9, 0, 0, 0, 0, 0, t48, -t52, t18 + t53, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t33 - t11 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t50, t32, t3, -t4, pkin(5) * t47, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33, 0, -t31 * t25, -t33 * t25, t55, t10 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33, 0, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
