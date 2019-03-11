% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t32 = sin(pkin(9));
t33 = cos(pkin(9));
t36 = cos(qJ(3));
t57 = sin(qJ(3));
t17 = -t32 * t57 + t33 * t36;
t62 = -0.2e1 * t17;
t18 = -t32 * t36 - t33 * t57;
t61 = (t17 * t33 - t18 * t32) * pkin(3);
t15 = t17 ^ 2;
t16 = t18 ^ 2;
t60 = t16 + t15;
t59 = 2 * qJ(2);
t35 = cos(qJ(5));
t58 = t35 * pkin(5);
t34 = sin(qJ(5));
t56 = t34 * t17;
t55 = t34 * t18;
t54 = t34 * t35;
t37 = -pkin(1) - pkin(7);
t47 = t57 * t37;
t21 = -t57 * qJ(4) + t47;
t46 = (-qJ(4) + t37) * t36;
t10 = t33 * t21 + t32 * t46;
t53 = t35 * t10;
t52 = t35 * t17;
t11 = t35 * t18;
t30 = t34 ^ 2;
t31 = t35 ^ 2;
t51 = t30 + t31;
t50 = qJ(6) * t17;
t28 = t57 * pkin(3) + qJ(2);
t26 = t32 * pkin(3) + pkin(8);
t49 = qJ(6) + t26;
t27 = -t33 * pkin(3) - pkin(4);
t7 = -t18 * pkin(4) - t17 * pkin(8) + t28;
t3 = -t34 * t10 + t35 * t7;
t8 = t32 * t21 - t33 * t46;
t1 = -t18 * pkin(5) - t35 * t50 + t3;
t2 = t53 + (t7 - t50) * t34;
t45 = t1 * t35 + t2 * t34;
t44 = t1 * t34 - t2 * t35;
t43 = t10 * t18 + t8 * t17;
t12 = t49 * t34;
t13 = t49 * t35;
t42 = -t12 * t35 + t13 * t34;
t41 = -t12 * t34 - t13 * t35;
t40 = t17 * t27 + t18 * t26;
t22 = t27 - t58;
t5 = pkin(5) * t56 + t8;
t4 = t34 * t7 + t53;
t6 = [1, 0, 0, -2 * pkin(1), t59, pkin(1) ^ 2 + qJ(2) ^ 2, t36 ^ 2, -0.2e1 * t36 * t57, 0, 0, 0, t57 * t59, t36 * t59, 0.2e1 * t43, t10 ^ 2 + t28 ^ 2 + t8 ^ 2, t31 * t15, -0.2e1 * t15 * t54, t11 * t62, 0.2e1 * t17 * t55, t16, -0.2e1 * t3 * t18 + 0.2e1 * t8 * t56, 0.2e1 * t4 * t18 + 0.2e1 * t8 * t52, t45 * t62, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t60, -t43, 0, 0, 0, 0, 0, -t60 * t34, -t60 * t35, 0, -t5 * t17 + t44 * t18; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t16 + t15; 0, 0, 0, 0, 0, 0, 0, 0, t36, -t57, 0, t36 * t37, -t47, -t61 (t10 * t32 - t33 * t8) * pkin(3), t34 * t52 (-t30 + t31) * t17, -t55, -t11, 0, t40 * t34 - t8 * t35, t8 * t34 + t40 * t35, -t42 * t17 - t44, -t1 * t12 + t2 * t13 + t5 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t57, 0, t61, 0, 0, 0, 0, 0, t52, -t56, -t51 * t18, -t17 * t22 + t41 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t32 ^ 2 + t33 ^ 2) * pkin(3) ^ 2, t30, 0.2e1 * t54, 0, 0, 0, -0.2e1 * t27 * t35, 0.2e1 * t27 * t34, -0.2e1 * t41, t12 ^ 2 + t13 ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, -t11, t55, -t51 * t17, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t56, -t18, t3, -t4, -pkin(5) * t52, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t11, 0, pkin(5) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t35, 0, -t34 * t26, -t35 * t26, -t34 * pkin(5), -t12 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
