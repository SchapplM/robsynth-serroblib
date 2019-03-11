% Calculate minimal parameter regressor of joint inertia matrix for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t37 = sin(qJ(5));
t56 = 0.2e1 * t37;
t39 = cos(qJ(6));
t55 = pkin(5) * t39;
t28 = sin(pkin(13));
t22 = t28 * pkin(3) + pkin(9);
t36 = sin(qJ(6));
t54 = t22 * t36;
t30 = sin(pkin(7));
t38 = sin(qJ(3));
t53 = t30 * t38;
t41 = cos(qJ(3));
t52 = t30 * t41;
t33 = cos(pkin(12));
t34 = cos(pkin(7));
t51 = t33 * t34;
t50 = t36 * t37;
t49 = t36 * t39;
t40 = cos(qJ(5));
t48 = t36 * t40;
t47 = t39 * t37;
t46 = t39 * t40;
t45 = t40 * t22;
t44 = t40 * t56;
t32 = cos(pkin(13));
t23 = -t32 * pkin(3) - pkin(4);
t29 = sin(pkin(12));
t31 = sin(pkin(6));
t35 = cos(pkin(6));
t43 = t35 * t52 + (-t29 * t38 + t41 * t51) * t31;
t27 = t39 ^ 2;
t26 = t37 ^ 2;
t25 = t36 ^ 2;
t20 = -t40 * pkin(5) - t37 * pkin(10) + t23;
t19 = -t31 * t33 * t30 + t35 * t34;
t18 = (t28 * t41 + t32 * t38) * t30;
t16 = t28 * t53 - t32 * t52;
t15 = t40 * t18 + t37 * t34;
t14 = t37 * t18 - t40 * t34;
t13 = t36 * t20 + t39 * t45;
t12 = t39 * t20 - t36 * t45;
t11 = t35 * t53 + (t29 * t41 + t38 * t51) * t31;
t9 = t39 * t15 + t36 * t16;
t8 = -t36 * t15 + t39 * t16;
t7 = t32 * t11 + t28 * t43;
t5 = t28 * t11 - t32 * t43;
t4 = t19 * t37 + t40 * t7;
t3 = -t19 * t40 + t37 * t7;
t2 = t36 * t5 + t39 * t4;
t1 = -t36 * t4 + t39 * t5;
t6 = [1, t35 ^ 2 + (t29 ^ 2 + t33 ^ 2) * t31 ^ 2, 0, 0, 0, t19 ^ 2 + t5 ^ 2 + t7 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t35, 0, 0, 0, t5 * t16 + t7 * t18 + t19 * t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, t16 ^ 2 + t18 ^ 2 + t34 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t43, -t11 (t28 * t7 - t32 * t5) * pkin(3), 0, 0, 0, 0, 0, -t5 * t40, t5 * t37, 0, 0, 0, 0, 0, -t1 * t40 + t3 * t50, t2 * t40 + t3 * t47; 0, 0, 0, t52, -t53 (-t16 * t32 + t18 * t28) * pkin(3), 0, 0, 0, 0, 0, -t16 * t40, t16 * t37, 0, 0, 0, 0, 0, t14 * t50 - t8 * t40, t14 * t47 + t9 * t40; 0, 0, 1, 0, 0 (t28 ^ 2 + t32 ^ 2) * pkin(3) ^ 2, t26, t44, 0, 0, 0, -0.2e1 * t23 * t40, t23 * t56, t27 * t26, -0.2e1 * t26 * t49, -0.2e1 * t37 * t46, t36 * t44, t40 ^ 2, -0.2e1 * t12 * t40 + 0.2e1 * t26 * t54, 0.2e1 * t26 * t22 * t39 + 0.2e1 * t13 * t40; 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, -t3 * t39, t3 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t14 * t39, t14 * t36; 0, 0, 0, 0, 0, 0, 0, 0, t37, t40, 0, -t37 * t22, -t45, t36 * t47 (-t25 + t27) * t37, -t48, -t46, 0, -t22 * t47 + (-pkin(5) * t37 + pkin(10) * t40) * t36, pkin(10) * t46 + (t54 - t55) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t37, 0, 0, 0, 0, 0, t46, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t25, 0.2e1 * t49, 0, 0, 0, 0.2e1 * t55, -0.2e1 * pkin(5) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t50, -t40, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t39, 0, -t36 * pkin(10), -t39 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t6;
