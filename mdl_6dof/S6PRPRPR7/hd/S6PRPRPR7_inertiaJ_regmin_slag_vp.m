% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t32 = sin(qJ(4));
t35 = cos(qJ(4));
t11 = t32 * pkin(4) - t35 * qJ(5) + qJ(3);
t53 = -0.2e1 * t11;
t52 = 0.2e1 * qJ(3);
t51 = 0.2e1 * qJ(5);
t29 = sin(pkin(6));
t33 = sin(qJ(2));
t50 = t29 * t33;
t36 = cos(qJ(2));
t49 = t29 * t36;
t31 = sin(qJ(6));
t18 = t31 * t32;
t48 = t31 * t35;
t34 = cos(qJ(6));
t47 = t34 * t31;
t19 = t34 * t32;
t46 = t35 * t32;
t38 = -pkin(2) - pkin(8);
t45 = t35 * t38;
t25 = t32 ^ 2;
t27 = t35 ^ 2;
t16 = t25 + t27;
t44 = qJ(5) * t32;
t43 = 0.2e1 * t46;
t42 = t32 * t50;
t41 = t35 * t50;
t30 = cos(pkin(6));
t6 = t30 * t32 + t35 * t49;
t7 = t30 * t35 - t32 * t49;
t40 = t7 * t32 - t6 * t35;
t14 = pkin(4) * t35 + t44;
t37 = -pkin(4) - pkin(9);
t39 = -t35 * t37 + t44;
t26 = t34 ^ 2;
t24 = t31 ^ 2;
t23 = t29 ^ 2;
t21 = t32 * t38;
t20 = t34 * t35;
t17 = t23 * t33 ^ 2;
t13 = (pkin(5) - t38) * t35;
t12 = -t32 * pkin(5) + t21;
t10 = t16 * t38;
t9 = t32 * pkin(9) + t11;
t5 = t6 * t31 + t34 * t50;
t4 = -t31 * t50 + t6 * t34;
t3 = t31 * t13 + t34 * t9;
t2 = t34 * t13 - t31 * t9;
t1 = [1, 0, 0, 0, 0, 0, t23 * t36 ^ 2 + t30 ^ 2 + t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 ^ 2 + t7 ^ 2 + t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, t49, -t50, -t49, t50 (pkin(2) * t36 + qJ(3) * t33) * t29, 0, 0, 0, 0, 0, t42, t41, -t40, -t42, -t41, t11 * t50 + t40 * t38, 0, 0, 0, 0, 0, -t7 * t19 + t4 * t35, t7 * t18 - t5 * t35; 0, 1, 0, 0, -0.2e1 * pkin(2), t52, pkin(2) ^ 2 + qJ(3) ^ 2, t27, -0.2e1 * t46, 0, 0, 0, t32 * t52, t35 * t52, -0.2e1 * t10, t32 * t53, t35 * t53, t16 * t38 ^ 2 + t11 ^ 2, t24 * t25, 0.2e1 * t25 * t47, t31 * t43, t34 * t43, t27, -0.2e1 * t12 * t19 + 0.2e1 * t2 * t35, 0.2e1 * t12 * t18 - 0.2e1 * t3 * t35; 0, 0, 0, 0, 0, 0, -t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, t10, 0, 0, 0, 0, 0, -t16 * t34, t16 * t31; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, t6, t7, -t6 * pkin(4) + t7 * qJ(5), 0, 0, 0, 0, 0, t7 * t31, t7 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t32, 0, t45, -t21, -t14, -t45, t21, t14 * t38, t31 * t19 (-t24 + t26) * t32, t20, -t48, 0, t12 * t31 - t34 * t39, t12 * t34 + t31 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t32, 0, -t35, t32, t14, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t51, pkin(4) ^ 2 + qJ(5) ^ 2, t26, -0.2e1 * t47, 0, 0, 0, t31 * t51, t34 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, -t45, 0, 0, 0, 0, 0, t20, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19, t35, t2, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t31, 0, t34 * t37, -t31 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t1;
