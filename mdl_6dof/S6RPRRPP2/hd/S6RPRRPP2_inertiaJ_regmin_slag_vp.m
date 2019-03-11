% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t39 = sin(qJ(4));
t43 = pkin(4) + pkin(5);
t62 = t43 * t39;
t41 = cos(qJ(4));
t59 = t39 * qJ(5);
t74 = t41 * t43 + t59;
t13 = pkin(3) + t74;
t73 = 0.2e1 * t13;
t52 = -t41 * pkin(4) - t59;
t17 = -pkin(3) + t52;
t72 = -0.2e1 * t17;
t40 = sin(qJ(3));
t71 = 0.2e1 * t40;
t70 = pkin(3) * t39;
t69 = pkin(3) * t41;
t68 = pkin(4) * t39;
t67 = t39 * pkin(8);
t38 = cos(pkin(9));
t24 = -t38 * pkin(1) - pkin(2);
t42 = cos(qJ(3));
t12 = -t42 * pkin(3) - t40 * pkin(8) + t24;
t37 = sin(pkin(9));
t23 = t37 * pkin(1) + pkin(7);
t63 = t42 * t23;
t7 = t39 * t12 + t41 * t63;
t66 = t23 * t39;
t65 = t23 * t41;
t25 = t39 * t40;
t64 = t39 * t41;
t26 = t39 * t42;
t28 = t41 * t40;
t29 = t41 * t42;
t61 = -t41 * t12 + t39 * t63;
t33 = t39 ^ 2;
t35 = t41 ^ 2;
t60 = t33 + t35;
t58 = t41 * qJ(5);
t57 = t41 * qJ(6);
t56 = t42 * qJ(5);
t55 = t42 * t71;
t54 = -0.2e1 * t56 + t7;
t32 = t42 * pkin(4);
t4 = t32 + t61;
t14 = t60 * t40;
t3 = -t56 + t7;
t53 = t3 * t41 + t4 * t39;
t51 = t58 - t68;
t18 = (pkin(8) - qJ(6)) * t39;
t31 = t41 * pkin(8);
t19 = t31 - t57;
t50 = t18 * t39 + t19 * t41;
t49 = t40 * t57 - t4;
t46 = qJ(5) ^ 2;
t45 = 0.2e1 * qJ(5);
t36 = t42 ^ 2;
t34 = t40 ^ 2;
t27 = t35 * t34;
t22 = pkin(8) * t26;
t21 = t40 * t58;
t20 = qJ(6) * t25;
t11 = t33 * t34 + t27 + t36;
t8 = -t21 + (t23 + t68) * t40;
t5 = t21 + (-t23 - t62) * t40;
t2 = t20 + t3;
t1 = t42 * pkin(5) - t49;
t6 = [1, 0, 0 (t37 ^ 2 + t38 ^ 2) * pkin(1) ^ 2, t34, t55, 0, 0, 0, -0.2e1 * t24 * t42, t24 * t71, t27, -0.2e1 * t34 * t64, -0.2e1 * t40 * t29, t39 * t55, t36, 0.2e1 * t34 * t66 + 0.2e1 * t42 * t61, 0.2e1 * t34 * t65 + 0.2e1 * t7 * t42, 0.2e1 * t8 * t25 + 0.2e1 * t4 * t42 (-t3 * t39 + t4 * t41) * t71, -0.2e1 * t8 * t28 - 0.2e1 * t3 * t42, t3 ^ 2 + t4 ^ 2 + t8 ^ 2, 0.2e1 * t1 * t42 - 0.2e1 * t5 * t25, -0.2e1 * t2 * t42 + 0.2e1 * t5 * t28 (-t1 * t41 + t2 * t39) * t71, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t40 - t8 * t42, 0, 0, 0, t5 * t42 + (t1 * t39 + t2 * t41) * t40; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, t40, t42, 0, -t40 * t23, -t63, t39 * t28 (-t33 + t35) * t40, -t26, -t29, 0, t22 + (-t65 - t70) * t40, pkin(8) * t29 + (t66 - t69) * t40, t17 * t25 - t8 * t41 + t22, t53, -t8 * t39 + (-pkin(8) * t42 - t17 * t40) * t41, t53 * pkin(8) + t8 * t17, -t13 * t25 + t18 * t42 + t5 * t41, t13 * t28 - t19 * t42 + t5 * t39 (-t18 * t40 - t2) * t41 + (t19 * t40 - t1) * t39, t1 * t18 + t5 * t13 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t40, 0, 0, 0, 0, 0, t29, -t26, t29, t14, t26, pkin(8) * t14 - t42 * t17, t29, t26, -t14, t42 * t13 + t50 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t33, 0.2e1 * t64, 0, 0, 0, 0.2e1 * t69, -0.2e1 * t70, t41 * t72, 0.2e1 * t60 * pkin(8), t39 * t72, t60 * pkin(8) ^ 2 + t17 ^ 2, t41 * t73, t39 * t73, -0.2e1 * t50, t13 ^ 2 + t18 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t25, -t42, -t61, -t7, -0.2e1 * t32 - t61, t52 * t40, t54, -t4 * pkin(4) + t3 * qJ(5) (-pkin(5) - t43) * t42 + t49, t20 + t54, t74 * t40, t2 * qJ(5) - t1 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t28, -t25, 0, t28, -pkin(4) * t25 + t21, -t25, t28, 0, -t40 * t62 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t41, 0, -t67, -t31, -t67, t51, t31, t51 * pkin(8), -t18, t19, -t58 + t62, t19 * qJ(5) - t18 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, t45, pkin(4) ^ 2 + t46, 0.2e1 * t43, t45, 0, t43 ^ 2 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t28, 0, t4, t42, 0, -t28, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t67, 0, 0, -t39, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), -1, 0, 0, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t28, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t39, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
