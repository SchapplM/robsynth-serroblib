% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t53 = cos(qJ(4));
t81 = 0.2e1 * t53;
t80 = 2 * qJ(2);
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t34 = t52 * pkin(3) - t54 * pkin(8) + qJ(2);
t29 = t53 * t34;
t68 = qJ(5) * t54;
t51 = sin(qJ(4));
t55 = -pkin(1) - pkin(7);
t76 = t51 * t55;
t13 = -t53 * t68 + t29 + (pkin(4) - t76) * t52;
t73 = t53 * t55;
t65 = t52 * t73;
t16 = t65 + (t34 - t68) * t51;
t49 = sin(pkin(9));
t50 = cos(pkin(9));
t4 = t49 * t13 + t50 * t16;
t79 = t51 * t52;
t78 = t51 * t53;
t77 = t51 * t54;
t75 = t52 * t55;
t74 = t53 * t52;
t42 = t53 * t54;
t31 = t49 * t53 + t50 * t51;
t26 = t54 * t31;
t72 = t54 * t52;
t71 = t54 * t55;
t70 = -qJ(5) - pkin(8);
t46 = t52 ^ 2;
t48 = t54 ^ 2;
t69 = -t46 - t48;
t67 = -0.2e1 * t72;
t35 = t70 * t53;
t62 = t70 * t51;
t18 = -t49 * t35 - t50 * t62;
t20 = -t50 * t35 + t49 * t62;
t66 = t18 ^ 2 + t20 ^ 2;
t1 = t52 * qJ(6) + t4;
t43 = -t53 * pkin(4) - pkin(3);
t24 = t31 * t52;
t30 = t49 * t51 - t50 * t53;
t27 = t30 * t52;
t64 = t24 * t18 - t27 * t20;
t28 = t50 * t42 - t49 * t77;
t63 = t18 * t28 - t20 * t26;
t3 = t50 * t13 - t49 * t16;
t61 = t24 * t28 + t27 * t26;
t60 = t24 * t31 + t27 * t30;
t32 = pkin(4) * t77 - t71;
t59 = t24 ^ 2 + t27 ^ 2 + t48;
t58 = -pkin(3) * t54 - pkin(8) * t52;
t57 = 0.2e1 * t18 * t31 - 0.2e1 * t20 * t30;
t47 = t53 ^ 2;
t45 = t51 ^ 2;
t40 = t50 * pkin(4) + pkin(5);
t38 = t49 * pkin(4) + qJ(6);
t22 = t51 * t34 + t65;
t21 = -t51 * t75 + t29;
t12 = t30 * pkin(5) - t31 * qJ(6) + t43;
t5 = t26 * pkin(5) - t28 * qJ(6) + t32;
t2 = -t52 * pkin(5) - t3;
t6 = [1, 0, 0, -2 * pkin(1), t80, pkin(1) ^ 2 + qJ(2) ^ 2, t48, t67, 0, 0, 0, t52 * t80, t54 * t80, t47 * t48, -0.2e1 * t48 * t78, t72 * t81, t51 * t67, t46, 0.2e1 * t21 * t52 - 0.2e1 * t48 * t76, -0.2e1 * t22 * t52 - 0.2e1 * t48 * t73, -0.2e1 * t4 * t26 - 0.2e1 * t3 * t28, t3 ^ 2 + t32 ^ 2 + t4 ^ 2, -0.2e1 * t2 * t52 + 0.2e1 * t5 * t26, -0.2e1 * t1 * t26 + 0.2e1 * t2 * t28, 0.2e1 * t1 * t52 - 0.2e1 * t5 * t28, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t51, t69 * t53, t61, -t3 * t24 - t4 * t27 - t32 * t54, -t24 * t52 - t54 * t26, t61, -t27 * t52 + t54 * t28, -t1 * t27 + t2 * t24 - t5 * t54; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, t54, -t52, 0, t71, -t75, t51 * t42 (-t45 + t47) * t54, t79, t74, 0, t58 * t51 + t53 * t71, -t51 * t71 + t58 * t53, -t3 * t31 - t4 * t30 + t63, -t3 * t18 + t4 * t20 + t32 * t43, t12 * t26 - t18 * t52 + t5 * t30, -t1 * t30 + t2 * t31 + t63, -t12 * t28 + t20 * t52 - t5 * t31, t1 * t20 + t5 * t12 + t2 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t52, 0, 0, 0, 0, 0, t42, -t77, t60, -t54 * t43 + t64, -t54 * t30, t60, t26, -t54 * t12 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t45, 0.2e1 * t78, 0, 0, 0, pkin(3) * t81, -0.2e1 * pkin(3) * t51, t57, t43 ^ 2 + t66, 0.2e1 * t12 * t30, t57, -0.2e1 * t12 * t31, t12 ^ 2 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t77, t52, t21, -t22 (-t26 * t49 - t28 * t50) * pkin(4) (t3 * t50 + t4 * t49) * pkin(4) (pkin(5) + t40) * t52 + t3, -t38 * t26 - t40 * t28, t38 * t52 + t1, t1 * t38 - t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t74, 0 (-t24 * t50 - t27 * t49) * pkin(4), -t24, 0, -t27, -t24 * t40 - t27 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t53, 0, -t51 * pkin(8), -t53 * pkin(8) (-t30 * t49 - t31 * t50) * pkin(4) (-t18 * t50 + t20 * t49) * pkin(4), -t18, -t38 * t30 - t40 * t31, t20, -t18 * t40 + t20 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t49 ^ 2 + t50 ^ 2) * pkin(4) ^ 2, 0.2e1 * t40, 0, 0.2e1 * t38, t38 ^ 2 + t40 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t26, 0, -t28, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, 0, 0, 0, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t30, 0, -t31, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t28, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
