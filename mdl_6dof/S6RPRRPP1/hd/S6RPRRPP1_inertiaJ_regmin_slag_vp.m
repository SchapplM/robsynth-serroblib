% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t54 = cos(qJ(4));
t79 = pkin(3) * t54;
t51 = cos(pkin(9));
t42 = -t51 * pkin(1) - pkin(2);
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t29 = -t55 * pkin(3) - t53 * pkin(8) + t42;
t27 = t54 * t29;
t68 = qJ(5) * t53;
t49 = sin(pkin(9));
t39 = t49 * pkin(1) + pkin(7);
t52 = sin(qJ(4));
t77 = t39 * t52;
t10 = -t54 * t68 + t27 + (-pkin(4) - t77) * t55;
t71 = t55 * t39;
t65 = t54 * t71;
t13 = t65 + (t29 - t68) * t52;
t48 = sin(pkin(10));
t50 = cos(pkin(10));
t4 = t48 * t10 + t50 * t13;
t76 = t52 * t53;
t75 = t52 * t54;
t74 = t52 * t55;
t73 = t54 * t53;
t72 = t54 * t55;
t70 = t55 * t53;
t69 = -qJ(5) - pkin(8);
t34 = t53 * t39;
t28 = pkin(4) * t76 + t34;
t67 = 0.2e1 * t70;
t33 = t69 * t54;
t62 = t69 * t52;
t20 = -t48 * t33 - t50 * t62;
t22 = -t50 * t33 + t48 * t62;
t66 = t20 ^ 2 + t22 ^ 2;
t43 = -t54 * pkin(4) - pkin(3);
t3 = t50 * t10 - t48 * t13;
t31 = t48 * t54 + t50 * t52;
t24 = t31 * t53;
t26 = -t48 * t76 + t50 * t73;
t64 = t24 * t20 + t26 * t22;
t63 = t20 * t26 - t22 * t24;
t30 = t48 * t52 - t50 * t54;
t60 = t24 * t31 - t26 * t30;
t47 = t55 ^ 2;
t59 = t24 ^ 2 + t26 ^ 2 + t47;
t58 = 0.2e1 * t20 * t31 - 0.2e1 * t22 * t30;
t46 = t54 ^ 2;
t45 = t53 ^ 2;
t44 = t52 ^ 2;
t40 = t50 * pkin(4) + pkin(5);
t37 = t48 * pkin(4) + qJ(6);
t18 = t52 * t29 + t65;
t17 = -t52 * t71 + t27;
t15 = t30 * pkin(5) - t31 * qJ(6) + t43;
t5 = t24 * pkin(5) - t26 * qJ(6) + t28;
t2 = t55 * pkin(5) - t3;
t1 = -t55 * qJ(6) + t4;
t6 = [1, 0, 0 (t49 ^ 2 + t51 ^ 2) * pkin(1) ^ 2, t45, t67, 0, 0, 0, -0.2e1 * t42 * t55, 0.2e1 * t42 * t53, t46 * t45, -0.2e1 * t45 * t75, -0.2e1 * t54 * t70, t52 * t67, t47, -0.2e1 * t17 * t55 + 0.2e1 * t45 * t77, 0.2e1 * t45 * t39 * t54 + 0.2e1 * t18 * t55, -0.2e1 * t4 * t24 - 0.2e1 * t3 * t26, t28 ^ 2 + t3 ^ 2 + t4 ^ 2, 0.2e1 * t2 * t55 + 0.2e1 * t5 * t24, -0.2e1 * t1 * t24 + 0.2e1 * t2 * t26, -0.2e1 * t1 * t55 - 0.2e1 * t5 * t26, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t24 + t4 * t26 - t28 * t55, 0, 0, 0, t1 * t26 + t2 * t24 - t5 * t55; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, t53, t55, 0, -t34, -t71, t52 * t73 (-t44 + t46) * t53, -t74, -t72, 0, -t39 * t73 + (-pkin(3) * t53 + pkin(8) * t55) * t52, pkin(8) * t72 + (t77 - t79) * t53, -t3 * t31 - t4 * t30 + t63, -t3 * t20 + t4 * t22 + t28 * t43, t15 * t24 + t20 * t55 + t5 * t30, -t1 * t30 + t2 * t31 + t63, -t15 * t26 - t22 * t55 - t5 * t31, t1 * t22 + t5 * t15 + t2 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t53, 0, 0, 0, 0, 0, t72, -t74, t60, -t55 * t43 + t64, -t55 * t30, t60, t55 * t31, -t55 * t15 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t44, 0.2e1 * t75, 0, 0, 0, 0.2e1 * t79, -0.2e1 * pkin(3) * t52, t58, t43 ^ 2 + t66, 0.2e1 * t15 * t30, t58, -0.2e1 * t15 * t31, t15 ^ 2 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t76, -t55, t17, -t18 (-t24 * t48 - t26 * t50) * pkin(4) (t3 * t50 + t4 * t48) * pkin(4) (-pkin(5) - t40) * t55 + t3, -t37 * t24 - t40 * t26 (-qJ(6) - t37) * t55 + t4, t1 * t37 - t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t73, 0 (-t24 * t50 + t26 * t48) * pkin(4), -t24, 0, t26, -t24 * t40 + t26 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t54, 0, -t52 * pkin(8), -t54 * pkin(8) (-t30 * t48 - t31 * t50) * pkin(4) (-t20 * t50 + t22 * t48) * pkin(4), -t20, -t37 * t30 - t40 * t31, t22, -t20 * t40 + t22 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t48 ^ 2 + t50 ^ 2) * pkin(4) ^ 2, 0.2e1 * t40, 0, 0.2e1 * t37, t37 ^ 2 + t40 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t24, 0, -t26, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, 0, 0, 0, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t30, 0, -t31, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t26, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
