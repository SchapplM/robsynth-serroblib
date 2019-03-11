% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t58 = sin(pkin(9));
t60 = cos(pkin(9));
t62 = sin(qJ(2));
t85 = cos(qJ(2));
t41 = t58 * t62 - t60 * t85;
t43 = t58 * t85 + t60 * t62;
t54 = -t85 * pkin(2) - pkin(1);
t26 = t41 * pkin(3) - t43 * pkin(8) + t54;
t61 = sin(qJ(4));
t79 = qJ(5) * t43;
t46 = (-qJ(3) - pkin(7)) * t62;
t75 = t85 * pkin(7);
t47 = t85 * qJ(3) + t75;
t30 = t58 * t46 + t60 * t47;
t63 = cos(qJ(4));
t81 = t63 * t30;
t10 = t81 + (t26 - t79) * t61;
t57 = sin(pkin(10));
t59 = cos(pkin(10));
t12 = t63 * t26 - t61 * t30;
t8 = t41 * pkin(4) - t63 * t79 + t12;
t4 = t59 * t10 + t57 * t8;
t84 = t61 * t41;
t83 = t61 * t43;
t82 = t61 * t63;
t80 = t63 * t43;
t78 = 0.2e1 * t85;
t1 = t41 * qJ(6) + t4;
t74 = t58 * pkin(2) + pkin(8);
t69 = qJ(5) + t74;
t36 = t69 * t63;
t67 = t69 * t61;
t23 = t57 * t36 + t59 * t67;
t25 = t59 * t36 - t57 * t67;
t77 = t23 ^ 2 + t25 ^ 2;
t39 = t57 * t61 - t59 * t63;
t42 = t57 * t63 + t59 * t61;
t76 = t39 ^ 2 + t42 ^ 2;
t53 = -t60 * pkin(2) - pkin(3);
t3 = -t57 * t10 + t59 * t8;
t19 = t42 * t43;
t20 = -t57 * t83 + t59 * t80;
t73 = -t25 * t19 + t23 * t20;
t72 = -t42 * t19 + t39 * t20;
t71 = t23 * t39 + t25 * t42;
t28 = -t60 * t46 + t58 * t47;
t17 = pkin(4) * t83 + t28;
t45 = -t63 * pkin(4) + t53;
t68 = 0.2e1 * t23 * t42 - 0.2e1 * t25 * t39;
t66 = -t74 * t41 + t53 * t43;
t56 = t63 ^ 2;
t55 = t61 ^ 2;
t51 = t59 * pkin(4) + pkin(5);
t48 = t57 * pkin(4) + qJ(6);
t38 = t43 ^ 2;
t34 = t63 * t41;
t18 = t39 * pkin(5) - t42 * qJ(6) + t45;
t13 = t61 * t26 + t81;
t5 = t19 * pkin(5) - t20 * qJ(6) + t17;
t2 = -t41 * pkin(5) - t3;
t6 = [1, 0, 0, t62 ^ 2, t62 * t78, 0, 0, 0, pkin(1) * t78, -0.2e1 * pkin(1) * t62, 0.2e1 * t28 * t43 - 0.2e1 * t30 * t41, t28 ^ 2 + t30 ^ 2 + t54 ^ 2, t56 * t38, -0.2e1 * t38 * t82, 0.2e1 * t41 * t80, -0.2e1 * t41 * t83, t41 ^ 2, 0.2e1 * t12 * t41 + 0.2e1 * t28 * t83, -0.2e1 * t13 * t41 + 0.2e1 * t28 * t80, -0.2e1 * t4 * t19 - 0.2e1 * t3 * t20, t17 ^ 2 + t3 ^ 2 + t4 ^ 2, 0.2e1 * t5 * t19 - 0.2e1 * t2 * t41, -0.2e1 * t1 * t19 + 0.2e1 * t2 * t20, 0.2e1 * t1 * t41 - 0.2e1 * t5 * t20, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t62, t85, 0, -t62 * pkin(7), -t75 (-t41 * t58 - t43 * t60) * pkin(2) (-t28 * t60 + t30 * t58) * pkin(2), t61 * t80 (-t55 + t56) * t43, t84, t34, 0, -t28 * t63 + t66 * t61, t28 * t61 + t66 * t63, -t3 * t42 - t4 * t39 + t73, t17 * t45 - t3 * t23 + t4 * t25, t18 * t19 - t23 * t41 + t5 * t39, -t1 * t39 + t2 * t42 + t73, -t18 * t20 + t25 * t41 - t5 * t42, t1 * t25 + t5 * t18 + t2 * t23; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t58 ^ 2 + t60 ^ 2) * pkin(2) ^ 2, t55, 0.2e1 * t82, 0, 0, 0, -0.2e1 * t53 * t63, 0.2e1 * t53 * t61, t68, t45 ^ 2 + t77, 0.2e1 * t18 * t39, t68, -0.2e1 * t18 * t42, t18 ^ 2 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, t34, -t84, t72, -t3 * t39 + t4 * t42, -t39 * t41, t72, t42 * t41, t1 * t42 + t2 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, 0, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t83, t41, t12, -t13 (-t19 * t57 - t20 * t59) * pkin(4) (t3 * t59 + t4 * t57) * pkin(4) (pkin(5) + t51) * t41 + t3, -t48 * t19 - t51 * t20, t48 * t41 + t1, t1 * t48 - t2 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t63, 0, -t61 * t74, -t63 * t74 (-t39 * t57 - t42 * t59) * pkin(4) (-t23 * t59 + t25 * t57) * pkin(4), -t23, -t48 * t39 - t51 * t42, t25, -t23 * t51 + t25 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t61, 0 (-t39 * t59 + t42 * t57) * pkin(4), -t39, 0, t42, -t39 * t51 + t42 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t57 ^ 2 + t59 ^ 2) * pkin(4) ^ 2, 0.2e1 * t51, 0, 0.2e1 * t48, t48 ^ 2 + t51 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t19, 0, -t20, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t39, 0, -t42, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t20, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
