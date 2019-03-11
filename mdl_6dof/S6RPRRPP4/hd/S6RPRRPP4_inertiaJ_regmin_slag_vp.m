% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPP4
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
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t60 = cos(pkin(9));
t51 = -t60 * pkin(2) - pkin(1);
t86 = 0.2e1 * t51;
t58 = sin(pkin(9));
t62 = sin(qJ(3));
t85 = cos(qJ(3));
t39 = t62 * t58 - t85 * t60;
t41 = t85 * t58 + t62 * t60;
t22 = t39 * pkin(3) - t41 * pkin(8) + t51;
t61 = sin(qJ(4));
t76 = qJ(5) * t41;
t79 = pkin(7) + qJ(2);
t43 = t79 * t58;
t44 = t79 * t60;
t24 = -t62 * t43 + t85 * t44;
t63 = cos(qJ(4));
t81 = t63 * t24;
t10 = t81 + (t22 - t76) * t61;
t57 = sin(pkin(10));
t59 = cos(pkin(10));
t11 = t63 * t22 - t61 * t24;
t8 = t39 * pkin(4) - t63 * t76 + t11;
t4 = t59 * t10 + t57 * t8;
t84 = t61 * t39;
t83 = t61 * t41;
t82 = t61 * t63;
t80 = t63 * t41;
t78 = -qJ(5) - pkin(8);
t77 = t58 ^ 2 + t60 ^ 2;
t75 = -0.2e1 * t41 * t39;
t1 = t39 * qJ(6) + t4;
t45 = t78 * t63;
t72 = t78 * t61;
t27 = -t57 * t45 - t59 * t72;
t29 = -t59 * t45 + t57 * t72;
t74 = t27 ^ 2 + t29 ^ 2;
t37 = t57 * t61 - t59 * t63;
t40 = t57 * t63 + t59 * t61;
t73 = t37 ^ 2 + t40 ^ 2;
t52 = -t63 * pkin(4) - pkin(3);
t3 = -t57 * t10 + t59 * t8;
t18 = t40 * t41;
t19 = -t57 * t83 + t59 * t80;
t71 = -t29 * t18 + t27 * t19;
t70 = -t40 * t18 + t37 * t19;
t69 = t37 * t27 + t40 * t29;
t67 = -pkin(3) * t41 - pkin(8) * t39;
t23 = t85 * t43 + t62 * t44;
t66 = 0.2e1 * t27 * t40 - 0.2e1 * t29 * t37;
t15 = pkin(4) * t83 + t23;
t56 = t63 ^ 2;
t55 = t61 ^ 2;
t49 = t59 * pkin(4) + pkin(5);
t46 = t57 * pkin(4) + qJ(6);
t35 = t41 ^ 2;
t32 = t63 * t39;
t21 = t37 * pkin(5) - t40 * qJ(6) + t52;
t12 = t61 * t22 + t81;
t5 = t18 * pkin(5) - t19 * qJ(6) + t15;
t2 = -t39 * pkin(5) - t3;
t6 = [1, 0, 0, 0.2e1 * pkin(1) * t60, -0.2e1 * pkin(1) * t58, 0.2e1 * t77 * qJ(2), t77 * qJ(2) ^ 2 + pkin(1) ^ 2, t35, t75, 0, 0, 0, t39 * t86, t41 * t86, t56 * t35, -0.2e1 * t35 * t82, 0.2e1 * t39 * t80, t61 * t75, t39 ^ 2, 0.2e1 * t11 * t39 + 0.2e1 * t23 * t83, -0.2e1 * t12 * t39 + 0.2e1 * t23 * t80, -0.2e1 * t4 * t18 - 0.2e1 * t3 * t19, t15 ^ 2 + t3 ^ 2 + t4 ^ 2, 0.2e1 * t5 * t18 - 0.2e1 * t2 * t39, -0.2e1 * t1 * t18 + 0.2e1 * t2 * t19, 0.2e1 * t1 * t39 - 0.2e1 * t5 * t19, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -t60, t58, 0, -pkin(1), 0, 0, 0, 0, 0, t39, t41, 0, 0, 0, 0, 0, t32, -t84, t70, -t3 * t37 + t4 * t40, -t37 * t39, t70, t40 * t39, t1 * t40 + t2 * t37; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, 0, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t39, 0, -t23, -t24, t61 * t80 (-t55 + t56) * t41, t84, t32, 0, -t23 * t63 + t67 * t61, t23 * t61 + t67 * t63, -t3 * t40 - t4 * t37 + t71, t15 * t52 - t3 * t27 + t4 * t29, t21 * t18 - t27 * t39 + t5 * t37, -t1 * t37 + t2 * t40 + t71, -t21 * t19 + t29 * t39 - t5 * t40, t1 * t29 + t2 * t27 + t5 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t55, 0.2e1 * t82, 0, 0, 0, 0.2e1 * pkin(3) * t63, -0.2e1 * pkin(3) * t61, t66, t52 ^ 2 + t74, 0.2e1 * t21 * t37, t66, -0.2e1 * t21 * t40, t21 ^ 2 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t83, t39, t11, -t12 (-t18 * t57 - t19 * t59) * pkin(4) (t3 * t59 + t4 * t57) * pkin(4) (pkin(5) + t49) * t39 + t3, -t46 * t18 - t49 * t19, t46 * t39 + t1, t1 * t46 - t2 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t61, 0 (-t37 * t59 + t40 * t57) * pkin(4), -t37, 0, t40, -t37 * t49 + t40 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t63, 0, -t61 * pkin(8), -t63 * pkin(8) (-t37 * t57 - t40 * t59) * pkin(4) (-t27 * t59 + t29 * t57) * pkin(4), -t27, -t46 * t37 - t49 * t40, t29, -t27 * t49 + t29 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t57 ^ 2 + t59 ^ 2) * pkin(4) ^ 2, 0.2e1 * t49, 0, 0.2e1 * t46, t46 ^ 2 + t49 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t18, 0, -t19, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t37, 0, -t40, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t19, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
