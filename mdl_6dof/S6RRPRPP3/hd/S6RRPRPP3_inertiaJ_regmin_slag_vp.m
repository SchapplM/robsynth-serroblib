% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t46 = sin(pkin(9));
t47 = cos(pkin(9));
t50 = sin(qJ(4));
t74 = cos(qJ(4));
t82 = -t50 * t46 + t74 * t47;
t51 = sin(qJ(2));
t25 = t82 * t51;
t81 = -0.2e1 * t25;
t30 = t74 * t46 + t50 * t47;
t80 = -0.2e1 * t30;
t39 = -t47 * pkin(3) - pkin(2);
t79 = 0.2e1 * t39;
t52 = cos(qJ(2));
t78 = 0.2e1 * t52;
t77 = pkin(7) * t46;
t24 = t30 * t51;
t76 = t24 * pkin(5);
t41 = t51 * pkin(7);
t75 = t52 * pkin(7);
t68 = pkin(8) + qJ(3);
t33 = t68 * t46;
t34 = t68 * t47;
t19 = t74 * t33 + t50 * t34;
t73 = t19 * t52;
t20 = -t50 * t33 + t74 * t34;
t72 = t20 * t52;
t71 = t46 * t51;
t70 = t47 * t51;
t48 = pkin(4) + qJ(6);
t32 = -t52 * pkin(2) - t51 * qJ(3) - pkin(1);
t27 = t47 * t32;
t15 = -pkin(8) * t70 + t27 + (-pkin(3) - t77) * t52;
t23 = t46 * t32 + t47 * t75;
t18 = -pkin(8) * t71 + t23;
t67 = -t74 * t15 + t50 * t18;
t7 = t50 * t15 + t74 * t18;
t31 = pkin(3) * t71 + t41;
t66 = t46 ^ 2 + t47 ^ 2;
t65 = qJ(5) * t24;
t64 = qJ(5) * t82;
t63 = t52 * qJ(5);
t42 = t52 * pkin(4);
t4 = t42 + t67;
t62 = -0.2e1 * t63 + t7;
t3 = t63 - t7;
t60 = -t25 * qJ(5) + t31;
t59 = -pkin(2) * t51 + qJ(3) * t52;
t22 = -t46 * t75 + t27;
t58 = -t22 * t46 + t23 * t47;
t57 = -t25 * pkin(5) - t4;
t56 = -t30 * qJ(5) + t39;
t54 = qJ(5) ^ 2;
t53 = 0.2e1 * qJ(5);
t45 = t51 ^ 2;
t14 = -pkin(4) * t82 + t56;
t11 = pkin(5) * t82 + t20;
t10 = t30 * pkin(5) + t19;
t9 = -t48 * t82 + t56;
t8 = t24 * pkin(4) + t60;
t5 = t48 * t24 + t60;
t2 = -t3 - t76;
t1 = t52 * qJ(6) - t57;
t6 = [1, 0, 0, t45, t51 * t78, 0, 0, 0, pkin(1) * t78, -0.2e1 * pkin(1) * t51, -0.2e1 * t22 * t52 + 0.2e1 * t45 * t77, 0.2e1 * t45 * pkin(7) * t47 + 0.2e1 * t23 * t52, 0.2e1 * (-t22 * t47 - t23 * t46) * t51, t45 * pkin(7) ^ 2 + t22 ^ 2 + t23 ^ 2, t25 ^ 2, t24 * t81, t52 * t81, t24 * t78, t52 ^ 2, 0.2e1 * t31 * t24 + 0.2e1 * t52 * t67, 0.2e1 * t31 * t25 + 0.2e1 * t7 * t52, 0.2e1 * t3 * t24 + 0.2e1 * t4 * t25, -0.2e1 * t8 * t24 - 0.2e1 * t4 * t52, -0.2e1 * t8 * t25 + 0.2e1 * t3 * t52, t3 ^ 2 + t4 ^ 2 + t8 ^ 2, 0.2e1 * t1 * t25 - 0.2e1 * t2 * t24, -0.2e1 * t2 * t52 - 0.2e1 * t5 * t25, 0.2e1 * t1 * t52 + 0.2e1 * t5 * t24, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t51, t52, 0, -t41, -t75, -pkin(7) * t70 + t59 * t46, pkin(7) * t71 + t59 * t47, t58, -pkin(2) * t41 + t58 * qJ(3), t25 * t30, -t30 * t24 + t25 * t82, -t30 * t52, -t82 * t52, 0, t39 * t24 - t31 * t82 + t73, t39 * t25 + t31 * t30 + t72, t19 * t25 - t20 * t24 - t3 * t82 + t4 * t30, -t14 * t24 + t8 * t82 - t73, -t14 * t25 - t8 * t30 - t72, t8 * t14 + t4 * t19 - t3 * t20, t1 * t30 + t10 * t25 - t11 * t24 + t2 * t82, -t11 * t52 - t9 * t25 - t5 * t30, t10 * t52 + t9 * t24 - t5 * t82, t1 * t10 + t2 * t11 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t47, -0.2e1 * pkin(2) * t46, 0.2e1 * t66 * qJ(3), t66 * qJ(3) ^ 2 + pkin(2) ^ 2, t30 ^ 2, -t82 * t80, 0, 0, 0, -t82 * t79, t30 * t79, 0.2e1 * t19 * t30 + 0.2e1 * t20 * t82, 0.2e1 * t14 * t82, t14 * t80, t14 ^ 2 + t19 ^ 2 + t20 ^ 2, 0.2e1 * t10 * t30 + 0.2e1 * t11 * t82, t9 * t80, -0.2e1 * t9 * t82, t10 ^ 2 + t11 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t70, 0, t41, 0, 0, 0, 0, 0, t24, t25, 0, -t24, -t25, t8, 0, -t25, t24, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t46, 0, -pkin(2), 0, 0, 0, 0, 0, -t82, t30, 0, t82, -t30, t14, 0, -t30, -t82, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, -t52, -t67, -t7, -pkin(4) * t25 - t65, 0.2e1 * t42 + t67, t62, -t4 * pkin(4) - t3 * qJ(5), -t48 * t25 - t65, t62 - t76 (-qJ(6) - t48) * t52 + t57, t2 * qJ(5) - t1 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t82, 0, -t19, -t20, -pkin(4) * t30 + t64, t19, t20, -t19 * pkin(4) + t20 * qJ(5), -t48 * t30 + t64, t11, -t10, t11 * qJ(5) - t10 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t53, pkin(4) ^ 2 + t54, 0, t53, 0.2e1 * t48, t48 ^ 2 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t52, 0, t4, t25, 0, t52, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, t19, t30, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, -1, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t52, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
