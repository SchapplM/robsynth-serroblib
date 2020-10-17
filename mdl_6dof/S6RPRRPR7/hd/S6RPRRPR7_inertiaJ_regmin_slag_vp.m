% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:01:48
% EndTime: 2019-05-05 23:01:49
% DurationCPUTime: 0.58s
% Computational Cost: add. (696->79), mult. (1198->131), div. (0->0), fcn. (1483->8), ass. (0->62)
t47 = sin(qJ(4));
t48 = sin(qJ(3));
t50 = cos(qJ(4));
t51 = cos(qJ(3));
t30 = t47 * t51 + t50 * t48;
t31 = -t47 * t48 + t50 * t51;
t44 = sin(pkin(10));
t45 = cos(pkin(10));
t16 = t44 * t30 - t45 * t31;
t55 = t45 * t30 + t44 * t31;
t82 = (t16 * t45 - t44 * t55) * pkin(4);
t15 = t16 ^ 2;
t46 = sin(qJ(6));
t81 = t16 * t46;
t49 = cos(qJ(6));
t80 = t16 * t49;
t75 = t55 ^ 2;
t79 = t75 + t15;
t52 = -pkin(1) - pkin(7);
t39 = t51 * t52;
t59 = -t51 * pkin(8) + t39;
t60 = (-pkin(8) + t52) * t48;
t22 = -t47 * t59 - t50 * t60;
t10 = -t30 * qJ(5) - t22;
t21 = -t47 * t60 + t50 * t59;
t54 = -t31 * qJ(5) + t21;
t4 = t44 * t10 - t45 * t54;
t6 = t45 * t10 + t44 * t54;
t78 = t4 * t16 + t55 * t6;
t41 = t50 * pkin(3);
t38 = t41 + pkin(4);
t66 = t47 * pkin(3);
t26 = t45 * t38 - t44 * t66;
t27 = t44 * t38 + t45 * t66;
t77 = t16 * t26 - t27 * t55;
t12 = t46 * t55;
t13 = t49 * t55;
t37 = t48 * pkin(3) + qJ(2);
t71 = 0.2e1 * t37;
t70 = 0.2e1 * t46;
t69 = -0.2e1 * t49;
t68 = 0.2e1 * qJ(2);
t67 = t4 * t49;
t63 = t46 * t49;
t24 = -pkin(5) - t26;
t36 = -t45 * pkin(4) - pkin(5);
t61 = t24 + t36;
t23 = t30 * pkin(4) + t37;
t25 = pkin(9) + t27;
t57 = -t16 * t24 - t25 * t55;
t35 = t44 * pkin(4) + pkin(9);
t56 = -t16 * t36 - t35 * t55;
t43 = t49 ^ 2;
t42 = t46 ^ 2;
t34 = 0.2e1 * t63;
t11 = t46 * t80;
t8 = (-t42 + t43) * t16;
t7 = pkin(5) * t55 + pkin(9) * t16 + t23;
t3 = t4 * t46;
t2 = t46 * t7 + t49 * t6;
t1 = -t46 * t6 + t49 * t7;
t5 = [1, 0, 0, -2 * pkin(1), t68 (pkin(1) ^ 2) + qJ(2) ^ 2, t51 ^ 2, -0.2e1 * t51 * t48, 0, 0, 0, t48 * t68, t51 * t68, t31 ^ 2, -0.2e1 * t31 * t30, 0, 0, 0, t30 * t71, t31 * t71, -0.2e1 * t78, t23 ^ 2 + t4 ^ 2 + t6 ^ 2, t43 * t15, -0.2e1 * t15 * t63, -0.2e1 * t55 * t80, 0.2e1 * t55 * t81, t75, 0.2e1 * t1 * t55 - 0.2e1 * t4 * t81, -0.2e1 * t2 * t55 - 0.2e1 * t4 * t80; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t78, 0, 0, 0, 0, 0, -t79 * t46, -t79 * t49; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48, 0, t39, -t48 * t52, 0, 0, t31, -t30, 0, t21, t22, t77, -t4 * t26 + t6 * t27, -t11, -t8, t12, t13, 0, t57 * t46 - t67, t57 * t49 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48, 0, 0, 0, 0, 0, t31, -t30, 0, -t77, 0, 0, 0, 0, 0, -t80, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t66, 0, t26 ^ 2 + t27 ^ 2, t42, t34, 0, 0, 0, t24 * t69, t24 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, t21, t22, t82 (-t4 * t45 + t44 * t6) * pkin(4), -t11, -t8, t12, t13, 0, t56 * t46 - t67, t56 * t49 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, -t82, 0, 0, 0, 0, 0, -t80, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t41, -t66, 0 (t26 * t45 + t27 * t44) * pkin(4), t42, t34, 0, 0, 0, -t61 * t49, t61 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t44 ^ 2 + t45 ^ 2) * pkin(4) ^ 2, t42, t34, 0, 0, 0, t36 * t69, t36 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, t13, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t81, t55, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t49, 0, -t46 * t25, -t49 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t49, 0, -t46 * t35, -t49 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
