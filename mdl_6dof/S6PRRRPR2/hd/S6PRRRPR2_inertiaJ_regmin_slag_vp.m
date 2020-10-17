% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:16:23
% EndTime: 2019-05-05 07:16:26
% DurationCPUTime: 0.84s
% Computational Cost: add. (795->119), mult. (1718->208), div. (0->0), fcn. (2122->12), ass. (0->87)
t69 = sin(pkin(12));
t71 = cos(pkin(12));
t84 = t69 ^ 2 + t71 ^ 2;
t85 = t84 * qJ(5);
t74 = sin(qJ(4));
t75 = sin(qJ(3));
t78 = cos(qJ(3));
t94 = cos(qJ(4));
t47 = t74 * t75 - t94 * t78;
t103 = -0.2e1 * t47;
t83 = t94 * pkin(3);
t62 = -t83 - pkin(4);
t97 = t71 * pkin(5);
t50 = t62 - t97;
t102 = 0.2e1 * t50;
t60 = -pkin(4) - t97;
t101 = 0.2e1 * t60;
t63 = -t78 * pkin(3) - pkin(2);
t100 = 0.2e1 * t63;
t99 = 0.2e1 * t78;
t98 = -pkin(9) - pkin(8);
t96 = t74 * pkin(3);
t95 = pkin(4) - t62;
t72 = cos(pkin(6));
t70 = sin(pkin(6));
t90 = t70 * sin(qJ(2));
t40 = t72 * t78 - t75 * t90;
t41 = t72 * t75 + t78 * t90;
t21 = -t94 * t40 + t74 * t41;
t93 = t21 * t71;
t55 = t98 * t75;
t56 = t98 * t78;
t37 = -t94 * t55 - t74 * t56;
t92 = t37 * t71;
t48 = t74 * t78 + t94 * t75;
t91 = t69 * t48;
t89 = t70 * cos(qJ(2));
t88 = t71 * t48;
t29 = t47 * pkin(4) - t48 * qJ(5) + t63;
t38 = t74 * t55 - t94 * t56;
t11 = t69 * t29 + t71 * t38;
t87 = t50 + t60;
t59 = qJ(5) + t96;
t86 = t84 * t59;
t10 = t71 * t29 - t69 * t38;
t82 = -pkin(4) * t48 - qJ(5) * t47;
t3 = -t10 * t69 + t11 * t71;
t22 = t74 * t40 + t94 * t41;
t16 = -t69 * t22 - t71 * t89;
t17 = t71 * t22 - t69 * t89;
t4 = -t16 * t69 + t17 * t71;
t81 = -t47 * t59 + t48 * t62;
t73 = sin(qJ(6));
t77 = cos(qJ(6));
t46 = t77 * t69 + t73 * t71;
t45 = t73 * t69 - t77 * t71;
t66 = t71 * pkin(10);
t52 = t71 * qJ(5) + t66;
t51 = (-pkin(10) - qJ(5)) * t69;
t44 = t46 ^ 2;
t43 = t71 * t59 + t66;
t42 = (-pkin(10) - t59) * t69;
t36 = t46 * t47;
t35 = t45 * t47;
t34 = t73 * t51 + t77 * t52;
t33 = t77 * t51 - t73 * t52;
t32 = -0.2e1 * t46 * t45;
t31 = t37 * t69;
t28 = t73 * t42 + t77 * t43;
t27 = t77 * t42 - t73 * t43;
t24 = t45 * t48;
t23 = t46 * t48;
t20 = pkin(5) * t91 + t37;
t19 = t21 * t69;
t18 = t24 * t46;
t15 = t21 * t46;
t14 = t21 * t45;
t13 = t20 * t46;
t12 = t20 * t45;
t9 = -pkin(10) * t91 + t11;
t8 = t47 * pkin(5) - pkin(10) * t88 + t10;
t7 = -t46 * t23 + t24 * t45;
t6 = t73 * t16 + t77 * t17;
t5 = t77 * t16 - t73 * t17;
t2 = t73 * t8 + t77 * t9;
t1 = -t73 * t9 + t77 * t8;
t25 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 ^ 2 + t17 ^ 2 + t21 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t89, -t90, 0, 0, 0, 0, 0, t78 * t89, -t75 * t89, 0, 0, 0, 0, 0, -t47 * t89, -t48 * t89, t16 * t47 + t21 * t91, -t17 * t47 + t21 * t88 (-t16 * t71 - t17 * t69) * t48, t16 * t10 + t17 * t11 + t21 * t37, 0, 0, 0, 0, 0, t21 * t23 + t5 * t47, -t21 * t24 - t6 * t47; 0, 1, 0, 0, t75 ^ 2, t75 * t99, 0, 0, 0, pkin(2) * t99, -0.2e1 * pkin(2) * t75, t48 ^ 2, t48 * t103, 0, 0, 0, t47 * t100, t48 * t100, 0.2e1 * t10 * t47 + 0.2e1 * t37 * t91, -0.2e1 * t11 * t47 + 0.2e1 * t37 * t88, 0.2e1 * (-t10 * t71 - t11 * t69) * t48, t10 ^ 2 + t11 ^ 2 + t37 ^ 2, t24 ^ 2, 0.2e1 * t24 * t23, t24 * t103, t23 * t103, t47 ^ 2, 0.2e1 * t1 * t47 + 0.2e1 * t20 * t23, -0.2e1 * t2 * t47 - 0.2e1 * t20 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t41, 0, 0, 0, 0, 0, -t21, -t22, -t93, t19, t4, t21 * t62 + t4 * t59, 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, t75, t78, 0, -t75 * pkin(8), -t78 * pkin(8), 0, 0, t48, -t47, 0, -t37, -t38, t81 * t69 - t92, t81 * t71 + t31, t3, t3 * t59 + t37 * t62, -t18, t7, t36, -t35, 0, t50 * t23 + t27 * t47 + t12, -t50 * t24 - t28 * t47 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t83, -0.2e1 * t96, -0.2e1 * t62 * t71, 0.2e1 * t62 * t69, 0.2e1 * t86, t84 * t59 ^ 2 + t62 ^ 2, t44, t32, 0, 0, 0, t45 * t102, t46 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, -t93, t19, t4, -t21 * pkin(4) + t4 * qJ(5), 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t47, 0, -t37, -t38, t82 * t69 - t92, t82 * t71 + t31, t3, -t37 * pkin(4) + t3 * qJ(5), -t18, t7, t36, -t35, 0, t60 * t23 + t33 * t47 + t12, -t60 * t24 - t34 * t47 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t83, -t96, t95 * t71, -t95 * t69, t85 + t86, -t62 * pkin(4) + t59 * t85, t44, t32, 0, 0, 0, t87 * t45, t87 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t71, -0.2e1 * pkin(4) * t69, 0.2e1 * t85, t84 * qJ(5) ^ 2 + pkin(4) ^ 2, t44, t32, 0, 0, 0, t45 * t101, t46 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t88, 0, t37, 0, 0, 0, 0, 0, t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t69, 0, t62, 0, 0, 0, 0, 0, t45, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t69, 0, -pkin(4), 0, 0, 0, 0, 0, t45, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, t47, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, 0, t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, 0, t33, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t25;
