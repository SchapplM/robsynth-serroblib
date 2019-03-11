% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t114 = cos(qJ(5));
t76 = sin(pkin(11));
t78 = cos(pkin(11));
t80 = sin(qJ(5));
t56 = -t114 * t78 + t80 * t76;
t77 = sin(pkin(6));
t84 = cos(qJ(2));
t104 = t77 * t84;
t82 = sin(qJ(2));
t105 = t77 * t82;
t79 = cos(pkin(6));
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t50 = t83 * t105 + t79 * t81;
t35 = t78 * t104 + t50 * t76;
t36 = -t76 * t104 + t50 * t78;
t18 = t114 * t35 + t80 * t36;
t126 = -0.2e1 * t18;
t48 = t56 * t81;
t125 = 0.2e1 * t48;
t124 = -0.2e1 * t50;
t57 = t114 * t76 + t80 * t78;
t123 = -0.2e1 * t57;
t70 = -t78 * pkin(4) - pkin(3);
t122 = 0.2e1 * t70;
t121 = 0.2e1 * t83;
t119 = pkin(1) * t84;
t65 = pkin(8) * t105;
t43 = t65 + (-pkin(2) - t119) * t79;
t49 = t81 * t105 - t79 * t83;
t23 = t49 * pkin(3) - t50 * qJ(4) + t43;
t120 = pkin(1) * t82;
t94 = pkin(8) * t104;
t44 = t94 + (pkin(9) + t120) * t79;
t45 = (-pkin(2) * t84 - pkin(9) * t82 - pkin(1)) * t77;
t28 = t83 * t44 + t81 * t45;
t24 = -qJ(4) * t104 + t28;
t12 = t76 * t23 + t78 * t24;
t10 = -pkin(10) * t35 + t12;
t11 = t78 * t23 - t76 * t24;
t8 = t49 * pkin(4) - t36 * pkin(10) + t11;
t4 = t114 * t10 + t80 * t8;
t118 = pkin(9) * t76;
t117 = t49 * pkin(5);
t71 = t81 * pkin(9);
t116 = t83 * pkin(5);
t115 = t83 * pkin(9);
t27 = -t81 * t44 + t83 * t45;
t25 = pkin(3) * t104 - t27;
t113 = t25 * t76;
t112 = t25 * t78;
t100 = pkin(10) + qJ(4);
t60 = t100 * t78;
t89 = t100 * t76;
t38 = t114 * t89 + t80 * t60;
t111 = t38 * t49;
t110 = t38 * t83;
t39 = t114 * t60 - t80 * t89;
t109 = t39 * t49;
t108 = t39 * t83;
t73 = t77 ^ 2;
t107 = t73 * t84;
t106 = t76 * t81;
t103 = t78 * t81;
t102 = t79 * t82;
t59 = -t83 * pkin(3) - t81 * qJ(4) - pkin(2);
t54 = t78 * t59;
t32 = -pkin(10) * t103 + t54 + (-pkin(4) - t118) * t83;
t42 = t78 * t115 + t76 * t59;
t37 = -pkin(10) * t106 + t42;
t17 = t114 * t37 + t80 * t32;
t58 = pkin(4) * t106 + t71;
t99 = t76 ^ 2 + t78 ^ 2;
t98 = qJ(4) * t49;
t97 = qJ(6) * t49;
t96 = t83 * qJ(6);
t95 = 0.2e1 * t104;
t93 = t81 * t104;
t92 = t83 * t104;
t90 = t80 * t10 - t114 * t8;
t16 = t114 * t32 - t80 * t37;
t88 = -pkin(3) * t81 + qJ(4) * t83;
t87 = -t11 * t76 + t12 * t78;
t41 = -t76 * t115 + t54;
t86 = -t41 * t76 + t42 * t78;
t13 = t35 * pkin(4) + t25;
t75 = t81 ^ 2;
t52 = pkin(1) * t102 + t94;
t51 = t79 * t119 - t65;
t47 = t57 * t81;
t31 = t56 * pkin(5) - t57 * qJ(6) + t70;
t26 = t47 * pkin(5) + t48 * qJ(6) + t58;
t19 = t114 * t36 - t80 * t35;
t15 = -t16 + t116;
t14 = -t96 + t17;
t5 = t18 * pkin(5) - t19 * qJ(6) + t13;
t2 = t90 - t117;
t1 = t97 + t4;
t3 = [1, 0, 0, t73 * t82 ^ 2, 0.2e1 * t82 * t107, 0.2e1 * t77 * t102, t79 * t95, t79 ^ 2, 0.2e1 * pkin(1) * t107 + 0.2e1 * t51 * t79, -0.2e1 * t73 * t120 - 0.2e1 * t52 * t79, t50 ^ 2, t49 * t124, t104 * t124, t49 * t95, t73 * t84 ^ 2, -0.2e1 * t27 * t104 + 0.2e1 * t43 * t49, 0.2e1 * t28 * t104 + 0.2e1 * t43 * t50, 0.2e1 * t11 * t49 + 0.2e1 * t25 * t35, -0.2e1 * t12 * t49 + 0.2e1 * t25 * t36, -0.2e1 * t11 * t36 - 0.2e1 * t12 * t35, t11 ^ 2 + t12 ^ 2 + t25 ^ 2, t19 ^ 2, t19 * t126, 0.2e1 * t19 * t49, t49 * t126, t49 ^ 2, 0.2e1 * t13 * t18 - 0.2e1 * t49 * t90, 0.2e1 * t13 * t19 - 0.2e1 * t4 * t49, 0.2e1 * t18 * t5 - 0.2e1 * t2 * t49, -0.2e1 * t1 * t18 + 0.2e1 * t19 * t2, 0.2e1 * t1 * t49 - 0.2e1 * t19 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t105, t104, t79, t51, -t52, t50 * t81, -t81 * t49 + t50 * t83, -t93, -t92, 0, -pkin(2) * t49 + pkin(9) * t93 - t43 * t83, -pkin(2) * t50 + pkin(9) * t92 + t43 * t81, -t11 * t83 + t41 * t49 + (pkin(9) * t35 + t113) * t81, t12 * t83 - t42 * t49 + (pkin(9) * t36 + t112) * t81, -t42 * t35 - t41 * t36 + (-t11 * t78 - t12 * t76) * t81, t11 * t41 + t12 * t42 + t25 * t71, -t19 * t48, t18 * t48 - t19 * t47, -t19 * t83 - t48 * t49, t18 * t83 - t47 * t49, -t49 * t83, t13 * t47 + t16 * t49 + t58 * t18 + t83 * t90, -t13 * t48 - t17 * t49 + t58 * t19 + t4 * t83, -t15 * t49 + t26 * t18 + t2 * t83 + t5 * t47, -t1 * t47 - t14 * t18 + t15 * t19 - t2 * t48, -t1 * t83 + t14 * t49 - t26 * t19 + t5 * t48, t1 * t14 + t15 * t2 + t26 * t5; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t75, t81 * t121, 0, 0, 0, pkin(2) * t121, -0.2e1 * pkin(2) * t81, 0.2e1 * t75 * t118 - 0.2e1 * t41 * t83, 0.2e1 * t75 * pkin(9) * t78 + 0.2e1 * t42 * t83, 0.2e1 * (-t41 * t78 - t42 * t76) * t81, t75 * pkin(9) ^ 2 + t41 ^ 2 + t42 ^ 2, t48 ^ 2, t47 * t125, t83 * t125, t47 * t121, t83 ^ 2, -0.2e1 * t16 * t83 + 0.2e1 * t58 * t47, 0.2e1 * t17 * t83 - 0.2e1 * t58 * t48, 0.2e1 * t15 * t83 + 0.2e1 * t26 * t47, -0.2e1 * t14 * t47 - 0.2e1 * t15 * t48, -0.2e1 * t14 * t83 + 0.2e1 * t26 * t48, t14 ^ 2 + t15 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t49, -t104, t27, -t28, -pkin(3) * t35 - t76 * t98 - t112, -pkin(3) * t36 - t78 * t98 + t113 (-t35 * t78 + t36 * t76) * qJ(4) + t87, -t25 * pkin(3) + qJ(4) * t87, t19 * t57, -t18 * t57 - t19 * t56, t57 * t49, -t56 * t49, 0, t13 * t56 + t18 * t70 - t111, t13 * t57 + t19 * t70 - t109, t31 * t18 + t5 * t56 - t111, -t1 * t56 - t39 * t18 + t38 * t19 + t2 * t57, -t31 * t19 - t5 * t57 + t109, t1 * t39 + t2 * t38 + t31 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t83, 0, -t71, -t115, -pkin(9) * t103 + t88 * t76, pkin(9) * t106 + t88 * t78, t86, -pkin(3) * t71 + qJ(4) * t86, -t48 * t57, -t47 * t57 + t48 * t56, -t57 * t83, t56 * t83, 0, t70 * t47 + t58 * t56 + t110, -t70 * t48 + t58 * t57 + t108, t26 * t56 + t31 * t47 + t110, -t14 * t56 + t15 * t57 - t38 * t48 - t39 * t47, -t26 * t57 + t31 * t48 - t108, t14 * t39 + t15 * t38 + t26 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t78, -0.2e1 * pkin(3) * t76, 0.2e1 * t99 * qJ(4), qJ(4) ^ 2 * t99 + pkin(3) ^ 2, t57 ^ 2, t56 * t123, 0, 0, 0, t56 * t122, t57 * t122, 0.2e1 * t31 * t56, 0.2e1 * t38 * t57 - 0.2e1 * t39 * t56, t31 * t123, t31 ^ 2 + t38 ^ 2 + t39 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t36, 0, t25, 0, 0, 0, 0, 0, t18, t19, t18, 0, -t19, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t103, 0, t71, 0, 0, 0, 0, 0, t47, -t48, t47, 0, t48, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t76, 0, -pkin(3), 0, 0, 0, 0, 0, t56, t57, t56, 0, -t57, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t49, -t90, -t4, -t90 + 0.2e1 * t117, -pkin(5) * t19 - qJ(6) * t18, 0.2e1 * t97 + t4, -pkin(5) * t2 + qJ(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, -t83, t16, -t17, t16 - 0.2e1 * t116, pkin(5) * t48 - qJ(6) * t47, -0.2e1 * t96 + t17, -pkin(5) * t15 + qJ(6) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t56, 0, -t38, -t39, -t38, -pkin(5) * t57 - qJ(6) * t56, t39, -pkin(5) * t38 + qJ(6) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t19, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t48, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
