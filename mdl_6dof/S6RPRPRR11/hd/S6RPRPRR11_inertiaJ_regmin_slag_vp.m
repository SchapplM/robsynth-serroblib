% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t77 = sin(pkin(6));
t96 = cos(pkin(12));
t97 = cos(pkin(7));
t87 = t97 * t96;
t76 = sin(pkin(7));
t98 = cos(pkin(6));
t91 = t98 * t76;
t117 = t77 * t87 + t91;
t75 = sin(pkin(12));
t112 = pkin(1) * t75;
t92 = t77 * t96;
t55 = qJ(2) * t92 + t98 * t112;
t40 = t117 * pkin(9) + t55;
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t109 = t75 * t77;
t94 = pkin(1) * t96;
t65 = t98 * t94;
t43 = t98 * pkin(2) + t65 + (-t97 * pkin(9) - qJ(2)) * t109;
t47 = (-pkin(9) * t75 * t76 - t96 * pkin(2) - pkin(1)) * t77;
t85 = t97 * t43 + t47 * t76;
t24 = -t81 * t40 + t85 * t83;
t111 = cos(qJ(5));
t42 = t81 * t91 + (t75 * t83 + t81 * t87) * t77;
t49 = t76 * t92 - t98 * t97;
t74 = sin(pkin(13));
t78 = cos(pkin(13));
t30 = t42 * t74 + t49 * t78;
t31 = t42 * t78 - t49 * t74;
t80 = sin(qJ(5));
t22 = t111 * t30 + t80 * t31;
t116 = -0.2e1 * t22;
t41 = t81 * t109 - t117 * t83;
t115 = 0.2e1 * t41;
t114 = -0.2e1 * t42;
t68 = -t78 * pkin(4) - pkin(3);
t113 = 0.2e1 * t68;
t23 = t111 * t31 - t80 * t30;
t79 = sin(qJ(6));
t82 = cos(qJ(6));
t14 = t82 * t23 + t41 * t79;
t110 = t14 * t79;
t108 = t76 * t81;
t107 = t76 * t83;
t106 = t79 * t22;
t58 = -t111 * t78 + t80 * t74;
t105 = t79 * t58;
t59 = t111 * t74 + t80 * t78;
t104 = t79 * t59;
t103 = t79 * t82;
t21 = t82 * t22;
t102 = t82 * t59;
t101 = pkin(10) + qJ(4);
t28 = -t76 * t43 + t97 * t47;
t17 = t41 * pkin(3) - t42 * qJ(4) + t28;
t25 = t83 * t40 + t85 * t81;
t19 = -t49 * qJ(4) + t25;
t11 = t74 * t17 + t78 * t19;
t100 = t74 ^ 2 + t78 ^ 2;
t99 = qJ(4) * t41;
t95 = -0.2e1 * t59 * t58;
t93 = t101 * t74;
t10 = t78 * t17 - t74 * t19;
t90 = -pkin(5) * t59 - pkin(11) * t58;
t89 = -t10 * t74 + t11 * t78;
t51 = -t74 * t108 + t78 * t97;
t52 = t78 * t108 + t74 * t97;
t88 = -t51 * t74 + t52 * t78;
t8 = t41 * pkin(4) - t31 * pkin(10) + t10;
t9 = -t30 * pkin(10) + t11;
t5 = t111 * t8 - t80 * t9;
t6 = t111 * t9 + t80 * t8;
t20 = t49 * pkin(3) - t24;
t12 = t30 * pkin(4) + t20;
t73 = t82 ^ 2;
t72 = t79 ^ 2;
t70 = t77 ^ 2;
t61 = t101 * t78;
t56 = t59 ^ 2;
t54 = -qJ(2) * t109 + t65;
t53 = t82 * t58;
t45 = t111 * t61 - t80 * t93;
t44 = t111 * t93 + t80 * t61;
t37 = t58 * pkin(5) - t59 * pkin(11) + t68;
t35 = t111 * t52 + t80 * t51;
t34 = -t111 * t51 + t80 * t52;
t33 = -t79 * t107 + t82 * t35;
t32 = -t82 * t107 - t79 * t35;
t27 = t79 * t37 + t82 * t45;
t26 = t82 * t37 - t79 * t45;
t13 = t79 * t23 - t41 * t82;
t7 = t22 * pkin(5) - t23 * pkin(11) + t12;
t4 = t41 * pkin(11) + t6;
t3 = -t41 * pkin(5) - t5;
t2 = t82 * t4 + t79 * t7;
t1 = -t79 * t4 + t82 * t7;
t15 = [1, 0, 0, 0.2e1 * t54 * t98 + 0.2e1 * t70 * t94, -0.2e1 * t70 * t112 - 0.2e1 * t55 * t98, 0.2e1 * (-t54 * t75 + t96 * t55) * t77, t70 * pkin(1) ^ 2 + t54 ^ 2 + t55 ^ 2, t42 ^ 2, t41 * t114, t49 * t114, t49 * t115, t49 ^ 2, -0.2e1 * t24 * t49 + 0.2e1 * t28 * t41, 0.2e1 * t25 * t49 + 0.2e1 * t28 * t42, 0.2e1 * t10 * t41 + 0.2e1 * t20 * t30, -0.2e1 * t11 * t41 + 0.2e1 * t20 * t31, -0.2e1 * t10 * t31 - 0.2e1 * t11 * t30, t10 ^ 2 + t11 ^ 2 + t20 ^ 2, t23 ^ 2, t23 * t116, t23 * t115, t41 * t116, t41 ^ 2, 0.2e1 * t12 * t22 + 0.2e1 * t5 * t41, 0.2e1 * t12 * t23 - 0.2e1 * t6 * t41, t14 ^ 2, -0.2e1 * t14 * t13, 0.2e1 * t14 * t22, t13 * t116, t22 ^ 2, 0.2e1 * t1 * t22 + 0.2e1 * t3 * t13, 0.2e1 * t3 * t14 - 0.2e1 * t2 * t22; 0, 0, 0, -t92, t109, 0, -t77 * pkin(1), 0, 0, 0, 0, 0, -t107 * t49 + t41 * t97, t108 * t49 + t42 * t97, -t107 * t30 + t51 * t41, -t107 * t31 - t52 * t41, -t52 * t30 - t51 * t31, t10 * t51 - t107 * t20 + t11 * t52, 0, 0, 0, 0, 0, -t107 * t22 - t34 * t41, -t107 * t23 - t35 * t41, 0, 0, 0, 0, 0, t34 * t13 + t32 * t22, t34 * t14 - t33 * t22; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76 ^ 2 * t83 ^ 2 + t51 ^ 2 + t52 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t41, -t49, t24, -t25, -pkin(3) * t30 - t20 * t78 - t74 * t99, -pkin(3) * t31 + t20 * t74 - t78 * t99 (-t30 * t78 + t31 * t74) * qJ(4) + t89, -t20 * pkin(3) + qJ(4) * t89, t23 * t59, -t59 * t22 - t23 * t58, t59 * t41, -t58 * t41, 0, t12 * t58 + t68 * t22 - t44 * t41, t12 * t59 + t68 * t23 - t45 * t41, t14 * t102 (-t13 * t82 - t110) * t59, t102 * t22 + t14 * t58, -t104 * t22 - t13 * t58, t22 * t58, t1 * t58 + t104 * t3 + t44 * t13 + t26 * t22, t102 * t3 + t44 * t14 - t2 * t58 - t27 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, -t108, t78 * t107, -t74 * t107, t88, pkin(3) * t107 + qJ(4) * t88, 0, 0, 0, 0, 0, -t58 * t107, -t59 * t107, 0, 0, 0, 0, 0, t104 * t34 + t32 * t58, t102 * t34 - t33 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t78, -0.2e1 * pkin(3) * t74, 0.2e1 * t100 * qJ(4), qJ(4) ^ 2 * t100 + pkin(3) ^ 2, t56, t95, 0, 0, 0, t58 * t113, t59 * t113, t73 * t56, -0.2e1 * t56 * t103, 0.2e1 * t58 * t102, t79 * t95, t58 ^ 2, 0.2e1 * t104 * t44 + 0.2e1 * t26 * t58, 0.2e1 * t102 * t44 - 0.2e1 * t27 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t31, 0, t20, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, t21, -t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t74, 0, -pkin(3), 0, 0, 0, 0, 0, t58, t59, 0, 0, 0, 0, 0, t53, -t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, t41, t5, -t6, t110, -t79 * t13 + t14 * t82, t106, t21, 0, -pkin(5) * t13 - pkin(11) * t106 - t3 * t82, -pkin(5) * t14 - pkin(11) * t21 + t3 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t35, 0, 0, 0, 0, 0, -t34 * t82, t34 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t58, 0, -t44, -t45, t79 * t102 (-t72 + t73) * t59, t105, t53, 0, -t44 * t82 + t79 * t90, t44 * t79 + t82 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t72, 0.2e1 * t103, 0, 0, 0, 0.2e1 * pkin(5) * t82, -0.2e1 * pkin(5) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, -t104, t58, t26, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t82, 0, -t79 * pkin(11), -t82 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t15;
