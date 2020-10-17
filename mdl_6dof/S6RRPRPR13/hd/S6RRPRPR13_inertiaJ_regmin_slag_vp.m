% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR13_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:38:26
% EndTime: 2019-05-06 16:38:32
% DurationCPUTime: 1.22s
% Computational Cost: add. (1373->169), mult. (3116->334), div. (0->0), fcn. (3464->10), ass. (0->101)
t75 = sin(qJ(4));
t114 = -0.2e1 * t75;
t118 = -2 * pkin(2);
t71 = sin(pkin(6));
t79 = cos(qJ(2));
t100 = t71 * t79;
t73 = cos(pkin(6));
t78 = cos(qJ(4));
t40 = t78 * t100 + t73 * t75;
t117 = -0.2e1 * t40;
t72 = cos(pkin(11));
t63 = -t72 * pkin(5) - pkin(4);
t116 = 0.2e1 * t63;
t115 = 0.2e1 * t71;
t113 = 2 * qJ(3);
t80 = -pkin(2) - pkin(9);
t76 = sin(qJ(2));
t112 = pkin(1) * t76;
t111 = pkin(1) * t79;
t110 = t78 * pkin(4);
t60 = t71 * t76;
t55 = pkin(8) * t60;
t88 = -pkin(2) - t111;
t22 = pkin(3) * t60 + t55 + (-pkin(9) + t88) * t73;
t87 = -qJ(3) * t76 - pkin(1);
t30 = (t80 * t79 + t87) * t71;
t18 = t75 * t22 + t78 * t30;
t11 = qJ(5) * t60 + t18;
t43 = pkin(8) * t100 + t73 * t112;
t64 = t73 * qJ(3);
t33 = -t64 - t43;
t29 = pkin(3) * t100 - t33;
t41 = -t75 * t100 + t73 * t78;
t19 = t40 * pkin(4) - t41 * qJ(5) + t29;
t70 = sin(pkin(11));
t6 = t72 * t11 + t70 * t19;
t17 = t78 * t22 - t75 * t30;
t12 = -pkin(4) * t60 - t17;
t109 = t12 * t70;
t108 = t12 * t72;
t107 = t40 * t75;
t74 = sin(qJ(6));
t77 = cos(qJ(6));
t47 = t74 * t70 - t77 * t72;
t106 = t47 * t75;
t48 = t77 * t70 + t74 * t72;
t105 = t48 * t75;
t66 = t71 ^ 2;
t104 = t66 * t79;
t69 = t78 ^ 2;
t103 = t69 * t80;
t102 = t70 * t78;
t101 = t70 * t80;
t61 = t72 * t78;
t99 = t73 * t79;
t98 = t75 * t80;
t97 = t78 * t47;
t37 = t78 * t48;
t96 = t78 * t80;
t95 = pkin(10) + qJ(5);
t49 = t75 * pkin(4) - t78 * qJ(5) + qJ(3);
t32 = t70 * t49 + t72 * t98;
t94 = t70 ^ 2 + t72 ^ 2;
t68 = t75 ^ 2;
t93 = -t68 - t69;
t92 = qJ(5) * t40;
t91 = 0.2e1 * t60;
t90 = t75 * t60;
t89 = t80 * t60;
t5 = -t70 * t11 + t72 * t19;
t86 = t94 * qJ(5);
t85 = -t5 * t70 + t6 * t72;
t84 = -qJ(5) * t75 - t110;
t25 = t41 * t70 - t72 * t60;
t26 = t41 * t72 + t70 * t60;
t83 = -t25 * t72 + t26 * t70;
t45 = t72 * t49;
t31 = -t70 * t98 + t45;
t82 = -t31 * t70 + t32 * t72;
t59 = t66 * t76 ^ 2;
t53 = t78 * t60;
t51 = t95 * t72;
t50 = t95 * t70;
t46 = (pkin(5) * t70 - t80) * t78;
t42 = pkin(1) * t99 - t55;
t35 = t88 * t73 + t55;
t34 = (-pkin(2) * t79 + t87) * t71;
t28 = -t74 * t50 + t77 * t51;
t27 = -t77 * t50 - t74 * t51;
t23 = -pkin(10) * t102 + t32;
t20 = -pkin(10) * t61 + t45 + (pkin(5) - t101) * t75;
t14 = -t74 * t25 + t77 * t26;
t13 = t77 * t25 + t74 * t26;
t10 = t74 * t20 + t77 * t23;
t9 = t77 * t20 - t74 * t23;
t7 = t25 * pkin(5) + t12;
t4 = -t25 * pkin(10) + t6;
t3 = t40 * pkin(5) - t26 * pkin(10) + t5;
t2 = t74 * t3 + t77 * t4;
t1 = t77 * t3 - t74 * t4;
t8 = [1, 0, 0, t59, 0.2e1 * t76 * t104, t73 * t91, t99 * t115, t73 ^ 2, 0.2e1 * pkin(1) * t104 + 0.2e1 * t42 * t73, -0.2e1 * t66 * t112 - 0.2e1 * t43 * t73 (-t33 * t79 + t35 * t76) * t115, 0.2e1 * t34 * t100 + 0.2e1 * t35 * t73, -0.2e1 * t33 * t73 - 0.2e1 * t34 * t60, t33 ^ 2 + t34 ^ 2 + t35 ^ 2, t41 ^ 2, t41 * t117, t41 * t91, t60 * t117, t59, 0.2e1 * t17 * t60 + 0.2e1 * t29 * t40, -0.2e1 * t18 * t60 + 0.2e1 * t29 * t41, 0.2e1 * t12 * t25 + 0.2e1 * t5 * t40, 0.2e1 * t12 * t26 - 0.2e1 * t6 * t40, -0.2e1 * t6 * t25 - 0.2e1 * t5 * t26, t12 ^ 2 + t5 ^ 2 + t6 ^ 2, t14 ^ 2, -0.2e1 * t14 * t13, 0.2e1 * t14 * t40, t13 * t117, t40 ^ 2, 0.2e1 * t1 * t40 + 0.2e1 * t7 * t13, 0.2e1 * t7 * t14 - 0.2e1 * t2 * t40; 0, 0, 0, 0, 0, t60, t100, t73, t42, -t43 (-pkin(2) * t76 + qJ(3) * t79) * t71, t55 + (t118 - t111) * t73, 0.2e1 * t64 + t43, -t35 * pkin(2) - t33 * qJ(3), t41 * t78, -t78 * t40 - t41 * t75, t53, -t90, 0, qJ(3) * t40 + t29 * t75 + t78 * t89, qJ(3) * t41 + t29 * t78 - t75 * t89, t31 * t40 + t5 * t75 + (-t25 * t80 + t109) * t78, -t32 * t40 - t6 * t75 + (-t26 * t80 + t108) * t78, -t32 * t25 - t31 * t26 + (-t5 * t72 - t6 * t70) * t78, -t12 * t96 + t5 * t31 + t6 * t32, -t14 * t97, t13 * t97 - t14 * t37, t14 * t75 - t40 * t97, -t13 * t75 - t37 * t40, t107, t1 * t75 + t46 * t13 + t7 * t37 + t9 * t40, -t10 * t40 + t46 * t14 - t2 * t75 - t7 * t97; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t118, t113, pkin(2) ^ 2 + qJ(3) ^ 2, t69, t78 * t114, 0, 0, 0, t75 * t113, t78 * t113, -0.2e1 * t69 * t101 + 0.2e1 * t31 * t75, -0.2e1 * t72 * t103 - 0.2e1 * t32 * t75, 0.2e1 * (-t31 * t72 - t32 * t70) * t78, t69 * (t80 ^ 2) + t31 ^ 2 + t32 ^ 2, t97 ^ 2, 0.2e1 * t97 * t37, t97 * t114, t37 * t114, t68, 0.2e1 * t46 * t37 + 0.2e1 * t9 * t75, -0.2e1 * t10 * t75 - 0.2e1 * t46 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t73, 0, t35, 0, 0, 0, 0, 0, t53, -t90, -t70 * t107 - t78 * t25, -t72 * t107 - t78 * t26, t83 * t75, -t12 * t78 + t85 * t75, 0, 0, 0, 0, 0, -t105 * t40 - t78 * t13, t106 * t40 - t78 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, t93 * t70, t93 * t72, 0, t82 * t75 + t103, 0, 0, 0, 0, 0, -t105 * t75 - t78 * t37, t106 * t75 + t78 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94 * t68 + t69, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, t60, t17, -t18, -pkin(4) * t25 - t70 * t92 - t108, -pkin(4) * t26 - t72 * t92 + t109, t83 * qJ(5) + t85, -t12 * pkin(4) + t85 * qJ(5), t14 * t48, -t48 * t13 - t14 * t47, t48 * t40, -t47 * t40, 0, t63 * t13 + t27 * t40 + t7 * t47, t63 * t14 - t28 * t40 + t7 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t75, 0, t96, -t98, t84 * t70 + t72 * t96, -t70 * t96 + t84 * t72, t82, pkin(4) * t96 + t82 * qJ(5), -t97 * t48, -t48 * t37 + t47 * t97, t105, -t106, 0, t27 * t75 + t63 * t37 + t46 * t47, -t28 * t75 + t46 * t48 - t63 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t75, t61, -t102, t94 * t75, t75 * t86 + t110, 0, 0, 0, 0, 0, -t97, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t72, -0.2e1 * pkin(4) * t70, 0.2e1 * t86, t94 * qJ(5) ^ 2 + pkin(4) ^ 2, t48 ^ 2, -0.2e1 * t48 * t47, 0, 0, 0, t47 * t116, t48 * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26, 0, t12, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t61, 0, -t96, 0, 0, 0, 0, 0, t37, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t70, 0, -pkin(4), 0, 0, 0, 0, 0, t47, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t40, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t37, t75, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t47, 0, t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
