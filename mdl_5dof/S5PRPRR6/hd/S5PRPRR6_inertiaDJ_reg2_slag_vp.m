% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:54
% EndTime: 2019-12-05 15:57:58
% DurationCPUTime: 1.14s
% Computational Cost: add. (1164->138), mult. (3051->273), div. (0->0), fcn. (3117->10), ass. (0->95)
t101 = pkin(7) + qJ(3);
t43 = cos(pkin(10));
t31 = t101 * t43;
t41 = sin(pkin(10));
t112 = cos(qJ(4));
t71 = t101 * t112;
t82 = t112 * qJD(3);
t46 = sin(qJ(4));
t92 = t46 * qJD(3);
t95 = qJD(4) * t46;
t50 = (-qJD(4) * t71 - t92) * t41 - t31 * t95 + t43 * t82;
t29 = t112 * t41 + t43 * t46;
t35 = -pkin(3) * t43 - pkin(2);
t61 = t112 * t43 - t41 * t46;
t60 = -pkin(4) * t61 - pkin(8) * t29 + t35;
t117 = -qJD(5) * t60 - t50;
t45 = sin(qJ(5));
t39 = t45 ^ 2;
t48 = cos(qJ(5));
t40 = t48 ^ 2;
t99 = t39 - t40;
t81 = qJD(5) * t99;
t42 = sin(pkin(5));
t47 = sin(qJ(2));
t106 = t42 * t47;
t44 = cos(pkin(5));
t22 = -t106 * t41 + t43 * t44;
t23 = t106 * t43 + t41 * t44;
t62 = t112 * t22 - t23 * t46;
t58 = t62 * qJD(4);
t15 = t112 * t23 + t22 * t46;
t49 = cos(qJ(2));
t105 = t42 * t49;
t89 = t45 * t105;
t11 = t15 * t48 - t89;
t97 = qJD(2) * t47;
t32 = t42 * t97;
t59 = t42 * t61;
t96 = qJD(2) * t49;
t51 = t59 * t96 + t58;
t66 = t105 * t48 + t15 * t45;
t3 = qJD(5) * t66 - t32 * t45 - t48 * t51;
t93 = qJD(5) * t48;
t4 = qJD(5) * t89 - t15 * t93 + t32 * t48 - t45 * t51;
t116 = (-t11 * t48 - t45 * t66) * qJD(5) + t3 * t45 - t4 * t48;
t85 = t46 * t101;
t20 = t112 * t31 - t41 * t85;
t83 = qJD(4) * t112;
t24 = t41 * t95 - t43 * t83;
t25 = t29 * qJD(4);
t78 = pkin(4) * t25 + pkin(8) * t24;
t94 = qJD(5) * t45;
t1 = t117 * t48 + t20 * t94 - t45 * t78;
t2 = t117 * t45 - t20 * t93 + t48 * t78;
t6 = -t20 * t45 + t48 * t60;
t7 = t48 * t20 + t45 * t60;
t115 = t1 * t45 - t2 * t48 + (t45 * t6 - t48 * t7) * qJD(5);
t114 = 0.2e1 * t25;
t87 = t42 * t96;
t8 = qJD(4) * t15 + t29 * t87;
t113 = t62 * t8;
t13 = t31 * t83 + t43 * t92 + (-qJD(4) * t85 + t82) * t41;
t19 = t46 * t31 + t41 * t71;
t111 = t19 * t13;
t110 = t29 * t24;
t109 = t29 * t45;
t108 = t29 * t48;
t107 = t42 ^ 2 * t47;
t104 = t45 * t25;
t103 = t48 * t24;
t102 = t48 * t25;
t100 = t41 ^ 2 + t43 ^ 2;
t98 = qJD(2) * t42;
t91 = t61 * t114;
t90 = -0.2e1 * pkin(4) * qJD(5);
t88 = t45 * t103;
t86 = t45 * t93;
t84 = t100 * t49;
t26 = t29 ^ 2;
t80 = t26 * t86;
t79 = 0.2e1 * t100 * qJD(3);
t77 = pkin(4) * t24 - pkin(8) * t25;
t76 = pkin(4) * t29 - pkin(8) * t61;
t73 = t45 * t7 + t48 * t6;
t70 = -t13 * t62 + t8 * t19;
t69 = t11 * t45 - t48 * t66;
t67 = -t22 * t41 + t23 * t43;
t65 = -t24 * t45 + t29 * t93;
t64 = -t29 * t94 - t103;
t63 = -t61 * t93 + t104;
t53 = -qJD(5) * t73 - t1 * t48 - t2 * t45;
t52 = -qJD(5) * t69 - t3 * t48 - t4 * t45;
t18 = t61 * t94 + t102;
t9 = t29 * t81 + t88;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t42 * t67 - t107) * t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t113 + 0.2e1 * t15 * t58 + 0.2e1 * (t15 * t59 - t107) * t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t11 * t3 - 0.2e1 * t4 * t66 - 0.2e1 * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t87, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t32, t41 * t32, t84 * t98, t67 * qJD(3) + (-pkin(2) * t47 + qJ(3) * t84) * t98, 0, 0, 0, 0, 0, 0, (-t25 * t49 - t61 * t97) * t42, (t24 * t49 + t29 * t97) * t42, -t15 * t25 + t24 * t62 + t8 * t29 + t51 * t61, t15 * t50 + t20 * t51 + t32 * t35 + t70, 0, 0, 0, 0, 0, 0, t109 * t8 - t25 * t66 - t4 * t61 - t62 * t65, t108 * t8 - t11 * t25 - t3 * t61 - t62 * t64, t116 * t29 + t69 * t24, -t1 * t11 - t2 * t66 - t3 * t7 + t4 * t6 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, qJ(3) * t79, -0.2e1 * t110, -0.2e1 * t24 * t61 - 0.2e1 * t25 * t29, 0, -t91, 0, 0, t35 * t114, -0.2e1 * t35 * t24, 0.2e1 * t13 * t29 - 0.2e1 * t19 * t24 - 0.2e1 * t20 * t25 + 0.2e1 * t50 * t61, 0.2e1 * t20 * t50 + 0.2e1 * t111, -0.2e1 * t110 * t40 - 0.2e1 * t80, 0.2e1 * t26 * t81 + 0.4e1 * t29 * t88, 0.2e1 * t102 * t29 - 0.2e1 * t61 * t64, -0.2e1 * t110 * t39 + 0.2e1 * t80, -0.2e1 * t104 * t29 + 0.2e1 * t61 * t65, -t91, 0.2e1 * t109 * t13 + 0.2e1 * t19 * t65 - 0.2e1 * t2 * t61 + 0.2e1 * t6 * t25, -0.2e1 * t1 * t61 + 0.2e1 * t108 * t13 + 0.2e1 * t19 * t64 - 0.2e1 * t7 * t25, 0.2e1 * t115 * t29 + 0.2e1 * t73 * t24, -0.2e1 * t1 * t7 + 0.2e1 * t2 * t6 + 0.2e1 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t63, (t39 + t40) * t24, -t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t61 * t87 - t58, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t8 - t62 * t94, t45 * t8 - t62 * t93, t52, -t8 * pkin(4) + pkin(8) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, -t25, 0, -t13, -t50, 0, 0, -t9, t24 * t99 - 0.4e1 * t29 * t86, t63, t9, t18, 0, -t13 * t48 + t77 * t45 + (t19 * t45 - t48 * t76) * qJD(5), t13 * t45 + t77 * t48 + (t19 * t48 + t45 * t76) * qJD(5), t53, -t13 * pkin(4) + pkin(8) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t86, -0.2e1 * t81, 0, -0.2e1 * t86, 0, 0, t45 * t90, t48 * t90, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t65, t25, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t93, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, -t94, 0, -pkin(8) * t93, pkin(8) * t94, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
