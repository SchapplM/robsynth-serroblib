% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:43
% EndTime: 2019-03-09 15:26:46
% DurationCPUTime: 1.07s
% Computational Cost: add. (2535->162), mult. (5738->279), div. (0->0), fcn. (5412->8), ass. (0->111)
t75 = sin(qJ(6));
t71 = t75 ^ 2;
t78 = cos(qJ(6));
t116 = -t78 ^ 2 + t71;
t103 = t116 * qJD(6);
t127 = qJD(2) + qJD(3);
t125 = pkin(7) + pkin(8);
t104 = qJD(2) * t125;
t77 = sin(qJ(2));
t54 = t77 * t104;
t80 = cos(qJ(2));
t55 = t80 * t104;
t76 = sin(qJ(3));
t79 = cos(qJ(3));
t57 = t125 * t77;
t58 = t125 * t80;
t98 = -t57 * t79 - t58 * t76;
t26 = -qJD(3) * t98 + t79 * t54 + t76 * t55;
t81 = 2 * qJD(5);
t126 = pkin(4) + pkin(9);
t52 = t76 * t77 - t79 * t80;
t36 = t127 * t52;
t53 = t76 * t80 + t77 * t79;
t37 = t127 * t53;
t73 = sin(pkin(10));
t74 = cos(pkin(10));
t24 = -t36 * t73 + t74 * t37;
t20 = -qJ(4) * t37 - qJD(4) * t52 - t26;
t97 = t57 * t76 - t79 * t58;
t27 = qJD(3) * t97 + t76 * t54 - t79 * t55;
t82 = qJ(4) * t36 - qJD(4) * t53 + t27;
t9 = t74 * t20 + t73 * t82;
t6 = -t24 * pkin(5) + t9;
t3 = t6 * t75;
t4 = t6 * t78;
t112 = qJD(6) * t78;
t29 = -qJ(4) * t53 + t98;
t30 = -qJ(4) * t52 - t97;
t23 = t29 * t73 + t30 * t74;
t34 = t74 * t52 + t53 * t73;
t16 = -pkin(5) * t34 + t23;
t124 = t16 * t112 + t3;
t22 = -t74 * t29 + t30 * t73;
t115 = pkin(2) * qJD(3);
t119 = t74 * t76;
t45 = (t73 * t79 + t119) * t115;
t123 = t22 * t45;
t122 = t24 * t75;
t121 = t24 * t78;
t35 = -t52 * t73 + t53 * t74;
t120 = t35 * t45;
t106 = t79 * t115;
t107 = t76 * t115;
t46 = t106 * t74 - t73 * t107;
t41 = qJD(5) + t46;
t66 = pkin(2) * t79 + pkin(3);
t48 = pkin(2) * t119 + t66 * t73;
t43 = qJ(5) + t48;
t118 = t43 * t112 + t41 * t75;
t63 = pkin(3) * t73 + qJ(5);
t117 = qJD(5) * t75 + t63 * t112;
t114 = qJD(6) * t16;
t113 = qJD(6) * t75;
t111 = t77 * qJD(2);
t110 = t80 * qJD(2);
t109 = -0.2e1 * pkin(1) * qJD(2);
t108 = t75 * t121;
t68 = pkin(2) * t111;
t105 = t75 * t112;
t67 = -pkin(2) * t80 - pkin(1);
t64 = -pkin(3) * t74 - pkin(4);
t31 = pkin(3) * t37 + t68;
t8 = t20 * t73 - t74 * t82;
t47 = -t73 * t76 * pkin(2) + t66 * t74;
t44 = -pkin(4) - t47;
t102 = t22 * t8 + t23 * t9;
t95 = pkin(3) * t52 + t67;
t86 = -qJ(5) * t35 + t95;
t14 = t126 * t34 + t86;
t15 = pkin(5) * t35 + t22;
t101 = t14 * t78 + t15 * t75;
t100 = t14 * t75 - t15 * t78;
t96 = -qJD(5) * t34 - t24 * t63;
t94 = t112 * t34 + t122;
t93 = t113 * t34 - t121;
t25 = -t36 * t74 - t37 * t73;
t92 = t112 * t35 + t25 * t75;
t91 = t113 * t35 - t25 * t78;
t42 = -pkin(9) + t44;
t90 = qJD(6) * (t34 * t43 - t35 * t42);
t62 = -pkin(9) + t64;
t89 = qJD(6) * (t34 * t63 - t35 * t62);
t88 = -t24 * t43 - t34 * t41 + t120;
t87 = t25 * t62 + t96;
t85 = -qJ(5) * t25 - qJD(5) * t35 + t31;
t84 = t25 * t42 + t88;
t83 = 0.2e1 * t22 * t25 - 0.2e1 * t23 * t24 - 0.2e1 * t34 * t9 + 0.2e1 * t35 * t8;
t70 = qJD(5) * t78;
t60 = -0.2e1 * t105;
t51 = 0.2e1 * t103;
t40 = t41 * t78;
t32 = t34 ^ 2;
t21 = pkin(4) * t34 + t86;
t12 = -t103 * t34 + t108;
t11 = -0.4e1 * t105 * t34 - t116 * t24;
t10 = pkin(4) * t24 + t85;
t7 = t126 * t24 + t85;
t5 = pkin(5) * t25 + t8;
t2 = -qJD(6) * t101 + t5 * t78 - t7 * t75;
t1 = qJD(6) * t100 - t5 * t75 - t7 * t78;
t13 = [0, 0, 0, 0.2e1 * t77 * t110, 0.2e1 * (-t77 ^ 2 + t80 ^ 2) * qJD(2), 0, 0, 0, t77 * t109, t80 * t109, -0.2e1 * t53 * t36, 0.2e1 * t36 * t52 - 0.2e1 * t37 * t53, 0, 0, 0, 0.2e1 * t37 * t67 + 0.2e1 * t52 * t68, -0.2e1 * t36 * t67 + 0.2e1 * t53 * t68, t83, 0.2e1 * t31 * t95 + 0.2e1 * t102, t83, -0.2e1 * t10 * t34 - 0.2e1 * t21 * t24, -0.2e1 * t10 * t35 - 0.2e1 * t21 * t25, 0.2e1 * t10 * t21 + 0.2e1 * t102, 0.2e1 * t24 * t34 * t71 + 0.2e1 * t105 * t32, -0.2e1 * t103 * t32 + 0.4e1 * t108 * t34, 0.2e1 * t122 * t35 + 0.2e1 * t34 * t92, 0.2e1 * t121 * t35 - 0.2e1 * t34 * t91, 0.2e1 * t35 * t25, -0.2e1 * t100 * t25 + 0.2e1 * t16 * t93 + 0.2e1 * t2 * t35 - 0.2e1 * t34 * t4, 0.2e1 * t1 * t35 - 0.2e1 * t101 * t25 + 0.2e1 * t16 * t94 + 0.2e1 * t3 * t34; 0, 0, 0, 0, 0, t110, -t111, 0, -pkin(7) * t110, pkin(7) * t111, 0, 0, -t36, -t37, 0, t27, t26, -t24 * t48 - t25 * t47 - t34 * t46 + t120, t23 * t46 - t47 * t8 + t48 * t9 + t123, t25 * t44 + t88, t8, t9, t23 * t41 + t43 * t9 + t44 * t8 + t123, t12, t11, -t91, -t92, 0, t75 * t90 + t78 * t84 + t124, t4 + t78 * t90 + (-t84 - t114) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t107, -0.2e1 * t106, 0, -0.2e1 * t45 * t47 + 0.2e1 * t46 * t48, 0, 0.2e1 * t45, 0.2e1 * t41, 0.2e1 * t41 * t43 + 0.2e1 * t44 * t45, t60, t51, 0, 0, 0, 0.2e1 * t118, -0.2e1 * t113 * t43 + 0.2e1 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, t27, t26 (-t24 * t73 - t25 * t74) * pkin(3) (t73 * t9 - t74 * t8) * pkin(3), t25 * t64 + t96, t8, t9, qJD(5) * t23 + t63 * t9 + t64 * t8, t12, t11, -t91, -t92, 0, t75 * t89 + t78 * t87 + t124, t4 + t78 * t89 + (-t87 - t114) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t106, 0 (-t45 * t74 + t46 * t73) * pkin(3), 0, t45, t46 + t81, qJD(5) * t43 + t41 * t63 + t45 * t64, t60, t51, 0, 0, 0, t117 + t118, t40 + t70 + (-t43 - t63) * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t63 * t81, t60, t51, 0, 0, 0, 0.2e1 * t117, -0.2e1 * t113 * t63 + 0.2e1 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t24, -t25, t10, 0, 0, 0, 0, 0, -t92, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, t8, 0, 0, 0, 0, 0, -t91, -t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t93, t25, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t112, 0, -t113 * t42 + t45 * t78, -t112 * t42 - t75 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t112, 0, -t62 * t113, -t62 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t13;
