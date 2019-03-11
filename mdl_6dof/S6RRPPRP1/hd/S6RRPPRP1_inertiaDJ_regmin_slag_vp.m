% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:47
% EndTime: 2019-03-09 08:27:51
% DurationCPUTime: 1.30s
% Computational Cost: add. (3139->191), mult. (7050->350), div. (0->0), fcn. (6829->8), ass. (0->98)
t88 = cos(pkin(9));
t81 = -t88 * pkin(2) - pkin(3);
t87 = cos(pkin(10));
t72 = -t87 * pkin(4) + t81;
t138 = 0.2e1 * t72;
t131 = cos(qJ(5));
t91 = cos(qJ(2));
t111 = -t91 * pkin(2) - pkin(1);
t86 = sin(pkin(9));
t90 = sin(qJ(2));
t67 = t86 * t90 - t88 * t91;
t69 = t86 * t91 + t88 * t90;
t47 = t67 * pkin(3) - t69 * qJ(4) + t111;
t123 = -qJ(3) - pkin(7);
t73 = t123 * t90;
t74 = t123 * t91;
t52 = t86 * t73 - t88 * t74;
t85 = sin(pkin(10));
t24 = t87 * t47 - t85 * t52;
t18 = -t87 * t69 * pkin(8) + t67 * pkin(4) + t24;
t129 = t69 * t85;
t25 = t85 * t47 + t87 * t52;
t21 = -pkin(8) * t129 + t25;
t89 = sin(qJ(5));
t137 = t131 * t21 + t89 * t18;
t106 = qJD(5) * t131;
t115 = qJD(5) * t89;
t136 = t87 * t106 - t85 * t115;
t135 = t131 * t87 - t89 * t85;
t116 = qJD(2) * t91;
t117 = qJD(2) * t90;
t63 = t116 * t88 - t117 * t86;
t126 = t87 * t63;
t62 = t69 * qJD(2);
t82 = pkin(2) * t117;
t31 = t62 * pkin(3) - t63 * qJ(4) - t69 * qJD(4) + t82;
t104 = qJD(2) * t123;
t60 = t91 * qJD(3) + t104 * t90;
t61 = -t90 * qJD(3) + t104 * t91;
t39 = t88 * t60 + t86 * t61;
t15 = t87 * t31 - t85 * t39;
t10 = t62 * pkin(4) - pkin(8) * t126 + t15;
t127 = t85 * t63;
t16 = t85 * t31 + t87 * t39;
t13 = -pkin(8) * t127 + t16;
t4 = -qJD(5) * t137 + t131 * t10 - t89 * t13;
t134 = 2 * qJD(6);
t133 = t62 * pkin(5);
t78 = t86 * pkin(2) + qJ(4);
t132 = pkin(8) + t78;
t38 = t86 * t60 - t88 * t61;
t51 = -t88 * t73 - t86 * t74;
t130 = t51 * t38;
t109 = t131 * t85;
t124 = t89 * t87;
t70 = t109 + t124;
t128 = t70 * t136;
t23 = t136 * t69 + t70 * t63;
t41 = t70 * t69;
t121 = -t136 * t41 - t70 * t23;
t119 = t85 ^ 2 + t87 ^ 2;
t118 = t62 * qJ(6);
t114 = t67 * qJD(6);
t113 = -0.2e1 * pkin(1) * qJD(2);
t110 = t89 * t132;
t105 = t131 * qJD(4);
t103 = 0.2e1 * t119 * qJD(4);
t26 = pkin(4) * t127 + t38;
t37 = pkin(4) * t129 + t51;
t101 = t15 * t87 + t16 * t85;
t100 = -t15 * t85 + t16 * t87;
t65 = t70 * qJD(5);
t22 = -t135 * t63 + t69 * t65;
t42 = t135 * t69;
t99 = t135 * t22 + t65 * t42;
t66 = t132 * t87;
t95 = t132 * t109;
t33 = qJD(5) * t95 - t87 * t105 + (qJD(4) * t85 + qJD(5) * t66) * t89;
t46 = -t85 * t110 + t131 * t66;
t98 = t33 * t67 - t46 * t62;
t34 = t66 * t106 + qJD(4) * t124 + (-qJD(5) * t110 + t105) * t85;
t45 = t89 * t66 + t95;
t97 = -t34 * t67 - t45 * t62;
t96 = t38 * t69 + t51 * t63;
t29 = t136 * t67 + t70 * t62;
t30 = t135 * t62 - t65 * t67;
t94 = t131 * t18 - t89 * t21;
t3 = -t89 * t10 - t18 * t106 + t115 * t21 - t131 * t13;
t35 = t65 * pkin(5) - qJ(6) * t136 - t70 * qJD(6);
t92 = -qJD(4) * t67 - t62 * t78 + t63 * t81;
t40 = -pkin(5) * t135 - t70 * qJ(6) + t72;
t12 = t41 * pkin(5) - t42 * qJ(6) + t37;
t7 = -t67 * pkin(5) - t94;
t6 = t67 * qJ(6) + t137;
t5 = t23 * pkin(5) + t22 * qJ(6) - t42 * qJD(6) + t26;
t2 = -t133 - t4;
t1 = t114 - t3 + t118;
t8 = [0, 0, 0, 0.2e1 * t90 * t116, 0.2e1 * (-t90 ^ 2 + t91 ^ 2) * qJD(2), 0, 0, 0, t90 * t113, t91 * t113, -0.2e1 * t39 * t67 - 0.2e1 * t52 * t62 + 0.2e1 * t96, 0.2e1 * t111 * t82 + 0.2e1 * t52 * t39 + 0.2e1 * t130, 0.2e1 * t15 * t67 + 0.2e1 * t24 * t62 + 0.2e1 * t85 * t96, -0.2e1 * t16 * t67 - 0.2e1 * t25 * t62 + 0.2e1 * t87 * t96, -0.2e1 * t101 * t69 + 0.2e1 * (-t24 * t87 - t25 * t85) * t63, 0.2e1 * t24 * t15 + 0.2e1 * t25 * t16 + 0.2e1 * t130, -0.2e1 * t42 * t22, 0.2e1 * t22 * t41 - 0.2e1 * t42 * t23, -0.2e1 * t22 * t67 + 0.2e1 * t42 * t62, -0.2e1 * t23 * t67 - 0.2e1 * t41 * t62, 0.2e1 * t67 * t62, 0.2e1 * t37 * t23 + 0.2e1 * t26 * t41 + 0.2e1 * t4 * t67 + 0.2e1 * t62 * t94, -0.2e1 * t137 * t62 - 0.2e1 * t37 * t22 + 0.2e1 * t26 * t42 + 0.2e1 * t3 * t67, 0.2e1 * t12 * t23 - 0.2e1 * t2 * t67 + 0.2e1 * t5 * t41 - 0.2e1 * t7 * t62, -0.2e1 * t1 * t41 + 0.2e1 * t2 * t42 - 0.2e1 * t7 * t22 - 0.2e1 * t6 * t23, 0.2e1 * t1 * t67 + 0.2e1 * t12 * t22 - 0.2e1 * t5 * t42 + 0.2e1 * t6 * t62, 0.2e1 * t6 * t1 + 0.2e1 * t12 * t5 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, t116, -t117, 0, -pkin(7) * t116, pkin(7) * t117 (-t62 * t86 - t63 * t88) * pkin(2) (-t38 * t88 + t39 * t86) * pkin(2), -t38 * t87 + t85 * t92, t38 * t85 + t87 * t92, t100, t38 * t81 + t100 * t78 + (-t24 * t85 + t25 * t87) * qJD(4), t136 * t42 - t22 * t70, -t99 + t121, t29, t30, 0, -t135 * t26 + t72 * t23 + t37 * t65 + t97, t136 * t37 - t72 * t22 + t26 * t70 + t98, t12 * t65 - t135 * t5 + t40 * t23 + t35 * t41 + t97, t1 * t135 + t136 * t7 + t2 * t70 - t45 * t22 - t46 * t23 + t33 * t41 + t34 * t42 - t6 * t65, -t12 * t136 + t40 * t22 - t35 * t42 - t5 * t70 - t98, t1 * t46 + t12 * t35 + t2 * t45 - t6 * t33 + t7 * t34 + t5 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t78 * t103, 0.2e1 * t128, 0.2e1 * t135 * t136 - 0.2e1 * t65 * t70, 0, 0, 0, t65 * t138, t136 * t138, -0.2e1 * t135 * t35 + 0.2e1 * t40 * t65, -0.2e1 * t135 * t33 + 0.2e1 * t136 * t45 + 0.2e1 * t34 * t70 - 0.2e1 * t46 * t65, -0.2e1 * t136 * t40 - 0.2e1 * t35 * t70, -0.2e1 * t46 * t33 + 0.2e1 * t45 * t34 + 0.2e1 * t40 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t87 * t62, -t85 * t62, -t119 * t63, t101, 0, 0, 0, 0, 0, t30, -t29, t30, t99 + t121, t29, t1 * t70 - t135 * t2 + t136 * t6 + t7 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135 * t34 + t136 * t46 - t33 * t70 + t45 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t135 * t65 + 0.2e1 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t126, 0, t38, 0, 0, 0, 0, 0, t23, -t22, t23, 0, t22, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t136, t65, 0, -t136, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, t62, t4, t3, t4 + 0.2e1 * t133, pkin(5) * t22 - t23 * qJ(6) - t41 * qJD(6), 0.2e1 * t114 - t3 + 0.2e1 * t118, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -t65, 0, -t34, t33, -t34, -pkin(5) * t136 - t65 * qJ(6) + qJD(6) * t135, -t33, -t34 * pkin(5) - t33 * qJ(6) + t46 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t136, -t65, 0, t136, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, qJ(6) * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t22, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
