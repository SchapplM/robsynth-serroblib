% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:29:51
% EndTime: 2019-03-09 04:29:56
% DurationCPUTime: 1.42s
% Computational Cost: add. (2053->204), mult. (4781->354), div. (0->0), fcn. (4001->8), ass. (0->107)
t96 = cos(qJ(3));
t128 = t96 * qJD(3);
t93 = sin(qJ(4));
t119 = t93 * t128;
t95 = cos(qJ(4));
t132 = qJD(4) * t95;
t94 = sin(qJ(3));
t150 = t94 * t132 + t119;
t149 = -0.4e1 * t94;
t107 = -t96 * pkin(3) - t94 * pkin(8);
t82 = -cos(pkin(9)) * pkin(1) - pkin(2);
t59 = t107 + t82;
t47 = t93 * t59;
t142 = t95 * t96;
t80 = sin(pkin(9)) * pkin(1) + pkin(7);
t64 = t80 * t142;
t148 = -t47 - t64;
t134 = cos(pkin(10));
t112 = t134 * t95;
t133 = qJD(4) * t93;
t91 = sin(pkin(10));
t52 = qJD(4) * t112 - t91 * t133;
t88 = t94 ^ 2;
t109 = (-t96 ^ 2 + t88) * qJD(3);
t89 = t95 ^ 2;
t138 = t93 ^ 2 - t89;
t110 = t138 * qJD(4);
t147 = 2 * qJD(6);
t106 = pkin(3) * t94 - pkin(8) * t96;
t65 = t106 * qJD(3);
t140 = -t59 * t132 - t93 * t65;
t143 = t94 * t95;
t14 = (-qJ(5) * qJD(4) - qJD(3) * t80) * t143 + (-qJD(5) * t94 + (-qJ(5) * qJD(3) - qJD(4) * t80) * t96) * t93 - t140;
t129 = t95 * qJD(5);
t135 = qJ(5) * t95;
t136 = qJ(5) * t94;
t86 = t94 * qJD(3);
t116 = t80 * t86;
t139 = t93 * t116 + t95 * t65;
t9 = -t94 * t129 + (pkin(4) * t94 - t96 * t135) * qJD(3) + (-t64 + (-t59 + t136) * t93) * qJD(4) + t139;
t4 = t134 * t14 + t91 * t9;
t146 = t80 * t93;
t145 = t91 * t93;
t144 = t93 * t94;
t141 = -qJ(5) - pkin(8);
t48 = t95 * t59;
t33 = -t94 * t135 + t48 + (-pkin(4) - t146) * t96;
t37 = -t93 * t136 - t148;
t13 = t134 * t37 + t91 * t33;
t50 = pkin(4) * t144 + t94 * t80;
t131 = qJD(4) * t96;
t46 = t94 * t112 - t91 * t144;
t130 = t46 * qJD(6);
t127 = t96 * qJD(6);
t126 = -0.2e1 * pkin(3) * qJD(4);
t125 = qJ(6) * t86 + t4;
t124 = 0.2e1 * qJD(3) * t82;
t66 = t80 * t128;
t40 = t150 * pkin(4) + t66;
t85 = pkin(4) * t133;
t123 = t93 * t131;
t121 = t95 * t131;
t118 = t93 * t132;
t117 = t94 * t128;
t115 = t95 * t128;
t84 = -t95 * pkin(4) - pkin(3);
t3 = t134 * t9 - t91 * t14;
t111 = qJD(4) * t141;
t49 = t93 * t111 + t129;
t99 = -t93 * qJD(5) + t95 * t111;
t34 = -t134 * t99 + t91 * t49;
t35 = t134 * t49 + t91 * t99;
t113 = t134 * t93;
t71 = t141 * t95;
t41 = -t141 * t113 - t91 * t71;
t42 = -t134 * t71 + t141 * t145;
t114 = t41 * t34 + t42 * t35;
t108 = t93 * t115;
t105 = t134 * t128;
t12 = t134 * t33 - t91 * t37;
t62 = t91 * t95 + t113;
t31 = -t93 * t105 - t91 * t115 - t52 * t94;
t51 = t62 * qJD(4);
t32 = -t95 * t105 + t91 * t119 + t94 * t51;
t45 = t62 * t94;
t103 = -t31 * t41 - t32 * t42 + t45 * t34 + t46 * t35;
t102 = t42 * t31 - t41 * t32 + t34 * t46 - t35 * t45;
t61 = -t112 + t145;
t100 = -t31 * t62 + t32 * t61 + t45 * t52 - t46 * t51;
t54 = t95 * t86 + t123;
t98 = 0.2e1 * t34 * t62 - 0.2e1 * t35 * t61 + 0.2e1 * t41 * t52 - 0.2e1 * t42 * t51;
t97 = -0.2e1 * t45 * t31 - 0.2e1 * t46 * t32 - 0.2e1 * t117;
t81 = -t134 * pkin(4) - pkin(5);
t78 = t91 * pkin(4) + qJ(6);
t56 = t93 * t86 - t121;
t53 = t94 * t133 - t115;
t38 = t61 * pkin(5) - t62 * qJ(6) + t84;
t26 = t51 * pkin(5) - t52 * qJ(6) - t62 * qJD(6) + t85;
t25 = t45 * pkin(5) - t46 * qJ(6) + t50;
t19 = t148 * qJD(4) + t139;
t18 = t54 * t80 + t140;
t10 = t96 * pkin(5) - t12;
t8 = -t96 * qJ(6) + t13;
t5 = -t31 * pkin(5) + t32 * qJ(6) - t130 + t40;
t2 = -pkin(5) * t86 - t3;
t1 = t125 - t127;
t6 = [0, 0, 0, 0, 0.2e1 * t117, -0.2e1 * t109, 0, 0, 0, t94 * t124, t96 * t124, 0.2e1 * t89 * t117 - 0.2e1 * t88 * t118, t108 * t149 + 0.2e1 * t88 * t110, 0.2e1 * t95 * t109 + 0.2e1 * t94 * t123, -0.2e1 * t93 * t109 + 0.2e1 * t94 * t121, -0.2e1 * t117, 0.2e1 * t48 * t86 - 0.2e1 * t19 * t96 + 0.2e1 * (t93 * t117 + t88 * t132) * t80, -0.2e1 * t88 * t80 * t133 - 0.2e1 * t18 * t96 + 0.2e1 * (-t47 + t64) * t86, 0.2e1 * t12 * t32 + 0.2e1 * t13 * t31 - 0.2e1 * t3 * t46 - 0.2e1 * t4 * t45, 0.2e1 * t12 * t3 + 0.2e1 * t13 * t4 + 0.2e1 * t50 * t40, -0.2e1 * t10 * t86 + 0.2e1 * t2 * t96 - 0.2e1 * t25 * t31 + 0.2e1 * t5 * t45, -0.2e1 * t1 * t45 - 0.2e1 * t10 * t32 + 0.2e1 * t2 * t46 + 0.2e1 * t8 * t31, -0.2e1 * t1 * t96 + 0.2e1 * t25 * t32 - 0.2e1 * t5 * t46 + 0.2e1 * t8 * t86, 0.2e1 * t8 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t25 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t31 - t13 * t32 - t3 * t45 + t4 * t46 - t40 * t96 + t50 * t86, 0, 0, 0, t1 * t46 - t10 * t31 + t2 * t45 + t25 * t86 - t8 * t32 - t5 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, t97; 0, 0, 0, 0, 0, 0, t128, -t86, 0, -t66, t116, -t94 * t110 + t108, t118 * t149 - t138 * t128, t56, t54, 0 (pkin(8) * t142 + (-pkin(3) * t95 + t146) * t94) * qJD(4) + (t107 * t93 - t64) * qJD(3) (t106 * t93 + t80 * t143) * qJD(4) + (t107 * t95 + t96 * t146) * qJD(3), -t12 * t52 - t13 * t51 - t3 * t62 - t4 * t61 + t102, -t12 * t34 + t13 * t35 - t3 * t41 + t4 * t42 + t40 * t84 + t50 * t85, t25 * t51 + t26 * t45 - t38 * t31 + t34 * t96 - t41 * t86 + t5 * t61, -t1 * t61 + t10 * t52 + t2 * t62 - t8 * t51 + t102, -t25 * t52 - t26 * t46 + t38 * t32 - t35 * t96 + t42 * t86 - t5 * t62, t1 * t42 + t10 * t34 + t2 * t41 + t25 * t26 + t8 * t35 + t5 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t128, 0, 0, 0, 0, 0, -t54, t56, t100, -pkin(4) * t123 + t84 * t86 + t103, -t96 * t51 + t61 * t86, t100, t96 * t52 - t62 * t86, -t96 * t26 + t38 * t86 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t118, -0.2e1 * t110, 0, 0, 0, t93 * t126, t95 * t126, t98, 0.2e1 * t84 * t85 + 0.2e1 * t114, 0.2e1 * t26 * t61 + 0.2e1 * t38 * t51, t98, -0.2e1 * t26 * t62 - 0.2e1 * t38 * t52, 0.2e1 * t38 * t26 + 0.2e1 * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t150, t86, t19, t18 (t134 * t32 + t31 * t91) * pkin(4) (t134 * t3 + t4 * t91) * pkin(4) (pkin(5) - t81) * t86 + t3, -qJD(6) * t45 + t78 * t31 - t81 * t32, t78 * t86 + t125 - 0.2e1 * t127, t8 * qJD(6) + t1 * t78 + t2 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, t53, 0 (t134 * t31 - t32 * t91) * pkin(4), t31, 0, -t32, -t31 * t81 - t32 * t78 + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, -t133, 0, -pkin(8) * t132, pkin(8) * t133 (-t134 * t52 - t51 * t91) * pkin(4) (-t134 * t34 + t35 * t91) * pkin(4), -t34, -qJD(6) * t61 - t78 * t51 + t81 * t52, t35, t42 * qJD(6) + t34 * t81 + t35 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t78 * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t31, 0, t32, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t51, 0, -t52, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t32, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
