% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:54
% EndTime: 2022-01-20 10:05:57
% DurationCPUTime: 0.82s
% Computational Cost: add. (832->122), mult. (1947->179), div. (0->0), fcn. (1567->8), ass. (0->122)
t84 = cos(pkin(8));
t86 = sin(qJ(2));
t134 = t84 * t86;
t88 = cos(qJ(2));
t145 = t88 * pkin(1);
t78 = pkin(2) + t145;
t82 = sin(pkin(8));
t99 = pkin(1) * t134 + t82 * t78;
t53 = qJ(4) + t99;
t76 = t82 * pkin(2) + qJ(4);
t152 = t76 + t53;
t81 = sin(pkin(9));
t79 = t81 ^ 2;
t83 = cos(pkin(9));
t80 = t83 ^ 2;
t154 = (t80 / 0.2e1 + t79 / 0.2e1) * t152;
t115 = qJD(1) + qJD(2);
t85 = sin(qJ(5));
t87 = cos(qJ(5));
t62 = (t85 ^ 2 - t87 ^ 2) * t79;
t151 = t115 * t62;
t74 = t79 + t80;
t63 = t74 * t85;
t32 = t115 * t63;
t64 = t74 * t87;
t33 = t115 * t64;
t150 = t115 * t74;
t101 = t115 * t83;
t75 = t82 * t86 * pkin(1);
t59 = t84 * t145 - t75;
t149 = -t59 / 0.2e1;
t146 = t84 * pkin(2);
t136 = t83 * t85;
t103 = t84 * t78 - t75;
t93 = -t83 * pkin(4) - t81 * pkin(7) - pkin(3);
t36 = -t103 + t93;
t17 = t53 * t136 - t87 * t36;
t144 = t17 * t83;
t135 = t83 * t87;
t18 = -t53 * t135 - t85 * t36;
t143 = t18 * t83;
t57 = t93 - t146;
t28 = t76 * t136 - t87 * t57;
t142 = t28 * t83;
t29 = -t76 * t135 - t85 * t57;
t141 = t29 * t83;
t137 = t82 * t88;
t58 = (t134 + t137) * pkin(1);
t140 = t58 * t83;
t139 = t79 * t85;
t138 = t79 * t87;
t133 = t85 * t58;
t132 = t87 * t58;
t14 = (t59 * t135 + t133) * t83 + t59 * t138;
t56 = t64 * qJD(4);
t131 = t14 * qJD(2) + t56;
t25 = t74 * t59;
t70 = t74 * qJD(4);
t130 = t25 * qJD(2) + t70;
t13 = (-t59 * t136 + t132) * t83 - t59 * t139;
t55 = t63 * qJD(4);
t129 = -t13 * qJD(2) + t55;
t128 = pkin(1) * qJD(1);
t127 = pkin(1) * qJD(2);
t5 = (-pkin(3) - t103) * t58 + t53 * t25;
t126 = t5 * qJD(1);
t6 = -t53 * t139 - t144;
t125 = t6 * qJD(1);
t39 = t53 * t138;
t7 = -t39 + t143;
t124 = t7 * qJD(1);
t123 = qJD(4) * t83;
t122 = qJD(5) * t85;
t121 = qJD(5) * t87;
t12 = -t103 * t58 + t99 * t59;
t120 = t12 * qJD(1);
t119 = t13 * qJD(1);
t118 = t14 * qJD(1);
t21 = t74 * t53;
t117 = t21 * qJD(1);
t116 = t25 * qJD(1);
t114 = t85 * t138;
t113 = t83 * t122;
t112 = t83 * t121;
t111 = t85 * t123;
t110 = t87 * t123;
t109 = t85 * t149;
t108 = t87 * t149;
t106 = t132 / 0.2e1;
t65 = t76 * t138;
t105 = -t39 / 0.2e1 - t65 / 0.2e1;
t48 = t74 * t76;
t102 = pkin(1) * t115;
t100 = qJD(5) * t114;
t98 = t85 * t101;
t97 = t87 * t101;
t1 = t106 + (t109 + t29 / 0.2e1 + t18 / 0.2e1) * t83 + t105;
t20 = -t65 + t141;
t96 = -t1 * qJD(1) - t20 * qJD(2);
t19 = -t76 * t139 - t142;
t2 = (-t58 / 0.2e1 + (t76 / 0.2e1 + t53 / 0.2e1) * t79) * t85 + (t108 + t28 / 0.2e1 + t17 / 0.2e1) * t83;
t95 = -t2 * qJD(1) + t19 * qJD(2);
t89 = (t137 / 0.2e1 + t134 / 0.2e1) * pkin(1);
t8 = t89 - t154;
t94 = t8 * qJD(1) - t48 * qJD(2);
t92 = -qJD(5) + t101;
t91 = t92 * t85;
t90 = t92 * t87;
t69 = t81 * t112;
t68 = t81 * t113;
t54 = t62 * qJD(5);
t46 = t115 * t114;
t45 = t81 * t97;
t44 = t81 * t98;
t35 = t81 * t90;
t34 = t81 * t91;
t24 = t112 - t33;
t23 = t113 - t32;
t9 = t89 + t154;
t4 = -t141 / 0.2e1 - t143 / 0.2e1 + t83 * t109 + t106 - t105;
t3 = -t142 / 0.2e1 - t144 / 0.2e1 + t83 * t108 - t133 / 0.2e1 - t152 * t139 / 0.2e1;
t10 = [0, 0, 0, 0, -t86 * t127, -t88 * t127, t12 * qJD(2), -qJD(2) * t140, t130, t5 * qJD(2) + t21 * qJD(4), -t100, t54, t68, t69, 0, -t7 * qJD(5) + t129, t6 * qJD(5) + t131; 0, 0, 0, 0, -t86 * t102, -t88 * t102, t120 + (-t58 * t84 + t59 * t82) * qJD(2) * pkin(2), -t58 * t101, t116 + t130, t126 + (t58 * (-pkin(3) - t146) + t59 * t48) * qJD(2) + t9 * qJD(4), -t100, t54, t68, t69, 0, t4 * qJD(5) - t119 + t129, t3 * qJD(5) + t118 + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t150, t9 * qJD(2) + t117, 0, 0, 0, 0, 0, t32, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t151, t34, t35, 0, t4 * qJD(2) + t18 * qJD(5) - t124, t3 * qJD(2) + t17 * qJD(5) + t125; 0, 0, 0, 0, t86 * t128, t88 * t128, -t120, qJD(1) * t140, t70 - t116, -t8 * qJD(4) - t126, -t100, t54, t68, t69, 0, -t1 * qJD(5) + t119 + t55, -t2 * qJD(5) - t118 + t56; 0, 0, 0, 0, 0, 0, 0, 0, t70, t48 * qJD(4), -t100, t54, t68, t69, 0, -t20 * qJD(5) + t55, t19 * qJD(5) + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t150, -t94, 0, 0, 0, 0, 0, t32, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t151, t34, t35, 0, t29 * qJD(5) + t96, t28 * qJD(5) + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81 * t121, t81 * t122; 0, 0, 0, 0, 0, 0, 0, 0, -t150, t8 * qJD(2) - t117, 0, 0, 0, 0, 0, t23, t24; 0, 0, 0, 0, 0, 0, 0, 0, -t150, t94, 0, 0, 0, 0, 0, t23, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t151, -t44, -t45, 0, t1 * qJD(2) - t111 + t124, t2 * qJD(2) - t110 - t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t151, -t44, -t45, 0, -t96 - t111, -t95 - t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
