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
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:57:38
% EndTime: 2020-01-03 11:57:40
% DurationCPUTime: 0.77s
% Computational Cost: add. (838->124), mult. (1971->182), div. (0->0), fcn. (1587->8), ass. (0->125)
t85 = cos(pkin(8));
t87 = sin(qJ(2));
t138 = t85 * t87;
t89 = cos(qJ(2));
t148 = t89 * pkin(1);
t79 = pkin(2) + t148;
t83 = sin(pkin(8));
t100 = pkin(1) * t138 + t83 * t79;
t54 = qJ(4) + t100;
t77 = t83 * pkin(2) + qJ(4);
t155 = t77 + t54;
t82 = sin(pkin(9));
t80 = t82 ^ 2;
t84 = cos(pkin(9));
t81 = t84 ^ 2;
t157 = (t81 / 0.2e1 + t80 / 0.2e1) * t155;
t117 = qJD(1) + qJD(2);
t86 = sin(qJ(5));
t88 = cos(qJ(5));
t63 = (t86 ^ 2 - t88 ^ 2) * t80;
t154 = t117 * t63;
t75 = t80 + t81;
t64 = t75 * t86;
t32 = t117 * t64;
t65 = t75 * t88;
t33 = t117 * t65;
t153 = t117 * t75;
t102 = t117 * t84;
t76 = t83 * t87 * pkin(1);
t60 = t85 * t148 - t76;
t152 = -t60 / 0.2e1;
t149 = t85 * pkin(2);
t140 = t84 * t86;
t104 = t85 * t79 - t76;
t94 = -t84 * pkin(4) - t82 * pkin(7) - pkin(3);
t36 = -t104 + t94;
t17 = t54 * t140 - t88 * t36;
t147 = t17 * t84;
t139 = t84 * t88;
t18 = -t54 * t139 - t86 * t36;
t146 = t18 * t84;
t58 = t94 - t149;
t28 = t77 * t140 - t88 * t58;
t145 = t28 * t84;
t29 = -t77 * t139 - t86 * t58;
t144 = t29 * t84;
t143 = t80 * t86;
t142 = t80 * t88;
t141 = t83 * t89;
t59 = (t138 + t141) * pkin(1);
t137 = t86 * t59;
t136 = t88 * t59;
t14 = (t60 * t139 + t137) * t84 + t60 * t142;
t57 = t65 * qJD(4);
t135 = t14 * qJD(2) + t57;
t25 = t75 * t60;
t71 = t75 * qJD(4);
t134 = t25 * qJD(2) + t71;
t13 = (-t60 * t140 + t136) * t84 - t60 * t143;
t56 = t64 * qJD(4);
t133 = -t13 * qJD(2) + t56;
t132 = pkin(1) * qJD(1);
t131 = pkin(1) * qJD(2);
t5 = (-pkin(3) - t104) * t59 + t54 * t25;
t130 = t5 * qJD(1);
t6 = -t54 * t143 - t147;
t129 = t6 * qJD(1);
t39 = t54 * t142;
t7 = -t39 + t146;
t128 = t7 * qJD(1);
t127 = qJD(1) * t59;
t126 = qJD(2) * t59;
t125 = qJD(4) * t84;
t124 = qJD(5) * t86;
t123 = qJD(5) * t88;
t12 = t100 * t60 - t104 * t59;
t122 = t12 * qJD(1);
t121 = t13 * qJD(1);
t120 = t14 * qJD(1);
t21 = t75 * t54;
t119 = t21 * qJD(1);
t118 = t25 * qJD(1);
t116 = t86 * t142;
t115 = t84 * t124;
t114 = t84 * t123;
t113 = t82 * t127;
t112 = t86 * t125;
t111 = t88 * t125;
t110 = t86 * t152;
t109 = t88 * t152;
t107 = t136 / 0.2e1;
t66 = t77 * t142;
t106 = -t39 / 0.2e1 - t66 / 0.2e1;
t48 = t75 * t77;
t103 = pkin(1) * t117;
t101 = qJD(5) * t116;
t99 = t86 * t102;
t98 = t88 * t102;
t1 = t107 + (t110 + t29 / 0.2e1 + t18 / 0.2e1) * t84 + t106;
t20 = -t66 + t144;
t97 = -t1 * qJD(1) - t20 * qJD(2);
t19 = -t77 * t143 - t145;
t2 = (-t59 / 0.2e1 + (t77 / 0.2e1 + t54 / 0.2e1) * t80) * t86 + (t109 + t28 / 0.2e1 + t17 / 0.2e1) * t84;
t96 = -t2 * qJD(1) + t19 * qJD(2);
t90 = (t141 / 0.2e1 + t138 / 0.2e1) * pkin(1);
t8 = t90 - t157;
t95 = t8 * qJD(1) - t48 * qJD(2);
t93 = -qJD(5) + t102;
t92 = t93 * t86;
t91 = t93 * t88;
t70 = t82 * t114;
t69 = t82 * t115;
t55 = t63 * qJD(5);
t51 = t82 * t126;
t46 = t117 * t116;
t45 = t82 * t98;
t44 = t82 * t99;
t35 = t82 * t91;
t34 = t82 * t92;
t24 = t114 - t33;
t23 = t115 - t32;
t9 = t90 + t157;
t4 = -t144 / 0.2e1 - t146 / 0.2e1 + t84 * t110 + t107 - t106;
t3 = -t145 / 0.2e1 - t147 / 0.2e1 + t84 * t109 - t137 / 0.2e1 - t155 * t143 / 0.2e1;
t10 = [0, 0, 0, 0, -t87 * t131, -t89 * t131, t12 * qJD(2), -t84 * t126, t51, t134, t5 * qJD(2) + t21 * qJD(4), -t101, t55, t69, t70, 0, -t7 * qJD(5) + t133, t6 * qJD(5) + t135; 0, 0, 0, 0, -t87 * t103, -t89 * t103, t122 + (-t59 * t85 + t60 * t83) * qJD(2) * pkin(2), -t59 * t102, t51 + t113, t118 + t134, t130 + (t59 * (-pkin(3) - t149) + t60 * t48) * qJD(2) + t9 * qJD(4), -t101, t55, t69, t70, 0, t4 * qJD(5) - t121 + t133, t3 * qJD(5) + t120 + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t9 * qJD(2) + t119, 0, 0, 0, 0, 0, t32, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t154, t34, t35, 0, t4 * qJD(2) + t18 * qJD(5) - t128, t3 * qJD(2) + t17 * qJD(5) + t129; 0, 0, 0, 0, t87 * t132, t89 * t132, -t122, t84 * t127, -t113, t71 - t118, -t8 * qJD(4) - t130, -t101, t55, t69, t70, 0, -t1 * qJD(5) + t121 + t56, -t2 * qJD(5) - t120 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t48 * qJD(4), -t101, t55, t69, t70, 0, -t20 * qJD(5) + t56, t19 * qJD(5) + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, -t95, 0, 0, 0, 0, 0, t32, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t154, t34, t35, 0, t29 * qJD(5) + t97, t28 * qJD(5) + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82 * t123, t82 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t8 * qJD(2) - t119, 0, 0, 0, 0, 0, t23, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t95, 0, 0, 0, 0, 0, t23, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t154, -t44, -t45, 0, t1 * qJD(2) - t112 + t128, t2 * qJD(2) - t111 - t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t154, -t44, -t45, 0, -t97 - t112, -t96 - t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
