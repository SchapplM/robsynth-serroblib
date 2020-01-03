% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:36
% EndTime: 2020-01-03 12:00:39
% DurationCPUTime: 0.73s
% Computational Cost: add. (1185->152), mult. (2465->206), div. (0->0), fcn. (2070->8), ass. (0->123)
t119 = -qJD(2) - qJD(4);
t111 = qJD(1) - t119;
t90 = sin(qJ(5));
t93 = cos(qJ(5));
t81 = -t90 ^ 2 + t93 ^ 2;
t149 = t111 * t81;
t148 = pkin(4) / 0.2e1;
t88 = sin(pkin(9));
t92 = sin(qJ(2));
t139 = t88 * t92;
t95 = cos(qJ(2));
t84 = t95 * pkin(1) + pkin(2);
t89 = cos(pkin(9));
t109 = -pkin(1) * t139 + t89 * t84;
t66 = pkin(3) + t109;
t94 = cos(qJ(4));
t58 = t94 * t66;
t138 = t89 * t92;
t70 = pkin(1) * t138 + t88 * t84;
t91 = sin(qJ(4));
t33 = t91 * t70 - t58;
t31 = -pkin(4) + t33;
t147 = -t31 / 0.2e1;
t141 = t88 * pkin(2);
t83 = t89 * pkin(2) + pkin(3);
t77 = t94 * t83;
t69 = t91 * t141 - t77;
t67 = -pkin(4) + t69;
t146 = -t67 / 0.2e1;
t145 = -t90 / 0.2e1;
t144 = t90 / 0.2e1;
t143 = -t93 / 0.2e1;
t142 = t93 / 0.2e1;
t73 = (t89 * t95 - t139) * pkin(1);
t136 = t91 * t73;
t72 = (-t88 * t95 - t138) * pkin(1);
t65 = t94 * t72;
t42 = -t65 + t136;
t140 = t42 * t93;
t137 = t91 * t72;
t135 = t94 * t73;
t34 = t91 * t66 + t94 * t70;
t30 = t34 * qJD(4);
t40 = t42 * qJD(2);
t134 = -t40 - t30;
t133 = pkin(1) * qJD(1);
t132 = pkin(1) * qJD(2);
t131 = pkin(4) * qJD(4);
t71 = t94 * t141 + t91 * t83;
t113 = t71 / 0.2e1 + t34 / 0.2e1;
t104 = -t42 / 0.2e1 + t113;
t7 = t104 * t93;
t130 = t7 * qJD(1);
t129 = qJD(1) * t31;
t128 = qJD(2) * t67;
t127 = qJD(5) * t90;
t87 = qJD(5) * t93;
t17 = t109 * t72 + t70 * t73;
t126 = t17 * qJD(1);
t125 = t33 * qJD(1);
t124 = t34 * qJD(1);
t123 = t42 * qJD(1);
t43 = t135 + t137;
t122 = t43 * qJD(1);
t121 = t71 * qJD(2);
t64 = t71 * qJD(4);
t120 = -qJD(1) - qJD(2);
t118 = t90 * t129;
t117 = t93 * t129;
t116 = t90 * t124;
t115 = t90 * t123;
t114 = -t58 / 0.2e1 - t77 / 0.2e1;
t112 = -t83 / 0.2e1 - t66 / 0.2e1;
t110 = pkin(1) * t120;
t108 = t141 / 0.2e1 + t70 / 0.2e1;
t107 = t33 / 0.2e1 + t148 + t147;
t106 = t69 / 0.2e1 + t148 + t146;
t105 = -t43 / 0.2e1 + t146 + t147;
t9 = t135 / 0.2e1 + (t72 / 0.2e1 + t108) * t91 + t114;
t103 = -t9 * qJD(1) - t69 * qJD(2);
t101 = t108 * t94;
t11 = -t65 / 0.2e1 - t101 + (t73 / 0.2e1 + t112) * t91;
t102 = -t11 * qJD(1) + t121;
t1 = t105 * t90;
t100 = t1 * qJD(1) - t90 * t128;
t2 = t105 * t93;
t99 = t2 * qJD(1) - t93 * t128;
t6 = t104 * t90;
t98 = -t6 * qJD(1) - t90 * t121;
t13 = t107 * t90;
t18 = t106 * t90;
t97 = t13 * qJD(1) + t18 * qJD(2) + t90 * t131;
t14 = t107 * t93;
t19 = t106 * t93;
t96 = t14 * qJD(1) + t19 * qJD(2) + t93 * t131;
t86 = pkin(4) * t143;
t85 = pkin(4) * t145;
t82 = t90 * t87;
t78 = t81 * qJD(5);
t68 = pkin(8) + t71;
t63 = t69 * qJD(4);
t57 = t90 * t64;
t50 = t111 * t93 * t90;
t49 = t67 * t142;
t48 = t67 * t144;
t41 = t43 * qJD(2);
t39 = t90 * t40;
t32 = pkin(8) + t34;
t29 = t33 * qJD(4);
t28 = t90 * t30;
t23 = t31 * t142;
t22 = t31 * t144;
t21 = t69 * t142 + t49 + t86;
t20 = t69 * t144 + t48 + t85;
t16 = t33 * t142 + t23 + t86;
t15 = t33 * t144 + t22 + t85;
t12 = -t136 / 0.2e1 + t65 / 0.2e1 - t101 + t112 * t91;
t10 = -t135 / 0.2e1 - t137 / 0.2e1 + t108 * t91 + t114;
t8 = -t140 / 0.2e1 - t113 * t93;
t5 = t113 * t90 + t42 * t144;
t4 = t43 * t143 + t23 + t49;
t3 = t43 * t145 + t22 + t48;
t24 = [0, 0, 0, 0, -t92 * t132, -t95 * t132, t17 * qJD(2), 0, t134, -t41 + t29, t82, t78, 0, 0, 0, t31 * t127 + t134 * t93, t31 * t87 + t28 + t39; 0, 0, 0, 0, t92 * t110, t95 * t110, t126 + (t72 * t89 + t73 * t88) * qJD(2) * pkin(2), 0, t12 * qJD(4) - t123 - t40, t10 * qJD(4) - t122 - t41, t82, t78, 0, 0, 0, t8 * qJD(4) + t3 * qJD(5) + t120 * t140, t5 * qJD(4) + t4 * qJD(5) + t115 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t12 * qJD(2) - t124 - t30, t10 * qJD(2) + t125 + t29, t82, t78, 0, 0, 0, t8 * qJD(2) + t15 * qJD(5) + (-qJD(1) - qJD(4)) * t93 * t34, t5 * qJD(2) + t16 * qJD(5) + t116 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t149, t87, -t127, 0, t3 * qJD(2) + t15 * qJD(4) - t32 * t87 + t118, t4 * qJD(2) + t16 * qJD(4) + t32 * t127 + t117; 0, 0, 0, 0, t92 * t133, t95 * t133, -t126, 0, t11 * qJD(4) + t123, t9 * qJD(4) + t122, t82, t78, 0, 0, 0, -t7 * qJD(4) - t1 * qJD(5) + t93 * t123, t6 * qJD(4) - t2 * qJD(5) - t115; 0, 0, 0, 0, 0, 0, 0, 0, -t64, t63, t82, t78, 0, 0, 0, t67 * t127 - t93 * t64, t67 * t87 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t102 - t64, -t103 + t63, t82, t78, 0, 0, 0, t119 * t93 * t71 + t20 * qJD(5) - t130, t21 * qJD(5) + t57 - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t149, t87, -t127, 0, t20 * qJD(4) - t68 * t87 - t100, t21 * qJD(4) + t68 * t127 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t87; 0, 0, 0, 0, 0, 0, 0, 0, -t11 * qJD(2) + t124, -t9 * qJD(2) - t125, t82, t78, 0, 0, 0, t7 * qJD(2) - t13 * qJD(5) + t93 * t124, -t6 * qJD(2) - t14 * qJD(5) - t116; 0, 0, 0, 0, 0, 0, 0, 0, t102, t103, t82, t78, 0, 0, 0, -t18 * qJD(5) + t93 * t121 + t130, -t19 * qJD(5) + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t78, 0, 0, 0, -pkin(4) * t127, -pkin(4) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t149, t87, -t127, 0, -pkin(8) * t87 - t97, pkin(8) * t127 - t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t149, 0, 0, 0, t1 * qJD(2) + t13 * qJD(4) - t118, t2 * qJD(2) + t14 * qJD(4) - t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t149, 0, 0, 0, t18 * qJD(4) + t100, t19 * qJD(4) + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t149, 0, 0, 0, t97, t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t24;
