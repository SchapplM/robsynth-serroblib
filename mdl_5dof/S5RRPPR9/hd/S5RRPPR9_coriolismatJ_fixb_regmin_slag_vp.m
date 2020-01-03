% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:50
% EndTime: 2019-12-31 19:41:54
% DurationCPUTime: 1.30s
% Computational Cost: add. (644->166), mult. (1260->233), div. (0->0), fcn. (1010->4), ass. (0->141)
t84 = cos(qJ(2));
t122 = t84 * qJD(3);
t82 = sin(qJ(2));
t69 = t82 * qJ(3);
t152 = -t84 * pkin(2) - t69;
t164 = t152 * qJD(2) + t122;
t163 = pkin(6) - qJ(4);
t83 = cos(qJ(5));
t153 = t83 * t84;
t81 = sin(qJ(5));
t120 = t81 * t153;
t100 = 0.2e1 * t120;
t76 = t81 ^ 2;
t78 = t83 ^ 2;
t51 = t78 - t76;
t89 = qJD(1) * t100 - qJD(2) * t51;
t85 = -pkin(2) - pkin(3);
t47 = t163 * t82;
t162 = t47 * t81;
t49 = t163 * t84;
t160 = t49 * t84;
t73 = -pkin(7) + t85;
t159 = t73 * t82;
t77 = t82 ^ 2;
t158 = t77 * t81;
t80 = qJ(3) + pkin(4);
t157 = t80 * t84;
t21 = t157 + t159;
t156 = t81 * t21;
t79 = t84 ^ 2;
t61 = t81 * t79;
t155 = t83 * t21;
t154 = t83 * t47;
t50 = t77 + t79;
t52 = t79 - t77;
t44 = -pkin(1) + t152;
t28 = t84 * pkin(3) - t44;
t87 = pkin(4) * t82 + pkin(7) * t84 + t28;
t10 = -t83 * t87 + t162;
t4 = -t10 * t82 - t81 * t160;
t151 = qJD(1) * t4;
t11 = t81 * t87 + t154;
t5 = -t11 * t82 - t49 * t153;
t150 = qJD(1) * t5;
t1 = t82 * t155 + (-t10 + t162) * t84;
t149 = t1 * qJD(1);
t2 = t156 * t82 + (t11 - t154) * t84;
t148 = t2 * qJD(1);
t146 = t84 * qJ(3);
t30 = t85 * t82 + t146;
t3 = t28 * t30;
t147 = t3 * qJD(1);
t12 = t28 * t84 + t30 * t82;
t145 = qJD(1) * t12;
t13 = -t28 * t82 + t30 * t84;
t144 = qJD(1) * t13;
t18 = t47 * t82 + t160;
t143 = qJD(1) * t18;
t38 = -t61 + t158;
t142 = qJD(1) * t38;
t40 = t52 * t83;
t141 = qJD(1) * t40;
t140 = qJD(3) * t81;
t139 = qJD(3) * t82;
t138 = qJD(3) * t83;
t137 = qJD(4) * t84;
t136 = qJD(5) * t81;
t135 = qJD(5) * t82;
t134 = qJD(5) * t83;
t48 = t82 * pkin(2) - t146;
t14 = t44 * t84 + t48 * t82;
t133 = t14 * qJD(1);
t15 = -t44 * t82 + t48 * t84;
t132 = t15 * qJD(1);
t119 = -pkin(2) / 0.2e1 - pkin(3) / 0.2e1;
t22 = t146 + (t85 / 0.2e1 + t119) * t82;
t131 = t22 * qJD(1);
t37 = t61 + t158;
t130 = t37 * qJD(1);
t39 = t50 * t83;
t129 = t39 * qJD(1);
t128 = t50 * qJD(1);
t127 = t52 * qJD(1);
t65 = t77 * qJD(1);
t64 = t77 * qJD(3);
t126 = t81 * qJD(2);
t66 = t82 * qJD(1);
t125 = t82 * qJD(2);
t124 = t83 * qJD(2);
t123 = t84 * qJD(1);
t67 = t84 * qJD(2);
t121 = t84 * qJD(5);
t118 = pkin(1) * t66;
t117 = pkin(1) * t123;
t116 = pkin(6) * t125;
t115 = pkin(6) * t67;
t114 = t81 * t65;
t113 = t81 * t124;
t112 = t81 * t67;
t111 = t83 * t67;
t110 = t82 * t122;
t109 = t82 * t134;
t108 = t82 * t121;
t107 = t44 * t48 * qJD(1);
t106 = t44 * t66;
t105 = t81 * t134;
t104 = t81 * t66;
t103 = t81 * t123;
t56 = t82 * t67;
t55 = t82 * t123;
t102 = t83 * t123;
t101 = -0.2e1 * t120;
t99 = qJD(5) + t66;
t98 = qJD(2) * t100;
t96 = t99 * t84;
t95 = t159 / 0.2e1 + t157 / 0.2e1;
t88 = t21 / 0.2e1 + t95;
t8 = t88 * t81;
t94 = -qJD(1) * t8 + t80 * t124;
t9 = t88 * t83;
t93 = qJD(1) * t9 + t80 * t126;
t29 = (t76 / 0.2e1 - t78 / 0.2e1) * t84;
t92 = -qJD(1) * t29 + t113;
t91 = qJD(1) * t83 * t61 + qJD(2) * t29;
t36 = t51 * t79;
t90 = qJD(1) * t36 + t98;
t86 = (-t73 * t84 + t80 * t82) * qJD(2) - t122;
t75 = qJ(3) * qJD(3);
t74 = qJD(2) * qJ(3);
t60 = t67 / 0.2e1;
t58 = t81 * t135;
t57 = t83 * t66;
t41 = t49 * qJD(2);
t34 = -t57 - t134;
t33 = t99 * t81;
t32 = t55 + t121 / 0.2e1;
t25 = t29 * qJD(5);
t24 = (-t85 / 0.2e1 + t119) * t82;
t7 = -t49 * t81 + t155 / 0.2e1 - t95 * t83;
t6 = -t49 * t83 - t156 / 0.2e1 + t95 * t81;
t16 = [0, 0, 0, t56, t52 * qJD(2), 0, 0, 0, -pkin(1) * t125, -pkin(1) * t67, -qJD(2) * t15 + t110, 0, -qJD(2) * t14 + t64, (qJD(2) * t48 - t139) * t44, qJD(2) * t12 + t64, -qJD(2) * t13 - t110, qJD(4) * t50, qJD(2) * t3 - qJD(4) * t18 + t28 * t139, -t105 * t79 - t56 * t78, -qJD(5) * t36 + t82 * t98, -qJD(2) * t40 + t108 * t81, -qJD(2) * t38 + t108 * t83, t56, qJD(2) * t1 + qJD(4) * t37 + qJD(5) * t5 + t83 * t64, -qJD(2) * t2 + qJD(4) * t39 - qJD(5) * t4 - t81 * t64; 0, 0, 0, t55, t127, t67, -t125, 0, -t115 - t118, t116 - t117, -t115 - t132, t164, -t116 - t133, t164 * pkin(6) + t107, -qJD(2) * t47 + t145, t41 - t144, (-t84 * t85 + t69) * qJD(2) - t122, t147 + (-t47 * qJ(3) + t49 * t85) * qJD(2) + t49 * qJD(3) + t24 * qJD(4), -t25 + (-t78 * t123 - t113) * t82, qJD(5) * t101 + t89 * t82, -t112 - t141, -t111 - t142, t32, t7 * qJD(5) - t47 * t124 + t81 * t86 + t149, t6 * qJD(5) + t47 * t126 + t83 * t86 - t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t67, t65, -t106 + t115, t65, -t55, -t67, t28 * t66 + t41, 0, 0, 0, 0, 0, t83 * t65 - t112, -t111 - t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, qJD(2) * t24 - t143, 0, 0, 0, 0, 0, t130, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t90, t81 * t96, t83 * t96, t60, qJD(2) * t7 - qJD(5) * t11 + t150, qJD(2) * t6 + qJD(5) * t10 - t151; 0, 0, 0, -t55, -t127, 0, 0, 0, t118, t117, t132, 0, t133, -t107, -t137 - t145, -qJD(4) * t82 + t144, 0, -qJD(4) * t22 - t147, t55 * t78 - t25, t99 * t101, -t109 + t141, t58 + t142, -t32, -qJD(5) * t9 - t83 * t137 - t149, qJD(5) * t8 + t81 * t137 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t75, qJD(3), 0, 0, t75, t105, t51 * qJD(5), 0, 0, 0, -t80 * t136 + t138, -t80 * t134 - t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t74, qJD(2), 0, 0, t74, 0, 0, 0, 0, 0, t124, -t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t66, 0, -t131, 0, 0, 0, 0, 0, -t102, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, -t89, t34, t33, -t123 / 0.2e1, -t73 * t134 - t93, t73 * t136 - t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, 0, -t65, t106, -t65, t55, 0, (-qJD(1) * t28 - qJD(4)) * t82, 0, 0, 0, 0, 0, (-t65 - t135) * t83, t58 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t74, -qJD(2), 0, 0, -t74, 0, 0, 0, 0, 0, -t124, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t125, -t128, qJD(2) * t22 + t139 + t143, 0, 0, 0, 0, 0, t111 - t58 - t130, -t109 - t112 - t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t66, 0, t131, 0, 0, 0, 0, 0, t102, -t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t90, (-t103 + t124) * t82, (-t102 - t126) * t82, t60, -t150 + qJD(2) * t9 + (qJD(4) * t81 + t138) * t82, t151 - qJD(2) * t8 + (qJD(4) * t83 - t140) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t89, t57, -t104, t123 / 0.2e1, t93, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t16;
