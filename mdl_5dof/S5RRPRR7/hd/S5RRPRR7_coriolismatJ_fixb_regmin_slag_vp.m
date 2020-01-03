% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:43
% EndTime: 2019-12-31 20:15:46
% DurationCPUTime: 1.03s
% Computational Cost: add. (991->164), mult. (1937->175), div. (0->0), fcn. (1845->6), ass. (0->126)
t130 = qJD(4) + qJD(5);
t101 = sin(qJ(2));
t168 = t101 * pkin(1);
t86 = qJ(3) + t168;
t180 = -qJ(3) / 0.2e1 - t86 / 0.2e1;
t102 = cos(qJ(5));
t103 = cos(qJ(4));
t147 = t102 * t103;
t100 = sin(qJ(4));
t99 = sin(qJ(5));
t154 = t99 * t100;
t61 = -t147 + t154;
t148 = t102 * t100;
t153 = t99 * t103;
t62 = t148 + t153;
t20 = -t61 ^ 2 + t62 ^ 2;
t97 = qJD(1) + qJD(2);
t179 = t97 * t20;
t178 = t97 * t61;
t177 = t97 * t62;
t84 = t100 ^ 2 - t103 ^ 2;
t176 = t97 * t84;
t105 = -pkin(2) - pkin(7);
t95 = t100 * pkin(8);
t71 = t100 * t105 - t95;
t72 = (-pkin(8) + t105) * t103;
t175 = t130 * (-t102 * t71 - t99 * t72);
t174 = t130 * (-t102 * t72 + t99 * t71);
t104 = cos(qJ(2));
t167 = t104 * pkin(1);
t121 = -pkin(2) - t167;
t85 = -pkin(7) + t121;
t57 = t100 * t85 - t95;
t58 = (-pkin(8) + t85) * t103;
t173 = t130 * (-t102 * t57 - t99 * t58);
t172 = t130 * (-t102 * t58 + t99 * t57);
t169 = pkin(4) * t103;
t96 = t100 * pkin(4);
t73 = t86 + t96;
t166 = t73 * t61;
t165 = t73 * t62;
t87 = qJ(3) + t96;
t164 = t87 * t61;
t163 = t87 * t62;
t59 = t62 * qJD(3);
t157 = pkin(1) * qJD(2);
t92 = t104 * t157;
t162 = t62 * t92 + t59;
t60 = t61 * qJD(3);
t161 = -t61 * t92 - t60;
t93 = qJD(3) * t100;
t160 = t100 * t92 + t93;
t94 = qJD(3) * t103;
t159 = t103 * t92 + t94;
t158 = pkin(1) * qJD(1);
t156 = pkin(4) * qJD(5);
t155 = qJD(4) * pkin(4);
t152 = qJD(1) * t73;
t151 = qJD(2) * t87;
t150 = qJD(5) * t73;
t149 = qJD(5) * t87;
t50 = t62 * t169;
t18 = t50 - t166;
t138 = t18 * qJD(1);
t51 = t61 * t169;
t19 = -t51 - t165;
t137 = t19 * qJD(1);
t39 = -t121 * t168 - t86 * t167;
t136 = t39 * qJD(1);
t135 = t86 * qJD(1);
t134 = t92 + qJD(3);
t133 = qJ(3) * qJD(2);
t132 = t100 * qJD(4);
t131 = t103 * qJD(4);
t129 = t101 * t157;
t128 = t101 * t158;
t91 = t104 * t158;
t127 = t61 * t152;
t126 = t62 * t152;
t125 = t168 / 0.2e1;
t124 = t87 / 0.2e1 + t73 / 0.2e1;
t123 = t100 * t135;
t122 = t103 * t135;
t120 = t100 * t131;
t119 = pkin(4) * t130;
t27 = t130 * t61;
t118 = t61 * t91;
t117 = t62 * t91;
t75 = t97 * t103;
t116 = t100 * t91;
t115 = t103 * t91;
t23 = t50 - t164;
t106 = (-t154 / 0.2e1 + t147 / 0.2e1) * t168;
t8 = t124 * t61 + t106;
t4 = -t50 + t8;
t114 = t4 * qJD(1) - t23 * qJD(2);
t24 = -t51 - t163;
t107 = (-t148 / 0.2e1 - t153 / 0.2e1) * t168;
t9 = t124 * t62 + t107;
t5 = t51 + t9;
t113 = t5 * qJD(1) - t24 * qJD(2);
t112 = t125 + t180;
t111 = t8 * qJD(1) + t61 * t151;
t110 = t9 * qJD(1) + t62 * t151;
t35 = t112 * t100;
t109 = -t35 * qJD(1) + t100 * t133;
t36 = t112 * t103;
t108 = t36 * qJD(1) - t103 * t133;
t10 = -t163 / 0.2e1 - t165 / 0.2e1 + t107;
t11 = -t164 / 0.2e1 - t166 / 0.2e1 + t106;
t98 = qJ(3) * qJD(3);
t83 = t97 * qJ(3);
t82 = t86 * qJD(3);
t77 = t84 * qJD(4);
t74 = t97 * t100;
t64 = t97 * t168;
t52 = t100 * t75;
t38 = (t125 - t180) * t103;
t37 = t86 * t100;
t30 = t130 * t62;
t15 = t61 * t177;
t14 = t62 * t27;
t7 = -t51 + t10;
t6 = t50 + t11;
t1 = t130 * t20;
t2 = [0, 0, 0, 0, -t129, -t92, t129, t134, -t39 * qJD(2) + t82, -t120, t77, 0, 0, 0, t86 * t131 + t160, -t86 * t132 + t159, t14, t1, 0, 0, 0, t18 * qJD(4) - t61 * t150 + t162, t19 * qJD(4) - t62 * t150 + t161; 0, 0, 0, 0, -t64, -t91 - t92, t64, t91 + t134, -t136 + t82 + (-pkin(2) * t101 + qJ(3) * t104) * t157, -t120, t77, 0, 0, 0, t38 * qJD(4) + t116 + t160, -t37 * qJD(4) + t115 + t159, t14, t1, 0, 0, 0, t6 * qJD(4) + t11 * qJD(5) + t117 + t162, t7 * qJD(4) + t10 * qJD(5) - t118 + t161; 0, 0, 0, 0, 0, 0, 0, t97, t97 * t86, 0, 0, 0, 0, 0, t74, t75, 0, 0, 0, 0, 0, t177, -t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t176, -t132, -t131, 0, t38 * qJD(2) - t85 * t132 + t122, -t37 * qJD(2) - t85 * t131 - t123, t15, t179, -t30, t27, 0, t6 * qJD(2) + t138 + t173, t7 * qJD(2) + t137 + t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t179, -t30, t27, 0, t11 * qJD(2) - t127 + t173, t10 * qJD(2) - t126 + t172; 0, 0, 0, 0, t128, t91, -t128, -t91 + qJD(3), t98 + t136, -t120, t77, 0, 0, 0, -t36 * qJD(4) - t116 + t93, t35 * qJD(4) - t115 + t94, t14, t1, 0, 0, 0, -t4 * qJD(4) - t8 * qJD(5) - t117 + t59, -t5 * qJD(4) - t9 * qJD(5) + t118 - t60; 0, 0, 0, 0, 0, 0, 0, qJD(3), t98, -t120, t77, 0, 0, 0, qJ(3) * t131 + t93, -qJ(3) * t132 + t94, t14, t1, 0, 0, 0, t23 * qJD(4) - t61 * t149 + t59, t24 * qJD(4) - t62 * t149 - t60; 0, 0, 0, 0, 0, 0, 0, t97, t83, 0, 0, 0, 0, 0, t74, t75, 0, 0, 0, 0, 0, t177, -t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t176, -t132, -t131, 0, -t105 * t132 - t108, -t105 * t131 - t109, t15, t179, -t30, t27, 0, -t114 + t175, -t113 + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t179, -t30, t27, 0, -t111 + t175, -t110 + t174; 0, 0, 0, 0, 0, 0, 0, -t97, -t133 - t135, 0, 0, 0, 0, 0, -t74, -t75, 0, 0, 0, 0, 0, -t177, t178; 0, 0, 0, 0, 0, 0, 0, -t97, -t83, 0, 0, 0, 0, 0, -t74, -t75, 0, 0, 0, 0, 0, -t177, t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t131, 0, 0, 0, 0, 0, -t30, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t176, 0, 0, 0, t36 * qJD(2) - t122, -t35 * qJD(2) + t123, -t15, -t179, 0, 0, 0, t4 * qJD(2) - t138, t5 * qJD(2) - t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t176, 0, 0, 0, t108, t109, -t15, -t179, 0, 0, 0, t114, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99 * t156, -t102 * t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99 * t119, -t102 * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t179, 0, 0, 0, t8 * qJD(2) + t127, t9 * qJD(2) + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t179, 0, 0, 0, t111, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99 * t155, t102 * t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
