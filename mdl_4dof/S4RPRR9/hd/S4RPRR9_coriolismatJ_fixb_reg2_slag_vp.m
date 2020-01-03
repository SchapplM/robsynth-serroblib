% Calculate inertial parameters regressor of coriolis matrix for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR9_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:31
% EndTime: 2019-12-31 16:56:33
% DurationCPUTime: 1.13s
% Computational Cost: add. (839->149), mult. (1741->225), div. (0->0), fcn. (1384->4), ass. (0->124)
t90 = sin(qJ(3));
t86 = t90 ^ 2;
t170 = t86 / 0.2e1;
t89 = sin(qJ(4));
t169 = 0.2e1 * t89;
t91 = cos(qJ(4));
t92 = cos(qJ(3));
t154 = t91 * t92;
t110 = t154 * t169;
t85 = t89 ^ 2;
t87 = t91 ^ 2;
t67 = t87 - t85;
t94 = qJD(1) * t110 - qJD(3) * t67;
t93 = -pkin(1) - pkin(5);
t151 = t92 * t93;
t166 = t90 * pkin(3);
t106 = -pkin(6) * t92 + t166;
t101 = qJ(2) + t106;
t153 = t91 * t93;
t124 = t90 * t153;
t40 = t101 * t89 + t124;
t160 = t40 * t92;
t157 = t90 * t93;
t39 = -t101 * t91 + t157 * t89;
t162 = t39 * t92;
t123 = t89 * t151;
t164 = t92 * pkin(3);
t165 = t90 * pkin(6);
t59 = t164 + t165;
t156 = t91 * t59;
t41 = -t123 + t156;
t122 = t91 * t151;
t159 = t89 * t59;
t42 = t122 + t159;
t2 = -t90 * t151 + (t42 * t90 / 0.2e1 + t160 / 0.2e1) * t91 + (-t41 * t90 / 0.2e1 + t162 / 0.2e1) * t89;
t167 = t2 * qJD(3);
t163 = t39 * t91;
t161 = t40 * t90;
t88 = t92 ^ 2;
t80 = t88 * t91;
t158 = t89 * t93;
t155 = t91 * t86;
t152 = t92 * t90;
t150 = t85 + t87;
t66 = t86 - t88;
t3 = t90 * t163 - t89 * t161 + (t41 * t91 + t42 * t89) * t92;
t149 = t3 * qJD(1);
t6 = -t162 + (t41 + 0.2e1 * t123) * t90;
t148 = t6 * qJD(1);
t7 = t160 + (t42 - 0.2e1 * t122) * t90;
t147 = t7 * qJD(1);
t144 = qJD(2) * t90;
t143 = qJD(3) * t91;
t142 = qJD(4) * t89;
t141 = qJD(4) * t91;
t15 = t40 * t89 - t163;
t140 = t15 * qJD(1);
t18 = -t158 * t88 - t39 * t90;
t139 = t18 * qJD(1);
t19 = -t153 * t88 - t161;
t138 = t19 * qJD(1);
t121 = 0.1e1 / 0.2e1 + t170;
t43 = (-t88 / 0.2e1 - t121) * t89;
t137 = t43 * qJD(1);
t44 = t80 / 0.2e1 + t121 * t91;
t136 = t44 * qJD(1);
t49 = (t85 / 0.2e1 - t87 / 0.2e1) * t92;
t135 = t49 * qJD(4);
t55 = t150 * t92;
t134 = t55 * qJD(1);
t56 = t66 * t89;
t133 = t56 * qJD(1);
t57 = -t80 + t155;
t132 = t57 * qJD(1);
t131 = t66 * qJD(1);
t130 = t90 * qJD(1);
t129 = t90 * qJD(3);
t128 = t92 * qJD(1);
t127 = t92 * qJD(3);
t126 = t92 * qJD(4);
t125 = qJ(2) * qJD(3);
t83 = qJD(1) * qJ(2);
t120 = t90 * t142;
t119 = t89 * t126;
t118 = t90 * t141;
t117 = t91 * t126;
t116 = t89 * t127;
t115 = t89 * t141;
t114 = t89 * t143;
t113 = t91 * t127;
t78 = t90 * t127;
t77 = t90 * t128;
t112 = t90 * t83;
t111 = t92 * t83;
t109 = t88 * t115;
t107 = qJD(3) * t110;
t105 = -t41 * t89 + t42 * t91;
t4 = -t152 * t93 ^ 2 - t39 * t41 + t40 * t42;
t104 = t4 * qJD(1) + t2 * qJD(2);
t38 = (-0.1e1 + t150) * t152;
t103 = -t2 * qJD(1) - t38 * qJD(2);
t102 = (-qJD(4) - t130) * t92;
t100 = t165 / 0.2e1 + t164 / 0.2e1;
t96 = t59 / 0.2e1 + t100;
t24 = t96 * t89;
t99 = pkin(3) * t143 - qJD(1) * t24;
t25 = t96 * t91;
t98 = pkin(3) * qJD(3) * t89 + qJD(1) * t25;
t97 = t91 * t102;
t37 = -qJD(1) * t49 + t114;
t34 = qJD(1) * t80 * t89 + qJD(3) * t49;
t54 = t67 * t88;
t95 = qJD(1) * t54 + t107;
t84 = qJ(2) * qJD(2);
t79 = t127 / 0.2e1;
t76 = t91 * t130;
t75 = t89 * t129;
t74 = t89 * t130;
t53 = t77 + t126 / 0.2e1;
t46 = -t155 / 0.2e1 - t80 / 0.2e1 + t91 / 0.2e1;
t45 = (-0.1e1 / 0.2e1 + t170 + t88 / 0.2e1) * t89;
t17 = -t123 + t156 / 0.2e1 - t100 * t91;
t16 = -t122 - t159 / 0.2e1 + t100 * t89;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t84, -t78, t66 * qJD(3), 0, t78, 0, 0, t125 * t92 + t144, qJD(2) * t92 - t125 * t90, 0, t84, -t78 * t87 - t109, -qJD(4) * t54 + t107 * t90, -qJD(3) * t57 - t119 * t90, -t78 * t85 + t109, qJD(3) * t56 - t117 * t90, t78, qJD(3) * t6 + qJD(4) * t19 + t144 * t91, -qJD(3) * t7 - qJD(4) * t18 - t144 * t89, -qJD(2) * t55 - qJD(3) * t3, qJD(2) * t15 + qJD(3) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t83, 0, 0, 0, 0, 0, 0, t130, t128, 0, t83, 0, 0, 0, 0, 0, 0, qJD(4) * t46 + t76, qJD(4) * t45 - t74, -t134, t140 + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t131, -t129, t77, -t127, 0, -t129 * t93 + t111, -t127 * t93 - t112, 0, 0, -t135 + (-t128 * t87 - t114) * t90, -0.2e1 * t92 * t115 + t90 * t94, t116 - t132, t135 + (-t128 * t85 + t114) * t90, t113 + t133, t53, t148 + (t106 * t89 - t124) * qJD(3) + t17 * qJD(4), -t147 + (-pkin(6) * t154 + (pkin(3) * t91 + t158) * t90) * qJD(3) + t16 * qJD(4), qJD(3) * t105 - t149, (-pkin(3) * t157 + pkin(6) * t105) * qJD(3) + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t95, t89 * t102, t34, t97, t79, qJD(2) * t46 + qJD(3) * t17 - qJD(4) * t40 + t138, qJD(2) * t45 + qJD(3) * t16 + qJD(4) * t39 - t139, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t83, 0, 0, 0, 0, 0, 0, -t130, -t128, 0, -t83, 0, 0, 0, 0, 0, 0, -qJD(4) * t44 - t76, -qJD(4) * t43 + t74, t134, -t140 + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -t127, 0, 0, 0, 0, 0, 0, 0, 0, -t129 * t91 - t119, t75 - t117, t55 * qJD(3), (pkin(6) * t55 - t166) * qJD(3) - t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 - t118 - t136, -t113 + t120 - t137, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t131, 0, -t77, 0, 0, -t111, t112, 0, 0, t77 * t87 - t135, t97 * t169, t118 + t132, t77 * t85 + t135, -t120 - t133, -t53, -qJD(4) * t25 - t148, qJD(4) * t24 + t147, t149, -t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t67 * qJD(4), 0, -t115, 0, 0, -pkin(3) * t142, -pkin(3) * t141, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t94, t76 + t141, -t37, -t74 - t142, -t128 / 0.2e1, -pkin(6) * t141 - t98, pkin(6) * t142 - t99, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t95, (t128 * t89 - t143) * t90, -t34, t77 * t91 + t75, t79, qJD(2) * t44 + qJD(3) * t25 - t138, qJD(2) * t43 - qJD(3) * t24 + t139, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t137, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t94, -t76, t37, t74, t128 / 0.2e1, t98, t99, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
