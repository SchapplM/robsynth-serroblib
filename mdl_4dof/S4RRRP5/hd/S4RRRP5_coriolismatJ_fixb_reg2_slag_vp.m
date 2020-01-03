% Calculate inertial parameters regressor of coriolis matrix for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:08
% EndTime: 2019-12-31 17:17:10
% DurationCPUTime: 0.99s
% Computational Cost: add. (1376->160), mult. (2799->185), div. (0->0), fcn. (2719->4), ass. (0->121)
t100 = cos(qJ(2));
t98 = sin(qJ(3));
t147 = t98 * t100;
t165 = cos(qJ(3));
t99 = sin(qJ(2));
t82 = t165 * t99 + t147;
t79 = t82 ^ 2;
t116 = t165 * t100;
t80 = t98 * t99 - t116;
t28 = t80 ^ 2 - t79;
t96 = qJD(2) + qJD(3);
t177 = t96 * t28;
t174 = t96 * t82;
t176 = t80 * t174;
t175 = t28 * qJD(1);
t169 = -pkin(6) - pkin(5);
t125 = t98 * t169;
t115 = t99 * t125;
t112 = t169 * t165;
t85 = t100 * t112;
t62 = -t85 + t115;
t114 = t169 * t147;
t84 = t99 * t112;
t63 = -t84 - t114;
t113 = t100 * t125;
t64 = t84 + t113;
t107 = t169 * t116;
t65 = t115 - t107;
t102 = (t62 - t65) * t82 + (-t63 - t64) * t80;
t173 = t102 * qJD(1);
t172 = t102 * qJD(2);
t171 = t65 / 0.2e1;
t167 = t98 * pkin(2);
t90 = qJ(4) + t167;
t170 = -t90 / 0.2e1;
t168 = t80 * pkin(3);
t166 = t99 * pkin(2);
t149 = t82 * qJ(4);
t110 = -t149 + t168;
t93 = -t100 * pkin(2) - pkin(1);
t36 = t110 + t93;
t164 = t36 * t80;
t163 = t36 * t82;
t162 = t82 * t98;
t161 = t90 * t82;
t124 = t165 * pkin(2);
t92 = -t124 - pkin(3);
t160 = t92 * t80;
t103 = -t107 / 0.2e1;
t30 = -t85 / 0.2e1 + t115 + t103;
t159 = -t62 * qJD(2) - t30 * qJD(3);
t158 = -t30 * qJD(2) - t65 * qJD(3);
t155 = pkin(2) * qJD(3);
t154 = qJD(2) * pkin(2);
t109 = t63 * t62 + t65 * t64;
t46 = t82 * pkin(3) + t80 * qJ(4);
t38 = t46 + t166;
t3 = t36 * t38 + t109;
t153 = t3 * qJD(1);
t4 = t36 * t46;
t152 = t4 * qJD(1);
t9 = t93 * t166 + t109;
t148 = t9 * qJD(1);
t146 = qJD(1) * t82;
t145 = qJD(1) * t93;
t144 = qJD(3) * t93;
t11 = (t170 + t167 / 0.2e1 + qJ(4) / 0.2e1) * t82 + (-t124 / 0.2e1 - t92 / 0.2e1 - pkin(3) / 0.2e1) * t80;
t143 = t11 * qJD(1);
t13 = t38 * t80 + t163;
t142 = t13 * qJD(1);
t14 = -t38 * t82 + t164;
t141 = t14 * qJD(1);
t15 = t46 * t80 + t163;
t140 = t15 * qJD(1);
t16 = -t46 * t82 + t164;
t139 = t16 * qJD(1);
t32 = t80 * t166 + t93 * t82;
t136 = t32 * qJD(1);
t33 = t82 * t166 - t93 * t80;
t135 = t33 * qJD(1);
t61 = t85 / 0.2e1 + t103;
t134 = t61 * qJD(1);
t49 = t61 * qJD(2);
t133 = t79 * qJD(1);
t74 = t80 * qJD(4);
t88 = t100 ^ 2 - t99 ^ 2;
t132 = t88 * qJD(1);
t131 = t99 * qJD(2);
t95 = qJD(3) * t124;
t130 = t95 + qJD(4);
t129 = qJD(1) * t100;
t128 = t100 * qJD(2);
t127 = pkin(1) * t99 * qJD(1);
t126 = t98 * t155;
t123 = pkin(1) * t129;
t122 = t36 * t146;
t121 = t80 * t146;
t120 = t80 * t145;
t119 = t82 * t145;
t118 = t165 * t80;
t117 = t99 * t128;
t101 = (t165 * t171 + t63 * t98 / 0.2e1) * pkin(2) + t63 * t170 + t92 * t171;
t104 = t62 * pkin(3) / 0.2e1 - t64 * qJ(4) / 0.2e1;
t2 = t101 + t104;
t66 = (t165 * t90 + t92 * t98) * pkin(2);
t108 = t2 * qJD(1) + t66 * qJD(2);
t29 = -t84 - t113 / 0.2e1 - t114 / 0.2e1;
t106 = -t29 * qJD(2) - t63 * qJD(3);
t105 = t64 * qJD(2) - t29 * qJD(3);
t97 = qJ(4) * qJD(4);
t94 = qJD(2) * t124;
t89 = t99 * t129;
t87 = t96 * qJ(4);
t86 = t90 * qJD(4);
t51 = t61 * qJD(3);
t43 = t96 * t80;
t35 = t98 * t154 + t134;
t26 = -t96 * t167 - t134;
t10 = -t161 / 0.2e1 - t160 / 0.2e1 - t149 / 0.2e1 + t168 / 0.2e1 + (-t118 / 0.2e1 + t162 / 0.2e1) * pkin(2);
t1 = t101 - t104;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, t88 * qJD(2), 0, -t117, 0, 0, -pkin(1) * t131, -pkin(1) * t128, 0, 0, -t176, t177, 0, t176, 0, 0, t32 * qJD(2) + t82 * t144, t33 * qJD(2) - t80 * t144, t172, t9 * qJD(2), -t176, 0, -t177, 0, 0, t176, t13 * qJD(2) + t15 * qJD(3) - t82 * t74, t172, t14 * qJD(2) + t16 * qJD(3) + t79 * qJD(4), t3 * qJD(2) + t4 * qJD(3) - qJD(4) * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t132, t128, -t89, -t131, 0, -pkin(5) * t128 - t127, pkin(5) * t131 - t123, 0, 0, -t121, t175, -t43, t121, -t174, 0, t136 + t159, -t105 + t135, t173 + (t118 - t162) * t154, t148 + (-t165 * t62 + t64 * t98) * t154, -t121, -t43, -t175, 0, t174, t121, t142 + t159, t173 + (-t160 - t161) * qJD(2) + t10 * qJD(3) - t74, t105 + t141, t153 + (t62 * t92 + t64 * t90) * qJD(2) + t1 * qJD(3) + t30 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t175, -t43, t121, -t174, 0, t119 + t158, -t106 - t120, 0, 0, -t121, -t43, -t175, 0, t174, t121, t140 + t158, t10 * qJD(2) + t110 * qJD(3) - t74, t106 + t139, t152 + t1 * qJD(2) + (-t65 * pkin(3) - t63 * qJ(4)) * qJD(3) + t65 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t43, t133, -t122 - t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t132, 0, t89, 0, 0, t127, t123, 0, 0, t121, -t175, 0, -t121, 0, 0, -t51 - t136, -t135, -t173, -t148, t121, 0, t175, 0, 0, -t121, -t51 - t142, t11 * qJD(3) - t173, -t141, t2 * qJD(3) + t61 * qJD(4) - t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t95, 0, 0, 0, 0, 0, 0, 0, 0, -t126, 0, t130, t66 * qJD(3) + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t94 - t95, 0, 0, 0, 0, 0, 0, 0, 0, t26, t143, t130 + t94, (-pkin(3) * t98 + t165 * qJ(4)) * t155 + t86 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t96 * t90 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, -t175, 0, -t121, 0, 0, t49 - t119, t120, 0, 0, t121, 0, t175, 0, 0, -t121, t49 - t140, -t11 * qJD(2), -t139, -t2 * qJD(2) - t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t94, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t143, qJD(4) - t94, -t108 + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, -t133, t122 - t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -qJ(4) * qJD(3) - t90 * qJD(2) - t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
