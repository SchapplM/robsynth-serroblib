% Calculate inertial parameters regressor of coriolis matrix for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR7_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:08
% EndTime: 2019-12-31 16:54:10
% DurationCPUTime: 1.24s
% Computational Cost: add. (2018->150), mult. (4217->231), div. (0->0), fcn. (4444->6), ass. (0->139)
t116 = cos(qJ(4));
t114 = sin(qJ(4));
t112 = sin(pkin(7));
t194 = cos(qJ(3));
t140 = t194 * t112;
t113 = cos(pkin(7));
t115 = sin(qJ(3));
t177 = t115 * t113;
t98 = t140 + t177;
t55 = t114 * t98;
t136 = 0.2e1 * t116 * t55;
t96 = t112 * t115 - t113 * t194;
t93 = t96 ^ 2;
t94 = t98 ^ 2;
t69 = t93 + t94;
t148 = t94 - t93;
t197 = t96 * pkin(6);
t196 = t98 * pkin(3);
t111 = t116 ^ 2;
t195 = -t111 / 0.2e1;
t107 = -pkin(2) * t113 - pkin(1);
t132 = pkin(3) * t96 - pkin(6) * t98;
t118 = t107 + t132;
t190 = pkin(5) + qJ(2);
t101 = t190 * t113;
t137 = t190 * t112;
t71 = t101 * t194 - t115 * t137;
t187 = t114 * t71;
t31 = -t116 * t118 + t187;
t193 = t31 * t96;
t185 = t116 * t71;
t32 = t114 * t118 + t185;
t192 = t32 * t96;
t70 = t101 * t115 + t137 * t194;
t191 = t70 * t98;
t179 = t70 * t114;
t68 = t196 + t197;
t186 = t116 * t68;
t35 = t179 + t186;
t178 = t70 * t116;
t188 = t114 * t68;
t36 = -t178 + t188;
t1 = (t35 * t98 + t193) * t116 + (t36 * t98 - t192) * t114;
t189 = t1 * qJD(1);
t2 = -t31 * t35 + t32 * t36 + t70 * t71;
t184 = t2 * qJD(1);
t130 = -t70 * t96 + t71 * t98;
t4 = t114 * t130 - t31 * t98 + t35 * t96;
t183 = t4 * qJD(1);
t5 = t116 * t130 - t32 * t98 - t36 * t96;
t182 = t5 * qJD(1);
t6 = t191 + (-t114 * t31 - t116 * t32) * t96;
t181 = t6 * qJD(1);
t110 = t114 ^ 2;
t117 = (t195 - t110 / 0.2e1) * t197 - t196 / 0.2e1;
t126 = t35 * t116 / 0.2e1 + t36 * t114 / 0.2e1;
t7 = t117 - t126;
t180 = t7 * qJD(1);
t21 = -t179 * t98 + t193;
t176 = t21 * qJD(1);
t22 = t178 * t98 - t192;
t175 = t22 * qJD(1);
t24 = -t71 * t96 + t191;
t174 = t24 * qJD(1);
t39 = t148 * t114;
t173 = t39 * qJD(1);
t40 = t69 * t114;
t172 = t40 * qJD(1);
t41 = t148 * t116;
t171 = t41 * qJD(1);
t170 = t148 * qJD(1);
t52 = (t110 / 0.2e1 + t195) * t98;
t169 = t52 * qJD(4);
t53 = t114 * t96;
t168 = t53 * qJD(1);
t167 = t55 * qJD(1);
t58 = t116 * t96;
t166 = t58 * qJD(1);
t88 = t110 * t96;
t89 = t111 * t96;
t61 = t88 + t89;
t165 = t61 * qJD(1);
t63 = t69 * t116;
t164 = t63 * qJD(1);
t163 = t69 * qJD(1);
t92 = t140 / 0.2e1 + t177 / 0.2e1;
t162 = t92 * qJD(1);
t161 = t96 * qJD(1);
t160 = t98 * qJD(1);
t159 = t98 * qJD(3);
t102 = t112 ^ 2 + t113 ^ 2;
t105 = t111 - t110;
t158 = qJD(1) * t107;
t157 = qJD(1) * t116;
t156 = qJD(2) * t116;
t155 = qJD(3) * t107;
t154 = qJD(3) * t114;
t153 = qJD(3) * t116;
t152 = qJD(4) * t114;
t151 = qJD(4) * t116;
t100 = t102 * qJ(2);
t150 = t100 * qJD(1);
t149 = t102 * qJD(1);
t146 = t96 * t160;
t145 = t96 * t159;
t144 = t110 * t160;
t143 = t111 * t160;
t142 = t96 * t151;
t141 = t98 * t157;
t139 = t114 * t151;
t138 = t114 * t153;
t135 = -qJD(4) - t161;
t134 = qJD(2) + t158;
t133 = t94 * t139;
t131 = qJD(3) * t136;
t129 = -t114 * t35 + t116 * t36;
t128 = t135 * t116;
t127 = t197 / 0.2e1 + t196 / 0.2e1;
t121 = t68 / 0.2e1 + t127;
t19 = t121 * t116;
t125 = pkin(3) * t154 + qJD(1) * t19;
t17 = t121 * t114;
t124 = pkin(3) * t153 - qJD(1) * t17;
t123 = qJD(4) * t92 + t146;
t122 = t98 * t128;
t43 = -qJD(1) * t52 + t138;
t37 = t114 * t157 * t94 + qJD(3) * t52;
t62 = t105 * t94;
t120 = qJD(1) * t62 + t131;
t119 = qJD(1) * t136 - qJD(3) * t105;
t91 = t96 * qJD(3);
t87 = t92 * qJD(3);
t86 = t98 * t153;
t49 = t53 * qJD(4);
t44 = -t152 - t168;
t20 = t179 + t186 / 0.2e1 - t127 * t116;
t18 = t178 - t188 / 0.2e1 + t127 * t114;
t8 = t117 + t126;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * qJD(2), t100 * qJD(2), -t145, -t148 * qJD(3), 0, t145, 0, 0, t98 * t155, -t96 * t155, t69 * qJD(2), t24 * qJD(2), -t111 * t145 - t133, -qJD(4) * t62 + t131 * t96, -t152 * t96 * t98 + qJD(3) * t41, -t110 * t145 + t133, -qJD(3) * t39 - t142 * t98, t145, qJD(2) * t40 + qJD(3) * t4 + qJD(4) * t22, qJD(2) * t63 + qJD(3) * t5 + qJD(4) * t21, -t1 * qJD(3), qJD(2) * t6 + qJD(3) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t150, 0, 0, 0, 0, 0, 0, 0, 0, t163, t174, 0, 0, 0, 0, 0, 0, t172, t164, 0, qJD(3) * t8 + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, -t170, -t91, t146, -t159, 0, -qJD(3) * t71 + t158 * t98, qJD(3) * t70 - t158 * t96, 0, 0, -t169 + (-t138 - t143) * t96, (t88 - t89) * qJD(3) + (-qJD(4) + t161) * t136, t154 * t98 + t171, t169 + (t138 - t144) * t96, t86 - t173, t123, t183 + (t114 * t132 - t185) * qJD(3) + t20 * qJD(4), t182 + (t116 * t132 + t187) * qJD(3) + t18 * qJD(4), qJD(3) * t129 - t189, t184 + t8 * qJD(2) + (-t71 * pkin(3) + pkin(6) * t129) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t120, t135 * t55, t37, t122, t87, qJD(3) * t20 - qJD(4) * t32 + t175, qJD(3) * t18 + qJD(4) * t31 + t176, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, -t150, 0, 0, 0, 0, 0, 0, t159, -t91, -t163, -t174, 0, 0, 0, 0, 0, 0, -t49 + t86 - t172, -qJD(3) * t55 - t142 - t164, t61 * qJD(3), -qJD(3) * t7 - t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, -t161, 0, 0, 0, 0, 0, 0, 0, 0, t141, -t167, t165, -t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t128, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t170, 0, -t146, 0, 0, -t134 * t98, t134 * t96, 0, 0, t143 * t96 - t169, 0.2e1 * t114 * t122, qJD(4) * t58 - t171, t144 * t96 + t169, -t49 + t173, -t123, -qJD(4) * t19 - t156 * t98 - t183, qJD(2) * t55 + qJD(4) * t17 - t182, -qJD(2) * t61 + t189, qJD(2) * t7 - t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, t161, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t167, -t165, t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t105 * qJD(4), 0, -t139, 0, 0, -pkin(3) * t152, -pkin(3) * t151, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t119, t151 + t166, -t43, t44, -t162, -pkin(6) * t151 - t125, pkin(6) * t152 - t124, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t120, -qJD(3) * t58 + t114 * t146, -t37, qJD(3) * t53 + t141 * t96, t87, qJD(2) * t53 + qJD(3) * t19 - t175, -qJD(3) * t17 + t156 * t96 - t176, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t96 * t157, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t119, -t166, t43, t168, t162, t125, t124, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
