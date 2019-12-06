% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRR8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:30
% EndTime: 2019-12-05 16:04:36
% DurationCPUTime: 1.20s
% Computational Cost: add. (552->156), mult. (1610->268), div. (0->0), fcn. (1511->8), ass. (0->134)
t90 = sin(qJ(4));
t84 = t90 ^ 2;
t178 = t84 / 0.2e1;
t89 = sin(qJ(5));
t177 = 0.2e1 * t89;
t92 = cos(qJ(5));
t93 = cos(qJ(4));
t159 = t92 * t93;
t115 = t159 * t177;
t83 = t89 ^ 2;
t85 = t92 ^ 2;
t69 = t85 - t83;
t97 = qJD(2) * t115 - qJD(4) * t69;
t87 = sin(pkin(5));
t94 = cos(qJ(2));
t169 = t87 * t94;
t88 = cos(pkin(5));
t45 = -t169 * t90 + t88 * t93;
t176 = t45 / 0.2e1;
t174 = -t90 / 0.2e1;
t173 = -t93 / 0.2e1;
t172 = t90 * pkin(8);
t171 = t93 * pkin(4);
t86 = t93 ^ 2;
t80 = t86 * t92;
t91 = sin(qJ(2));
t170 = t87 * t91;
t61 = t171 + t172;
t168 = t89 * t61;
t167 = t89 * t91;
t166 = t89 * t93;
t165 = t89 * t94;
t95 = -pkin(2) - pkin(7);
t164 = t89 * t95;
t163 = t90 * t91;
t162 = t91 * t92;
t161 = t92 * t61;
t160 = t92 * t84;
t158 = t92 * t94;
t157 = t92 * t95;
t156 = t93 * t95;
t68 = t84 - t86;
t134 = t89 * t156;
t112 = pkin(4) * t90 - pkin(8) * t93;
t108 = qJ(3) + t112;
t33 = -t108 * t92 + t164 * t90;
t9 = -t33 * t93 + (t134 + t161) * t90;
t155 = t9 * qJD(2);
t154 = qJD(2) * t87;
t153 = qJD(3) * t90;
t152 = qJD(4) * t92;
t151 = qJD(5) * t89;
t150 = qJD(5) * t92;
t133 = t92 * t156;
t135 = t90 * t157;
t34 = t108 * t89 + t135;
t10 = t34 * t93 + (-t133 + t168) * t90;
t149 = t10 * qJD(2);
t132 = 0.1e1 / 0.2e1 + t178;
t39 = (-t86 / 0.2e1 - t132) * t89;
t148 = t39 * qJD(2);
t40 = t80 / 0.2e1 + t132 * t92;
t147 = t40 * qJD(2);
t146 = t45 * qJD(4);
t56 = t68 * t89;
t145 = t56 * qJD(2);
t57 = -t80 + t160;
t144 = t57 * qJD(2);
t143 = t68 * qJD(2);
t142 = t90 * qJD(2);
t141 = t90 * qJD(4);
t140 = t93 * qJD(2);
t139 = t93 * qJD(4);
t138 = t93 * qJD(5);
t137 = qJ(3) * qJD(4);
t136 = qJD(2) * qJ(3);
t131 = t90 * t151;
t130 = t89 * t138;
t129 = t90 * t150;
t128 = t92 * t138;
t67 = t91 * t154;
t127 = t94 * t154;
t126 = t89 * t139;
t125 = t89 * t150;
t124 = t89 * t152;
t123 = t92 * t139;
t122 = t90 * t139;
t121 = t90 * t140;
t120 = -t166 / 0.2e1;
t119 = -t163 / 0.2e1;
t118 = t159 / 0.2e1;
t117 = t90 * t136;
t116 = t93 * t136;
t113 = qJD(4) * t115;
t17 = -t164 * t86 - t33 * t90;
t44 = t169 * t93 + t88 * t90;
t96 = t169 / 0.2e1 + t90 * t176 + t44 * t173;
t6 = t96 * t89;
t111 = -qJD(1) * t6 + qJD(2) * t17;
t18 = -t157 * t86 - t34 * t90;
t5 = t96 * t92;
t110 = qJD(1) * t5 - qJD(2) * t18;
t109 = (-qJD(5) - t142) * t93;
t107 = t172 / 0.2e1 + t171 / 0.2e1;
t106 = t162 * t87 - t45 * t89;
t105 = t167 * t87 + t45 * t92;
t100 = t61 / 0.2e1 + t107;
t23 = t100 * t89;
t104 = pkin(4) * t152 - qJD(2) * t23;
t24 = t100 * t92;
t103 = pkin(4) * qJD(4) * t89 + qJD(2) * t24;
t102 = t92 * t109;
t46 = (t83 / 0.2e1 - t85 / 0.2e1) * t93;
t101 = -qJD(2) * t46 + t124;
t99 = qJD(2) * t80 * t89 + qJD(4) * t46;
t55 = t69 * t86;
t98 = qJD(2) * t55 + t113;
t79 = t139 / 0.2e1;
t78 = t92 * t142;
t77 = t89 * t141;
t76 = t89 * t142;
t50 = (t142 + qJD(5) / 0.2e1) * t93;
t43 = t46 * qJD(5);
t42 = -t160 / 0.2e1 - t80 / 0.2e1 + t92 / 0.2e1;
t41 = (-0.1e1 / 0.2e1 + t178 + t86 / 0.2e1) * t89;
t16 = -t134 + t161 / 0.2e1 - t107 * t92;
t15 = -t133 - t168 / 0.2e1 + t107 * t89;
t14 = t44 * t92;
t12 = t44 * t89;
t8 = t105 * t174 + t44 * t118 + (t89 * t119 + t158 / 0.2e1) * t87;
t7 = t106 * t174 + t44 * t120 + (t92 * t119 - t165 / 0.2e1) * t87;
t4 = t105 * t173 + t118 * t45 + t120 * t170;
t3 = t106 * t93 / 0.2e1 + t166 * t176 + t118 * t170;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t67, -t127, t67, t127, ((-pkin(2) * t91 + qJ(3) * t94) * qJD(2) + t91 * qJD(3)) * t87, 0, 0, 0, 0, 0, (t139 * t91 + t142 * t94) * t87, (t140 * t94 - t141 * t91) * t87, 0, 0, 0, 0, 0, ((-t163 * t89 + t158) * t90 - t86 * t167) * t154 + t3 * qJD(4) + t8 * qJD(5), (-(t162 * t90 + t165) * t90 - t86 * t162) * t154 + t4 * qJD(4) + t7 * qJD(5); 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t93 - t146, qJD(4) * t44 - t67 * t90, 0, 0, 0, 0, 0, qJD(2) * t3 + qJD(5) * t12 - t146 * t92, qJD(2) * t4 + qJD(5) * t14 + t146 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(2) + t12 * qJD(4) - qJD(5) * t105, t7 * qJD(2) + t14 * qJD(4) - qJD(5) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * qJD(5), t6 * qJD(5); 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), -t122, t68 * qJD(4), 0, 0, 0, t137 * t93 + t153, qJD(3) * t93 - t137 * t90, -t122 * t85 - t125 * t86, -qJD(5) * t55 + t113 * t90, -qJD(4) * t57 - t130 * t90, qJD(4) * t56 - t128 * t90, t122, qJD(4) * t9 + qJD(5) * t18 + t153 * t92, -qJD(4) * t10 - qJD(5) * t17 - t153 * t89; 0, 0, 0, 0, 0, qJD(2), t136, 0, 0, 0, 0, 0, t142, t140, 0, 0, 0, 0, 0, qJD(5) * t42 + t78, qJD(5) * t41 - t76; 0, 0, 0, 0, 0, 0, 0, -t121, t143, -t141, -t139, 0, -t141 * t95 + t116, -t139 * t95 - t117, -t43 + (-t140 * t85 - t124) * t90, -0.2e1 * t93 * t125 + t90 * t97, t126 - t144, t123 + t145, t50, t155 + (t112 * t89 - t135) * qJD(4) + t16 * qJD(5), -t149 + (-pkin(8) * t159 + (pkin(4) * t92 + t164) * t90) * qJD(4) + t15 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t98, t89 * t109, t102, t79, qJD(3) * t42 + qJD(4) * t16 - qJD(5) * t34 - t110, qJD(3) * t41 + qJD(4) * t15 + qJD(5) * t33 - t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(2), -t136, 0, 0, 0, 0, 0, -t142, -t140, 0, 0, 0, 0, 0, -qJD(5) * t40 - t78, -qJD(5) * t39 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t139, 0, 0, 0, 0, 0, -t141 * t92 - t130, t77 - t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 - t129 - t147, -t123 + t131 - t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t121, -t143, 0, 0, 0, -t116, t117, t121 * t85 - t43, t102 * t177, t129 + t144, -t131 - t145, -t50, -qJD(5) * t24 - t155, qJD(5) * t23 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t69 * qJD(5), 0, 0, 0, -pkin(4) * t151, -pkin(4) * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t97, t78 + t150, -t76 - t151, -t140 / 0.2e1, -pkin(8) * t150 - t103, pkin(8) * t151 - t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(2), -t6 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t98, (t140 * t89 - t152) * t90, t121 * t92 + t77, t79, qJD(3) * t40 + qJD(4) * t24 + t110, qJD(3) * t39 - qJD(4) * t23 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, t97, -t78, t76, t140 / 0.2e1, t103, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
