% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR12_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:24
% EndTime: 2019-12-31 18:07:28
% DurationCPUTime: 1.18s
% Computational Cost: add. (1226->145), mult. (2528->222), div. (0->0), fcn. (2735->6), ass. (0->133)
t111 = sin(pkin(8));
t115 = sin(qJ(4));
t177 = t115 * t111;
t112 = cos(pkin(8));
t186 = cos(qJ(4));
t99 = t186 * t112;
t87 = -t99 + t177;
t194 = 0.2e1 * t87;
t88 = t111 * t186 + t112 * t115;
t84 = t88 ^ 2;
t85 = t87 ^ 2;
t193 = t84 + t85;
t149 = t85 - t84;
t114 = sin(qJ(5));
t116 = cos(qJ(5));
t157 = qJD(1) * t116;
t140 = t114 * t157;
t109 = t114 ^ 2;
t110 = t116 ^ 2;
t97 = t110 - t109;
t117 = t97 * qJD(4) + t140 * t194;
t192 = t84 / 0.2e1 + t85 / 0.2e1;
t46 = t193 * t116;
t191 = t87 * pkin(4);
t190 = t88 * pkin(7);
t189 = -t114 / 0.2e1;
t187 = t116 / 0.2e1;
t113 = -pkin(1) - qJ(3);
t185 = -pkin(6) + t113;
t131 = pkin(4) * t88 + pkin(7) * t87;
t98 = pkin(3) * t111 + qJ(2);
t121 = t131 + t98;
t136 = t185 * t112;
t91 = t185 * t111;
t57 = t115 * t136 + t186 * t91;
t183 = t114 * t57;
t11 = -t116 * t121 + t183;
t37 = t87 * t114;
t40 = t116 * t88;
t58 = t190 - t191;
t1 = t11 * t87 - t37 * t57 + t40 * t58;
t184 = t1 * qJD(1);
t35 = t114 * t88;
t182 = t116 * t57;
t12 = t114 * t121 + t182;
t42 = t87 * t116;
t2 = t12 * t87 - t35 * t58 - t42 * t57;
t181 = t2 * qJD(1);
t56 = t115 * t91 - t136 * t186;
t7 = t11 * t88 + t37 * t56;
t180 = t7 * qJD(1);
t8 = -t12 * t88 - t42 * t56;
t179 = t8 * qJD(1);
t178 = qJD(2) * t88;
t128 = 0.1e1 / 0.2e1 + t192;
t13 = t128 * t114;
t176 = t13 * qJD(1);
t17 = t128 * t116;
t175 = t17 * qJD(1);
t19 = t149 * t114;
t174 = t19 * qJD(1);
t20 = t193 * t114;
t173 = t20 * qJD(1);
t21 = t149 * t116;
t172 = t21 * qJD(1);
t171 = t149 * qJD(1);
t170 = t35 * qJD(1);
t169 = t37 * qJD(1);
t168 = t40 * qJD(1);
t167 = t46 * qJD(1);
t83 = t99 / 0.2e1 - t177 / 0.2e1;
t166 = t83 * qJD(1);
t105 = t111 ^ 2;
t106 = t112 ^ 2;
t94 = t105 + t106;
t86 = t94 * t113;
t165 = t86 * qJD(1);
t164 = t87 * qJD(1);
t76 = t87 * qJD(4);
t163 = t88 * qJD(1);
t79 = t88 * qJD(4);
t137 = -t105 / 0.2e1 - t106 / 0.2e1;
t93 = -0.1e1 / 0.2e1 + t137;
t162 = t93 * qJD(1);
t161 = t94 * qJD(1);
t156 = qJD(3) * t116;
t155 = qJD(4) * t114;
t154 = qJD(4) * t116;
t153 = qJD(5) * t114;
t152 = qJD(5) * t116;
t151 = t111 * qJD(1);
t150 = t112 * qJD(1);
t148 = t87 * t163;
t147 = t87 * t79;
t146 = t110 * t164;
t145 = t88 * t153;
t144 = t88 * t152;
t143 = t114 * t163;
t142 = t87 * t157;
t141 = t88 * t157;
t139 = t114 * t152;
t138 = t114 * t154;
t135 = qJD(1) * t98 + qJD(3);
t134 = qJD(5) + t163;
t132 = t87 * t138;
t130 = qJD(4) * t37 - t144;
t129 = t134 * t87;
t127 = t190 / 0.2e1 - t191 / 0.2e1;
t120 = t58 / 0.2e1 + t127;
t5 = t120 * t116;
t126 = pkin(4) * t155 + qJD(1) * t5;
t3 = t120 * t114;
t125 = pkin(4) * t154 - qJD(1) * t3;
t124 = -qJD(5) * t83 + t148;
t123 = t116 * t129;
t34 = (-t109 / 0.2e1 + t110 / 0.2e1) * t87;
t122 = -qJD(1) * t34 + t138;
t119 = qJD(4) * t34 + t140 * t85;
t45 = t97 * t85;
t118 = -qJD(1) * t45 + 0.2e1 * t132;
t108 = qJ(2) * qJD(2);
t107 = qJD(1) * qJ(2);
t92 = 0.1e1 / 0.2e1 + t137;
t75 = t83 * qJD(4);
t74 = t87 * t154;
t26 = t35 * qJD(5);
t25 = t34 * qJD(5);
t22 = -t153 - t170;
t18 = t187 - t46 / 0.2e1;
t14 = t114 * t192 + t189;
t6 = t56 * t114 - t116 * t127 + t187 * t58;
t4 = t114 * t127 + t56 * t116 + t189 * t58;
t9 = [0, 0, 0, 0, qJD(2), t108, qJD(2) * t111, qJD(2) * t112, t94 * qJD(3), -qJD(3) * t86 + t108, t147, -t149 * qJD(4), 0, 0, 0, -t76 * t98 + t178, -qJD(2) * t87 - t79 * t98, t110 * t147 - t139 * t85, -qJD(5) * t45 - 0.2e1 * t132 * t88, qJD(4) * t21 + t145 * t87, -qJD(4) * t19 + t144 * t87, -t147, qJD(3) * t20 + qJD(4) * t1 + qJD(5) * t8 + t116 * t178, qJD(3) * t46 + qJD(4) * t2 + qJD(5) * t7 - t114 * t178; 0, 0, 0, 0, qJD(1), t107, t151, t150, 0, qJD(3) * t92 + t107, 0, 0, 0, 0, 0, t163, -t164, 0, 0, 0, 0, 0, qJD(5) * t18 + t141, qJD(5) * t14 - t143; 0, 0, 0, 0, 0, 0, 0, 0, t161, qJD(2) * t92 - t165, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t171, -t79, t76, 0, -qJD(4) * t57 - t164 * t98, qJD(4) * t56 - t163 * t98, -t25 + (-t138 + t146) * t88, -t117 * t88 + t139 * t194, -t155 * t87 + t172, -t74 - t174, -t124, t184 + (t114 * t131 - t182) * qJD(4) + t6 * qJD(5), t181 + (t116 * t131 + t183) * qJD(4) + t4 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t118, t114 * t129, t123, t75, qJD(2) * t18 + qJD(4) * t6 - qJD(5) * t12 + t179, qJD(2) * t14 + qJD(4) * t4 + qJD(5) * t11 + t180; 0, 0, 0, 0, -qJD(1), -t107, -t151, -t150, 0, qJD(3) * t93 - t107, 0, 0, 0, 0, 0, -t163, t164, 0, 0, 0, 0, 0, -qJD(5) * t17 - t141, qJD(5) * t13 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t76, 0, 0, 0, 0, 0, qJD(5) * t37 - t154 * t88, qJD(5) * t42 + t155 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130 - t175, qJD(4) * t42 + t145 + t176; 0, 0, 0, 0, 0, 0, 0, 0, -t161, -qJD(2) * t93 + t165, 0, 0, 0, 0, 0, -t76, -t79, 0, 0, 0, 0, 0, -t26 - t74 - t173, t130 - t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, -t163, 0, 0, 0, 0, 0, -t142, t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t134 * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t171, 0, 0, 0, t135 * t87, t135 * t88, -t146 * t88 - t25, 0.2e1 * t114 * t123, qJD(5) * t40 - t172, -t26 + t174, t124, -qJD(5) * t5 + t156 * t87 - t184, -qJD(3) * t37 + qJD(5) * t3 - t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t163, 0, 0, 0, 0, 0, t142, -t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t97 * qJD(5), 0, 0, 0, -pkin(4) * t153, -pkin(4) * t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t117, t152 + t168, t22, -t166, -pkin(7) * t152 - t126, pkin(7) * t153 - t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t118, -qJD(4) * t40 - t143 * t87, qJD(4) * t35 - t141 * t87, t75, qJD(2) * t17 + qJD(3) * t35 + qJD(4) * t5 - t179, -qJD(2) * t13 - qJD(4) * t3 + t156 * t88 - t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, -t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t117, -t168, t170, t166, t126, t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
