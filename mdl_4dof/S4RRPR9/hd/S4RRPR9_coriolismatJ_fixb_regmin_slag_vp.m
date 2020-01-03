% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:06
% EndTime: 2019-12-31 17:10:09
% DurationCPUTime: 1.21s
% Computational Cost: add. (1081->181), mult. (2650->295), div. (0->0), fcn. (2671->6), ass. (0->152)
t108 = sin(pkin(7));
t112 = cos(qJ(4));
t163 = t112 * t108;
t109 = cos(pkin(7));
t110 = sin(qJ(4));
t164 = t110 * t109;
t81 = t163 + t164;
t190 = -t81 / 0.2e1;
t162 = t112 * t109;
t165 = t110 * t108;
t192 = t162 - t165;
t191 = t192 / 0.2e1;
t102 = -t109 * pkin(3) - pkin(2);
t189 = -t102 / 0.2e1;
t113 = cos(qJ(2));
t188 = -t113 / 0.2e1;
t187 = t113 / 0.2e1;
t186 = t113 * pkin(5);
t185 = pkin(6) + qJ(3);
t111 = sin(qJ(2));
t93 = t111 * pkin(2) - t113 * qJ(3);
t169 = t108 * t111;
t97 = pkin(5) * t169;
t65 = t109 * t93 + t97;
t130 = -t113 * pkin(2) - t111 * qJ(3);
t88 = -pkin(1) + t130;
t166 = t109 * t113;
t98 = pkin(5) * t166;
t62 = t108 * t88 + t98;
t167 = t109 * t111;
t76 = t109 * t88;
t40 = -pkin(6) * t167 + t76 + (-pkin(5) * t108 - pkin(3)) * t113;
t49 = -pkin(6) * t169 + t62;
t16 = t110 * t49 - t112 * t40;
t67 = t81 * t111;
t136 = pkin(3) * t108 + pkin(5);
t84 = t136 * t111;
t8 = -t16 * t113 - t84 * t67;
t184 = qJD(1) * t8;
t17 = t110 * t40 + t112 * t49;
t70 = t192 * t111;
t9 = -t17 * t113 - t84 * t70;
t183 = qJD(1) * t9;
t43 = t111 * pkin(3) - pkin(6) * t166 + t65;
t179 = t112 * t43;
t168 = t108 * t113;
t66 = -pkin(5) * t167 + t108 * t93;
t50 = -pkin(6) * t168 + t66;
t180 = t110 * t50;
t69 = t81 * t113;
t85 = t136 * t113;
t1 = (t179 - t180) * t113 + t16 * t111 - t85 * t67 - t84 * t69;
t182 = t1 * qJD(1);
t181 = t110 * t43;
t178 = t112 * t50;
t71 = t192 * t113;
t2 = (t178 + t181) * t113 - t17 * t111 + t85 * t70 + t84 * t71;
t177 = t2 * qJD(1);
t91 = t185 * t108;
t92 = t185 * t109;
t52 = -t110 * t91 + t112 * t92;
t115 = t52 * t188 + t70 * t189 + t84 * t190;
t122 = -t180 / 0.2e1 + t179 / 0.2e1;
t3 = t115 + t122;
t176 = t3 * qJD(1);
t51 = t110 * t92 + t112 * t91;
t116 = -t67 * t189 + t51 * t187 - t84 * t192 / 0.2e1;
t123 = -t181 / 0.2e1 - t178 / 0.2e1;
t4 = t116 + t123;
t175 = t4 * qJD(1);
t61 = -pkin(5) * t168 + t76;
t7 = (t111 * t65 + t113 * t61) * t109 + (t111 * t66 + t113 * t62) * t108;
t174 = t7 * qJD(1);
t26 = (t108 * t62 + t109 * t61) * t111;
t173 = qJD(1) * t26;
t172 = qJD(1) * t70;
t171 = qJD(2) * t81;
t170 = qJD(4) * t81;
t12 = t111 * pkin(5) ^ 2 * t113 + t61 * t65 + t62 * t66;
t161 = t12 * qJD(1);
t14 = -t67 * t71 - t69 * t70;
t160 = t14 * qJD(1);
t21 = -t61 * t111 + (t65 - 0.2e1 * t97) * t113;
t159 = t21 * qJD(1);
t22 = t66 * t113 + (-t62 + 0.2e1 * t98) * t111;
t158 = t22 * qJD(1);
t27 = t67 * t111 - t69 * t113;
t157 = t27 * qJD(1);
t28 = -t70 * t111 + t71 * t113;
t156 = t28 * qJD(1);
t117 = -t164 / 0.2e1 - t163 / 0.2e1;
t30 = (t190 - t117) * t113;
t155 = t30 * qJD(1);
t31 = (t190 + t117) * t113;
t154 = t31 * qJD(1);
t32 = t162 * t187 + (t165 - t192) * t188;
t153 = t32 * qJD(1);
t33 = (t191 - t162 / 0.2e1 + t165 / 0.2e1) * t113;
t152 = t33 * qJD(1);
t107 = t111 ^ 2;
t95 = t108 ^ 2 + t109 ^ 2;
t78 = t95 * t107;
t151 = t78 * qJD(1);
t74 = t192 * qJD(4);
t150 = t95 * qJD(2);
t96 = t113 ^ 2 - t107;
t149 = t96 * qJD(1);
t148 = qJD(2) * t102;
t147 = qJD(3) * t113;
t146 = qJD(4) * t102;
t145 = qJD(4) * t113;
t144 = t111 * qJD(1);
t143 = t111 * qJD(2);
t142 = t113 * qJD(1);
t141 = t113 * qJD(2);
t140 = pkin(1) * t144;
t139 = pkin(1) * t142;
t138 = pkin(5) * t141;
t137 = t186 / 0.2e1;
t135 = t111 * t147;
t134 = t111 * t141;
t133 = t111 * t142;
t60 = t70 * t142;
t132 = qJD(2) * t31 - t60;
t131 = qJD(3) + t148;
t129 = -t65 * t108 + t66 * t109;
t11 = -t192 * t67 - t81 * t70;
t25 = t67 ^ 2 - t70 ^ 2;
t128 = qJD(1) * t25 + qJD(2) * t11;
t29 = t192 ^ 2 - t81 ^ 2;
t127 = qJD(1) * t11 + qJD(2) * t29;
t121 = t61 * t108 / 0.2e1 - t62 * t109 / 0.2e1;
t23 = t137 + t121;
t86 = t95 * qJ(3);
t126 = qJD(1) * t23 - qJD(2) * t86;
t125 = -qJD(1) * t67 + qJD(2) * t192;
t124 = t171 + t172;
t19 = t67 * t190 + t70 * t191;
t120 = qJD(2) * t19 - t67 * t172;
t119 = -qJD(1) * t19 - t171 * t192;
t118 = qJD(2) * t30 + qJD(4) * t70 - t60;
t114 = t130 * qJD(2) + t147;
t103 = t143 / 0.2e1;
t77 = (t142 - qJD(4) / 0.2e1) * t111;
t24 = t137 - t121;
t20 = t32 * qJD(2) - t67 * t142;
t18 = t19 * qJD(4);
t13 = -t33 * qJD(2) + (-qJD(4) + t142) * t67;
t10 = t11 * qJD(4);
t6 = -t115 + t122;
t5 = -t116 + t123;
t15 = [0, 0, 0, t134, t96 * qJD(2), 0, 0, 0, -pkin(1) * t143, -pkin(1) * t141, -t21 * qJD(2) + t109 * t135, t22 * qJD(2) - t108 * t135, -qJD(2) * t7 + qJD(3) * t78, qJD(2) * t12 - qJD(3) * t26, (qJD(2) * t71 - t67 * qJD(4)) * t70, qJD(2) * t14 + qJD(4) * t25, -t28 * qJD(2) + t67 * t145, -t27 * qJD(2) + t70 * t145, -t134, -t1 * qJD(2) - t9 * qJD(4) + t70 * t147, t2 * qJD(2) + t8 * qJD(4) - t147 * t67; 0, 0, 0, t133, t149, t141, -t143, 0, -t138 - t140, pkin(5) * t143 - t139, -t98 * qJD(2) + t114 * t108 - t159, t108 * t138 + t114 * t109 + t158, t129 * qJD(2) - t174, t161 + (-pkin(2) * t186 + t129 * qJ(3)) * qJD(2) + t24 * qJD(3), t124 * t71 + t18, t160 + (t192 * t71 - t69 * t81) * qJD(2) + t10, -t33 * qJD(4) + t81 * t143 - t156, -t30 * qJD(4) + t143 * t192 - t157, -t77, -t182 + (t102 * t69 - t51 * t111 - t192 * t85) * qJD(2) - t31 * qJD(3) + t6 * qJD(4), t177 + (t102 * t71 - t52 * t111 + t85 * t81) * qJD(2) + t32 * qJD(3) + t5 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (qJD(2) * t108 + t109 * t144) * t113, (qJD(2) * t109 - t108 * t144) * t113, t151, qJD(2) * t24 - t173, 0, 0, 0, 0, 0, -t132, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t128, t13, -t118, t103, qJD(2) * t6 - qJD(4) * t17 - t183, qJD(2) * t5 + qJD(4) * t16 + t184; 0, 0, 0, -t133, -t149, 0, 0, 0, t140, t139, t159, -t158, t174, -qJD(3) * t23 - t161, -t71 * t172 + t18, t10 - t160, -qJD(4) * t32 + t156, -qJD(4) * t31 + t157, t77, -qJD(3) * t30 - qJD(4) * t3 + t182, qJD(3) * t33 - qJD(4) * t4 - t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * qJD(3), t86 * qJD(3), t81 * t74, t29 * qJD(4), 0, 0, 0, t81 * t146, t192 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, -t126, 0, 0, 0, 0, 0, -t155, t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t127, t74 - t153, -t154 - t170, -t144 / 0.2e1, -t52 * qJD(4) + t148 * t81 - t176, t51 * qJD(4) + t148 * t192 - t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109 * t133, t108 * t133, -t151, qJD(2) * t23 + t173, 0, 0, 0, 0, 0, t118, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, t126, 0, 0, 0, 0, 0, t155 + t170, t74 - t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t128, t20, t132, t103, qJD(2) * t3 - qJD(3) * t70 + t183, qJD(2) * t4 + qJD(3) * t67 - t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t127, t153, t154, t144 / 0.2e1, -t131 * t81 + t176, -t131 * t192 + t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t15;
