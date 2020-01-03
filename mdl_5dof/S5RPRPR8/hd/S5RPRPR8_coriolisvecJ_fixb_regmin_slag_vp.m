% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:11
% EndTime: 2019-12-31 18:22:16
% DurationCPUTime: 1.33s
% Computational Cost: add. (1245->214), mult. (3077->331), div. (0->0), fcn. (2067->8), ass. (0->126)
t105 = cos(qJ(3));
t141 = t105 * qJD(1);
t88 = -qJD(5) + t141;
t168 = qJD(5) + t88;
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t100 = cos(pkin(9));
t142 = t100 * qJD(3);
t103 = sin(qJ(3));
t148 = qJD(1) * t103;
t98 = sin(pkin(9));
t72 = t98 * t148 - t142;
t149 = t98 * qJD(3);
t74 = t100 * t148 + t149;
t120 = t102 * t72 - t104 * t74;
t167 = t120 * t88;
t90 = sin(pkin(8)) * pkin(1) + pkin(6);
t86 = t90 * qJD(1);
t166 = t105 * qJD(2) - t103 * t86;
t165 = t120 * qJD(5);
t164 = pkin(7) + qJ(4);
t139 = qJD(1) * qJD(3);
t131 = t103 * t139;
t77 = t102 * t100 + t104 * t98;
t111 = t77 * t105;
t108 = qJD(3) * t111;
t143 = qJD(5) * t104;
t144 = qJD(5) * t103;
t154 = t100 * t103;
t158 = t102 * t98;
t21 = t143 * t154 - t144 * t158 + t108;
t55 = t77 * t103;
t163 = -t55 * t131 + t21 * t88;
t42 = (qJD(4) + t166) * qJD(3);
t123 = pkin(3) * t103 - qJ(4) * t105;
t63 = t123 * qJD(3) - t103 * qJD(4);
t50 = t63 * qJD(1);
t11 = t100 * t42 + t98 * t50;
t147 = qJD(3) * qJ(4);
t95 = t103 * qJD(2);
t60 = t105 * t86 + t95;
t46 = t147 + t60;
t91 = -cos(pkin(8)) * pkin(1) - pkin(2);
t69 = -t105 * pkin(3) - t103 * qJ(4) + t91;
t49 = t69 * qJD(1);
t14 = t100 * t46 + t98 * t49;
t81 = t123 * qJD(1);
t25 = t100 * t166 + t98 * t81;
t152 = t104 * t100;
t119 = t152 - t158;
t112 = t119 * t105;
t162 = -qJD(1) * t112 + t119 * qJD(5);
t161 = -qJD(1) * t111 + t77 * qJD(5);
t145 = qJD(3) * t105;
t51 = qJD(3) * t95 + t86 * t145;
t146 = qJD(3) * t103;
t136 = t90 * t146;
t27 = t100 * t63 + t98 * t136;
t153 = t100 * t105;
t80 = t90 * t153;
t33 = t98 * t69 + t80;
t132 = t105 * t139;
t160 = t132 * t152 - t72 * t143;
t96 = t103 ^ 2;
t159 = -t105 ^ 2 + t96;
t157 = t103 * t98;
t156 = t105 * t98;
t87 = qJD(1) * t91;
t155 = qJD(1) * t96;
t106 = qJD(3) ^ 2;
t151 = t106 * t103;
t150 = t106 * t105;
t138 = t90 * t156;
t137 = t98 * t141;
t135 = pkin(4) * t98 + t90;
t10 = t100 * t50 - t98 * t42;
t116 = pkin(4) * t103 - pkin(7) * t153;
t109 = t116 * qJD(3);
t5 = qJD(1) * t109 + t10;
t129 = t98 * t132;
t7 = -pkin(7) * t129 + t11;
t134 = -t102 * t7 + t104 * t5;
t8 = (-qJD(5) * t74 - t129) * t102 + t160;
t133 = -t8 * t105 - t120 * t146;
t13 = t100 * t49 - t98 * t46;
t24 = t100 * t81 - t166 * t98;
t130 = -qJD(3) * pkin(3) + qJD(4);
t44 = t130 - t166;
t128 = t130 - t44;
t4 = -pkin(4) * t141 - t74 * pkin(7) + t13;
t6 = -t72 * pkin(7) + t14;
t127 = t102 * t6 - t104 * t4;
t2 = t102 * t4 + t104 * t6;
t126 = t102 * t5 + t104 * t7;
t125 = -t10 * t98 + t11 * t100;
t124 = t100 * t14 - t13 * t98;
t58 = t100 * t69;
t19 = -pkin(7) * t154 + t58 + (-t90 * t98 - pkin(4)) * t105;
t23 = -pkin(7) * t157 + t33;
t122 = -t102 * t23 + t104 * t19;
t121 = t102 * t19 + t104 * t23;
t118 = 0.2e1 * qJD(3) * t87;
t64 = t104 * t72;
t29 = t102 * t74 + t64;
t9 = qJD(1) * t108 - t165;
t115 = t105 * t9 - t29 * t146;
t85 = t164 * t100;
t114 = t116 * qJD(1) + qJD(4) * t98 + qJD(5) * t85 + t24;
t84 = t164 * t98;
t113 = pkin(7) * t137 + qJD(4) * t100 - qJD(5) * t84 - t25;
t20 = qJD(3) * t112 - t77 * t144;
t56 = t119 * t103;
t110 = -t56 * t131 + t20 * t88;
t107 = qJD(1) ^ 2;
t92 = -t100 * pkin(4) - pkin(3);
t62 = t135 * t103;
t54 = t135 * t145;
t52 = t98 * t63;
t36 = pkin(4) * t137 + t60;
t34 = pkin(4) * t129 + t51;
t32 = t58 - t138;
t28 = -t100 * t136 + t52;
t26 = t72 * pkin(4) + t44;
t18 = t52 + (-pkin(7) * t156 - t90 * t154) * qJD(3);
t16 = t109 + t27;
t1 = [0, 0, 0, 0, 0.2e1 * t105 * t131, -0.2e1 * t159 * t139, t150, -t151, 0, t103 * t118 - t90 * t150, t105 * t118 + t90 * t151, t51 * t157 + (-qJD(1) * t27 - t10) * t105 + ((t44 * t98 + t72 * t90) * t105 + (t13 + (t32 + t138) * qJD(1)) * t103) * qJD(3), t51 * t154 + (qJD(1) * t28 + t11) * t105 + ((t100 * t44 + t74 * t90) * t105 + (-t14 + (-t33 + t80) * qJD(1)) * t103) * qJD(3), -t27 * t74 - t28 * t72 + (-t10 * t100 - t11 * t98) * t103 + (-t100 * t13 - t14 * t98 + (-t100 * t32 - t33 * t98) * qJD(1)) * t145, t10 * t32 + t11 * t33 + t13 * t27 + t14 * t28 + (t103 * t51 + t145 * t44) * t90, -t120 * t20 + t8 * t56, t120 * t21 - t20 * t29 - t8 * t55 - t56 * t9, -t110 + t133, t115 + t163, (-t88 - t141) * t146, -(-t102 * t18 + t104 * t16) * t88 - t134 * t105 + t54 * t29 + t62 * t9 + t34 * t55 + t26 * t21 + (t105 * t2 + t121 * t88) * qJD(5) + (qJD(1) * t122 - t127) * t146, (t102 * t16 + t104 * t18) * t88 + t126 * t105 - t54 * t120 + t62 * t8 + t34 * t56 + t26 * t20 + (-t105 * t127 + t122 * t88) * qJD(5) + (-qJD(1) * t121 - t2) * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, -t150, (t103 * t72 - t98 * t155) * qJD(3), (-t100 * t155 + t103 * t74) * qJD(3), (-t100 * t72 + t74 * t98) * t145, -t51 * t105 + t125 * t103 + (t103 * t44 + t105 * t124) * qJD(3), 0, 0, 0, 0, 0, -t115 + t163, t110 + t133; 0, 0, 0, 0, -t103 * t107 * t105, t159 * t107, 0, 0, 0, t60 * qJD(3) - t87 * t148 - t51, -t87 * t141, -t51 * t100 - t60 * t72 + ((-t147 * t98 - t13) * t103 + (t128 * t98 + t24) * t105) * qJD(1), t51 * t98 - t60 * t74 + (t103 * t14 - t105 * t25 + (-qJ(4) * t146 + t105 * t128) * t100) * qJD(1), t24 * t74 + t25 * t72 + (qJD(4) * t74 + t14 * t141 - t10) * t98 + (-qJD(4) * t72 + t13 * t141 + t11) * t100, -t51 * pkin(3) + t125 * qJ(4) + qJD(4) * t124 - t13 * t24 - t14 * t25 - t44 * t60, -t120 * t162 + t8 * t77, t119 * t8 + t120 * t161 - t162 * t29 - t77 * t9, -t162 * t88 + (qJD(3) * t77 + t120) * t148, t161 * t88 + (qJD(3) * t119 + t29) * t148, t88 * t148, -t36 * t29 - t34 * t119 + t92 * t9 + (t102 * t113 + t104 * t114) * t88 + t161 * t26 + ((-t102 * t85 - t104 * t84) * qJD(3) + t127) * t148, t36 * t120 + t34 * t77 + t92 * t8 + (-t102 * t114 + t104 * t113) * t88 + t162 * t26 + (-(-t102 * t84 + t104 * t85) * qJD(3) + t2) * t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t74 + t149) * t141, (t72 + t142) * t141, -t72 ^ 2 - t74 ^ 2, t13 * t74 + t14 * t72 + t51, 0, 0, 0, 0, 0, t9 + t167, t64 * t88 + (-t129 + (-qJD(5) + t88) * t74) * t102 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120 * t29, t120 ^ 2 - t29 ^ 2, -t29 * t88 + t8, -t132 * t77 + t165 + t167, t131, t26 * t120 - t168 * t2 + t134, t168 * t127 + t26 * t29 - t126;];
tauc_reg = t1;
