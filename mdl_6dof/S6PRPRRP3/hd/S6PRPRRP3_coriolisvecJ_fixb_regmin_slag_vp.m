% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:45
% EndTime: 2019-03-08 20:07:51
% DurationCPUTime: 1.85s
% Computational Cost: add. (2870->272), mult. (7448->386), div. (0->0), fcn. (5886->10), ass. (0->148)
t115 = sin(pkin(11));
t184 = pkin(8) + qJ(3);
t100 = t184 * t115;
t117 = cos(pkin(11));
t101 = t184 * t117;
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t66 = t100 * t123 + t120 * t101;
t170 = t115 * t120;
t96 = -t123 * t117 + t170;
t40 = -t96 * qJD(3) - qJD(4) * t66;
t116 = sin(pkin(6));
t124 = cos(qJ(2));
t167 = t116 * t124;
t127 = t96 * t167;
t69 = qJD(1) * t127;
t202 = t40 + t69;
t121 = sin(qJ(2));
t163 = qJD(1) * t116;
t151 = t121 * t163;
t91 = t96 * qJD(4);
t97 = t115 * t123 + t117 * t120;
t92 = t97 * qJD(4);
t201 = pkin(4) * t92 + pkin(9) * t91 - t151;
t118 = cos(pkin(6));
t162 = qJD(1) * t118;
t107 = t117 * t162;
t99 = qJD(2) * qJ(3) + t151;
t63 = t107 + (-pkin(8) * qJD(2) - t99) * t115;
t161 = qJD(2) * t117;
t76 = t115 * t162 + t117 * t99;
t64 = pkin(8) * t161 + t76;
t27 = t120 * t63 + t123 * t64;
t200 = qJD(4) * t27;
t150 = t124 * t163;
t95 = (qJD(3) + t150) * qJD(2);
t199 = t96 * t95;
t119 = sin(qJ(5));
t122 = cos(qJ(5));
t198 = -t119 * t69 + t201 * t122;
t108 = t123 * t161;
t149 = qJD(2) * t170;
t88 = -t108 + t149;
t83 = qJD(5) + t88;
t143 = t119 * t83;
t90 = qJD(2) * t97;
t73 = qJD(4) * t119 + t122 * t90;
t197 = t73 * t143;
t158 = qJD(5) * t122;
t110 = -pkin(3) * t117 - pkin(2);
t56 = pkin(4) * t96 - pkin(9) * t97 + t110;
t196 = t201 * t119 + t202 * t122 + t56 * t158;
t195 = -t120 * t64 + t123 * t63;
t138 = qJD(3) - t150;
t23 = qJD(4) * pkin(9) + t27;
t82 = t110 * qJD(2) + t138;
t30 = pkin(4) * t88 - pkin(9) * t90 + t82;
t13 = t119 * t30 + t122 * t23;
t16 = qJD(4) * t195 - t199;
t160 = qJD(2) * t121;
t148 = t116 * t160;
t104 = qJD(1) * t148;
t105 = qJD(4) * t108;
t133 = qJD(4) * t149 - t105;
t81 = qJD(2) * t92;
t36 = t81 * pkin(4) + t133 * pkin(9) + t104;
t32 = t122 * t36;
t126 = -t13 * qJD(5) - t119 * t16 + t32;
t157 = t122 * qJD(4);
t159 = qJD(5) * t119;
t37 = -qJD(5) * t157 + t122 * t133 + t90 * t159;
t1 = pkin(5) * t81 + qJ(6) * t37 - qJD(6) * t73 + t126;
t71 = t119 * t90 - t157;
t8 = -qJ(6) * t71 + t13;
t194 = t8 * t83 + t1;
t193 = t73 ^ 2;
t12 = -t119 * t23 + t122 * t30;
t7 = -qJ(6) * t73 + t12;
t5 = pkin(5) * t83 + t7;
t192 = t5 - t7;
t135 = qJ(6) * t91 - qJD(6) * t97;
t67 = -t100 * t120 + t101 * t123;
t62 = t122 * t67;
t191 = pkin(5) * t92 - t119 * t40 + t135 * t122 + (-t62 + (qJ(6) * t97 - t56) * t119) * qJD(5) + t198;
t154 = t97 * t158;
t190 = -qJ(6) * t154 + (-qJD(5) * t67 + t135) * t119 + t196;
t183 = -qJ(6) - pkin(9);
t147 = qJD(5) * t183;
t172 = qJ(6) * t122;
t54 = pkin(4) * t90 + pkin(9) * t88;
t49 = t122 * t54;
t189 = -pkin(5) * t90 + t122 * t147 - t88 * t172 - t49 + (-qJD(6) + t195) * t119;
t188 = t71 * t88;
t187 = t71 * t90;
t186 = t73 * t90;
t185 = t97 * t81;
t182 = t119 * t54 + t122 * t195;
t129 = t119 * t133;
t38 = t73 * qJD(5) - t129;
t181 = -t119 * t38 - t71 * t158;
t128 = t97 * t167;
t180 = -qJD(1) * t128 + t97 * qJD(3) + t67 * qJD(4);
t179 = t119 * t56 + t62;
t173 = qJ(6) * t119;
t178 = qJD(6) * t122 + t119 * t147 - t88 * t173 - t182;
t177 = qJD(2) * pkin(2);
t176 = t119 * t37;
t175 = t119 * t81;
t168 = t116 * t121;
t125 = qJD(2) ^ 2;
t166 = t116 * t125;
t164 = t115 ^ 2 + t117 ^ 2;
t155 = t97 * t159;
t153 = t121 * t166;
t152 = t124 * t166;
t146 = t164 * t95;
t17 = t97 * t95 + t200;
t144 = t122 * t83;
t142 = qJD(5) * pkin(9) * t83 + t17;
t130 = t119 * t36 + t122 * t16 + t30 * t158 - t23 * t159;
t2 = -qJ(6) * t38 - qJD(6) * t71 + t130;
t141 = -t83 * t5 + t2;
t22 = -qJD(4) * pkin(4) - t195;
t140 = t17 * t97 - t22 * t91;
t139 = -t83 * t91 + t185;
t137 = t115 * (-t115 * t99 + t107) - t117 * t76;
t86 = -t115 * t168 + t117 * t118;
t87 = t115 * t118 + t117 * t168;
t44 = t120 * t87 - t123 * t86;
t45 = t120 * t86 + t123 * t87;
t134 = t122 * t81 - t88 * t143 - t83 * t159;
t6 = pkin(5) * t38 + t17;
t33 = -t119 * t45 - t122 * t167;
t132 = t119 * t167 - t122 * t45;
t131 = -pkin(9) * t81 + t83 * t22;
t103 = t183 * t122;
t102 = t183 * t119;
t98 = t138 - t177;
t70 = t71 ^ 2;
t52 = t122 * t56;
t25 = qJD(2) * t128 + t45 * qJD(4);
t24 = -qJD(2) * t127 - t44 * qJD(4);
t20 = -t97 * t173 + t179;
t19 = pkin(5) * t71 + qJD(6) + t22;
t18 = pkin(5) * t96 - t119 * t67 - t97 * t172 + t52;
t11 = t132 * qJD(5) - t119 * t24 + t122 * t148;
t10 = t33 * qJD(5) + t119 * t148 + t122 * t24;
t3 = [0, 0, -t153, -t152, -t117 * t153, t115 * t153, t164 * t152 (-t115 * t86 + t117 * t87) * t95 + (t121 * t98 + (-t137 - t151) * t124) * t116 * qJD(2), 0, 0, 0, 0, 0, -qJD(4) * t25 + (-t124 * t81 + t88 * t160) * t116, -t24 * qJD(4) + (t124 * t133 + t90 * t160) * t116, 0, 0, 0, 0, 0, t11 * t83 + t25 * t71 + t33 * t81 + t38 * t44, -t10 * t83 + t132 * t81 + t25 * t73 - t37 * t44, -t10 * t71 - t11 * t73 + t132 * t38 + t33 * t37, t1 * t33 + t10 * t8 + t11 * t5 - t132 * t2 + t19 * t25 + t44 * t6; 0, 0, 0, 0, 0, 0, qJD(2) * t138 * t164 + t146, -t137 * qJD(3) + qJ(3) * t146 + (t137 * t124 + (-t98 - t177) * t121) * t163, -t133 * t97 - t90 * t91, t133 * t96 + t91 * t88 - t90 * t92 - t185, -t91 * qJD(4), -t92 * qJD(4), 0, t110 * t81 + t82 * t92 - t180 * qJD(4) + (qJD(2) * t96 - t88) * t151, t110 * t105 - t82 * t91 + (-t110 * t149 - t202) * qJD(4), -t73 * t155 + (-t37 * t97 - t73 * t91) * t122 -(-t119 * t73 - t122 * t71) * t91 + (t176 - t122 * t38 + (t119 * t71 - t122 * t73) * qJD(5)) * t97, t122 * t139 - t155 * t83 - t37 * t96 + t73 * t92, -t119 * t139 - t154 * t83 - t38 * t96 - t71 * t92, t81 * t96 + t83 * t92, t52 * t81 + (-t158 * t23 + t32) * t96 + t12 * t92 + t66 * t38 + t22 * t154 + (-t67 * t158 + t198) * t83 + t180 * t71 + ((-qJD(5) * t56 - t40) * t83 - t67 * t81 + (-qJD(5) * t30 - t16) * t96 + t140) * t119, -t179 * t81 - t130 * t96 - t13 * t92 - t66 * t37 - t22 * t155 + (t159 * t67 - t196) * t83 + t180 * t73 + t140 * t122, t18 * t37 - t20 * t38 - (-t119 * t8 - t122 * t5) * t91 - t191 * t73 - t190 * t71 + (-t1 * t122 - t119 * t2 + (t119 * t5 - t122 * t8) * qJD(5)) * t97, t2 * t20 + t1 * t18 + t6 * (pkin(5) * t119 * t97 + t66) + t190 * t8 + t191 * t5 + ((-t119 * t91 + t154) * pkin(5) + t180) * t19; 0, 0, 0, 0, 0, 0, -t164 * t125, t137 * qJD(2) + t104, 0, 0, 0, 0, 0, 0.2e1 * qJD(4) * t90, t105 + (-t88 - t149) * qJD(4), 0, 0, 0, 0, 0, t134 - t187, -t122 * t83 ^ 2 - t175 - t186 (t37 - t188) * t122 + t197 + t181, t141 * t119 + t122 * t194 - t19 * t90; 0, 0, 0, 0, 0, 0, 0, 0, t90 * t88, -t88 ^ 2 + t90 ^ 2, t105 + (t88 - t149) * qJD(4), 0, 0, -t82 * t90 - t17 + t200, t82 * t88 + t199, t144 * t73 - t176 (-t37 - t188) * t122 - t197 + t181, t144 * t83 + t175 - t186, t134 + t187, -t83 * t90, -pkin(4) * t38 - t12 * t90 - t27 * t71 - t49 * t83 - t142 * t122 + (t195 * t83 + t131) * t119, pkin(4) * t37 + t119 * t142 + t122 * t131 + t13 * t90 + t182 * t83 - t27 * t73, t102 * t37 + t103 * t38 - t119 * t194 + t141 * t122 - t178 * t71 - t189 * t73, -t2 * t103 + t1 * t102 + t6 * (-pkin(5) * t122 - pkin(4)) + t178 * t8 + t189 * t5 + (pkin(5) * t143 - t27) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t70 + t193, t71 * t83 - t37, t129 + (-qJD(5) + t83) * t73, t81, t13 * t83 - t22 * t73 + t126, t12 * t83 + t22 * t71 - t130, pkin(5) * t37 - t192 * t71, t192 * t8 + (-t19 * t73 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70 - t193, t5 * t73 + t71 * t8 + t6;];
tauc_reg  = t3;
