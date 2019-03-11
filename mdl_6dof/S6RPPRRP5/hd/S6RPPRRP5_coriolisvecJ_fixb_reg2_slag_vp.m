% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:47
% EndTime: 2019-03-09 02:08:52
% DurationCPUTime: 2.25s
% Computational Cost: add. (2576->313), mult. (5115->421), div. (0->0), fcn. (2799->4), ass. (0->174)
t103 = sin(qJ(5));
t105 = cos(qJ(5));
t158 = t105 * qJD(4);
t106 = cos(qJ(4));
t167 = qJD(1) * t106;
t74 = t103 * t167 - t158;
t104 = sin(qJ(4));
t168 = qJD(1) * t104;
t90 = qJD(5) + t168;
t120 = t74 * t90;
t159 = qJD(5) * t106;
t144 = t103 * t159;
t42 = (t104 * t158 + t144) * qJD(1) - qJD(5) * t158;
t210 = -t42 - t120;
t29 = -t42 + t120;
t147 = t105 * t167;
t165 = qJD(4) * t103;
t76 = t147 + t165;
t196 = t76 * t90;
t155 = qJD(1) * qJD(4);
t142 = t104 * t155;
t175 = qJD(5) * t76;
t43 = -t103 * t142 + t175;
t209 = t43 - t196;
t208 = t43 + t196;
t160 = qJD(5) * t105;
t162 = qJD(5) * t103;
t156 = qJD(1) * qJD(2);
t163 = qJD(4) * t106;
t94 = qJD(1) * qJ(2) + qJD(3);
t87 = -pkin(7) * qJD(1) + t94;
t55 = t104 * t156 + t87 * t163;
t102 = pkin(1) + qJ(3);
t80 = pkin(4) * t104 - pkin(8) * t106 + t102;
t56 = t80 * qJD(1) - qJD(2);
t128 = pkin(4) * t106 + pkin(8) * t104;
t72 = t128 * qJD(4) + qJD(3);
t57 = t72 * qJD(1);
t79 = t104 * t87;
t63 = qJD(4) * pkin(8) + t79;
t139 = -t103 * t57 - t105 * t55 - t56 * t160 + t63 * t162;
t32 = -t103 * t63 + t105 * t56;
t207 = -t32 * t90 - t139;
t101 = -pkin(7) + qJ(2);
t206 = qJD(2) * t104 + t101 * t163;
t33 = t103 * t56 + t105 * t63;
t11 = -qJD(5) * t33 - t103 * t55 + t105 * t57;
t205 = -t33 * t90 - t11;
t204 = qJD(1) * t102;
t100 = t106 ^ 2;
t99 = t104 ^ 2;
t176 = t100 + t99;
t203 = t176 * qJD(2);
t202 = t76 ^ 2;
t201 = pkin(5) * t74;
t23 = -qJ(6) * t74 + t33;
t200 = t23 * t90;
t197 = t76 * t74;
t195 = -qJ(6) - pkin(8);
t22 = -qJ(6) * t76 + t32;
t19 = pkin(5) * t90 + t22;
t194 = t19 - t22;
t140 = qJD(5) * t195;
t172 = t104 * t105;
t173 = t103 * t106;
t78 = t128 * qJD(1);
t39 = t105 * t78 - t87 * t173;
t193 = (pkin(5) * t106 + qJ(6) * t172) * qJD(1) + t39 + t103 * qJD(6) - t105 * t140;
t148 = t103 * t168;
t171 = t105 * t106;
t40 = t103 * t78 + t87 * t171;
t192 = qJ(6) * t148 - t105 * qJD(6) - t103 * t140 + t40;
t46 = t101 * t172 + t103 * t80;
t164 = qJD(4) * t104;
t54 = -t106 * t156 + t87 * t164;
t28 = pkin(5) * t43 + t54;
t191 = t103 * t28;
t190 = t103 * t54;
t64 = -qJD(4) * pkin(4) - t106 * t87;
t189 = t103 * t64;
t188 = t103 * t76;
t187 = t103 * t90;
t186 = t105 * t28;
t185 = t105 * t54;
t184 = t105 * t64;
t183 = t105 * t90;
t182 = t106 * t42;
t181 = t106 * t43;
t180 = t106 * t76;
t108 = qJD(1) ^ 2;
t179 = t108 * t99;
t178 = t42 * t103;
t177 = t43 * t105;
t174 = t103 * t104;
t88 = -qJD(2) + t204;
t170 = qJD(2) - t88;
t107 = qJD(4) ^ 2;
t169 = -t107 - t108;
t161 = qJD(5) * t104;
t157 = t106 * qJD(2);
t154 = t90 * t174;
t153 = t90 * t172;
t97 = 0.2e1 * t156;
t152 = 0.2e1 * qJD(3) * qJD(1);
t151 = t90 * t162;
t150 = t106 * t108 * t104;
t149 = t90 * t167;
t145 = t101 * t161;
t143 = t105 * t159;
t91 = t106 * t155;
t141 = -qJD(6) - t201;
t138 = t103 * t72 + t206 * t105 + t80 * t160;
t137 = qJD(5) * t74 - t42;
t136 = -t43 + t175;
t135 = t88 + t204;
t134 = t90 + t168;
t133 = t170 * qJD(1);
t132 = qJD(1) + t161;
t131 = pkin(5) * t91;
t130 = t104 * t91;
t129 = t103 * t91 + t90 * t160;
t127 = qJD(2) + t135;
t126 = t103 * t23 + t105 * t19;
t125 = t103 * t19 - t105 * t23;
t124 = t103 * t33 + t105 * t32;
t123 = t103 * t32 - t105 * t33;
t122 = qJD(1) * t100 - t104 * t90;
t119 = -t101 * t107 + t152;
t118 = -pkin(8) * t163 + t104 * t64;
t117 = qJ(6) * t43 + t139;
t114 = -qJD(6) * t106 + (qJ(6) * qJD(4) - qJD(5) * t101) * t104;
t109 = qJ(6) * t42 + t11;
t3 = -qJD(6) * t76 + t109 + t131;
t6 = -qJD(6) * t74 - t117;
t112 = t125 * qJD(5) - t103 * t6 - t105 * t3;
t111 = t123 * qJD(5) + t103 * t139 - t105 * t11;
t110 = -t124 * qJD(5) - t103 * t11 - t105 * t139;
t96 = t100 * t108;
t93 = -pkin(5) * t105 - pkin(4);
t83 = t195 * t105;
t82 = t195 * t103;
t71 = t74 ^ 2;
t70 = (pkin(5) * t103 - t101) * t106;
t68 = t105 * t80;
t61 = t105 * t72;
t52 = -pkin(5) * t148 + t79;
t51 = t134 * t163;
t45 = -t101 * t174 + t68;
t41 = t101 * t164 - t157 + (-t103 * t164 + t143) * pkin(5);
t38 = -qJ(6) * t173 + t46;
t37 = -t141 + t64;
t36 = -qJ(6) * t171 + t68 + (-t101 * t103 + pkin(5)) * t104;
t34 = -t71 + t202;
t27 = (t153 + t180) * qJD(1) + t129;
t26 = (t153 - t180) * qJD(1) + t129;
t25 = -t151 + (-t154 + (t74 + t158) * t106) * qJD(1);
t24 = t151 + (t154 + (t74 - t158) * t106) * qJD(1);
t21 = -t105 * t145 + t61 + (-qJD(5) * t80 - t206) * t103;
t20 = -t103 * t145 + t138;
t18 = t103 * t120 - t177;
t17 = t76 * t183 - t178;
t16 = t74 * t143 + (-t74 * t164 + t181) * t103;
t15 = -t76 * t144 + (-t76 * t164 - t182) * t105;
t14 = -t90 * t143 - t43 * t104 + (-t122 * t103 - t106 * t74) * qJD(4);
t13 = -t90 * t144 - t42 * t104 + (t122 * t105 + t180) * qJD(4);
t12 = -qJ(6) * t143 + t114 * t103 + t138;
t9 = pkin(5) * t163 + t61 + t114 * t105 + ((qJ(6) * t106 - t80) * qJD(5) - t206) * t103;
t8 = -t181 - t132 * t183 + (t104 * t74 - t134 * t173) * qJD(4);
t7 = t182 + t132 * t187 + (-t90 * t171 + (t76 - t147) * t104) * qJD(4);
t5 = t209 * t103 + t29 * t105;
t4 = -t208 * t103 + t210 * t105;
t2 = (t105 * t74 + t188) * t164 + (t178 - t177 + (t103 * t74 - t105 * t76) * qJD(5)) * t106;
t1 = (qJD(1) * t76 + t136 * t104 - t74 * t163) * t105 + (qJD(1) * t74 + t137 * t104 + t76 * t163) * t103;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, qJ(2) * t97, 0, 0, 0, 0, 0, 0, 0, t97, t152, t94 * qJD(2) + t88 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t102) * qJD(1), -0.2e1 * t130, 0.2e1 * (-t100 + t99) * t155, -t107 * t104, 0.2e1 * t130, -t107 * t106, 0, t119 * t104 + t127 * t163, t119 * t106 - t127 * t164, -t176 * t97, t135 * qJD(3) + (qJD(1) * t101 + t87) * t203, t15, t2, t13, t16, t14, t51, t21 * t90 + (t11 + (t101 * t74 - t189) * qJD(4)) * t104 + (t64 * t160 - qJD(2) * t74 - t101 * t43 + t190 + (qJD(1) * t45 + t32) * qJD(4)) * t106, -t20 * t90 + (t139 + (t101 * t76 - t184) * qJD(4)) * t104 + (-t64 * t162 - qJD(2) * t76 + t101 * t42 + t185 + (-qJD(1) * t46 - t33) * qJD(4)) * t106, t106 * t111 + t124 * t164 - t20 * t74 - t21 * t76 + t42 * t45 - t43 * t46, -t64 * t157 - t139 * t46 + t11 * t45 + t20 * t33 + t21 * t32 + (-t106 * t54 + t64 * t164) * t101, t15, t2, t13, t16, t14, t51, t41 * t74 + t43 * t70 + t9 * t90 + (-t37 * t165 + t3) * t104 + (t37 * t160 + t191 + (qJD(1) * t36 + t19) * qJD(4)) * t106, -t12 * t90 + t41 * t76 - t42 * t70 + (-t37 * t158 - t6) * t104 + (-t37 * t162 + t186 + (-qJD(1) * t38 - t23) * qJD(4)) * t106, t106 * t112 - t12 * t74 + t126 * t164 + t36 * t42 - t38 * t43 - t76 * t9, t12 * t23 + t19 * t9 + t28 * t70 + t3 * t36 + t37 * t41 + t38 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t108 * qJ(2), 0, 0, 0, 0, 0, 0, 0, -t108, 0 (-qJD(3) - t94) * qJD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t91, 0.2e1 * t142, t96 + t179 (-t176 * t87 - qJD(3)) * qJD(1), 0, 0, 0, 0, 0, 0, t24, t27, t5 (t104 * t123 + t106 * t64) * qJD(1) + t111, 0, 0, 0, 0, 0, 0, t24, t27, t5 (t104 * t125 + t106 * t37) * qJD(1) + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, t133, 0, 0, 0, 0, 0, 0, t169 * t104, t169 * t106, 0 (-t88 + t203) * qJD(1), 0, 0, 0, 0, 0, 0, t8, t7, t1, -t124 * qJD(1) + (-qJD(4) * t123 - t54) * t106 + (qJD(4) * t64 + t110) * t104, 0, 0, 0, 0, 0, 0, t8, t7, t1, -t126 * qJD(1) + (-qJD(4) * t125 - t28) * t106 + (qJD(4) * t37 - qJD(5) * t126 - t103 * t3 + t105 * t6) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t96 - t179, 0, -t150, 0, 0, t106 * t133, -t170 * t168, 0, 0, t17, t4, t26, t18, t25, -t149, -t74 * t79 - pkin(4) * t43 - t185 - t39 * t90 + (-pkin(8) * t183 + t189) * qJD(5) + (t103 * t118 - t106 * t32) * qJD(1), -t76 * t79 + pkin(4) * t42 + t190 + t40 * t90 + (pkin(8) * t187 + t184) * qJD(5) + (t105 * t118 + t106 * t33) * qJD(1), t39 * t76 + t40 * t74 + (t136 * pkin(8) + t207) * t105 + (t137 * pkin(8) + t205) * t103, -pkin(4) * t54 + pkin(8) * t110 - t32 * t39 - t33 * t40 - t64 * t79, t17, t4, t26, t18, t25, -t149, -t186 + t43 * t93 - t52 * t74 - t193 * t90 + (t37 + t201) * t162 + (t37 * t174 + (qJD(4) * t82 - t19) * t106) * qJD(1), t191 - t42 * t93 - t52 * t76 + t192 * t90 + (pkin(5) * t188 + t105 * t37) * qJD(5) + (t37 * t172 + (qJD(4) * t83 + t23) * t106) * qJD(1), t42 * t82 + t43 * t83 + t193 * t76 + t192 * t74 + (-t19 * t90 + t6) * t105 + (-t3 - t200) * t103, t28 * t93 + t3 * t82 - t6 * t83 + (pkin(5) * t162 - t52) * t37 - t192 * t23 - t193 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t34, t29, -t197, -t209, t91, -t64 * t76 - t205, t64 * t74 - t207, 0, 0, t197, t34, t29, -t197, -t209, t91, 0.2e1 * t131 + t200 + (t141 - t37) * t76 + t109, -pkin(5) * t202 + t22 * t90 + (qJD(6) + t37) * t74 + t117, t42 * pkin(5) - t194 * t74, t194 * t23 + (-t37 * t76 + t3) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t210, -t71 - t202, t19 * t76 + t23 * t74 + t28;];
tauc_reg  = t10;
