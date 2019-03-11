% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% tau_reg [6x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:50
% EndTime: 2019-03-09 02:08:54
% DurationCPUTime: 1.92s
% Computational Cost: add. (1840->315), mult. (3247->385), div. (0->0), fcn. (1837->6), ass. (0->167)
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t124 = t92 * pkin(4) - t95 * pkin(8);
t90 = pkin(1) + qJ(3);
t55 = t124 + t90;
t26 = qJD(1) * t55 - qJD(2);
t68 = qJ(2) * qJD(1) + qJD(3);
t61 = -pkin(7) * qJD(1) + t68;
t54 = t92 * t61;
t32 = qJD(4) * pkin(8) + t54;
t91 = sin(qJ(5));
t94 = cos(qJ(5));
t11 = t91 * t26 + t94 * t32;
t87 = qJDD(1) * pkin(1);
t150 = t87 - qJDD(2);
t81 = qJDD(1) * qJ(3);
t138 = -t81 - t150;
t125 = pkin(4) * t95 + pkin(8) * t92;
t45 = qJD(4) * t125 + qJD(3);
t15 = qJD(1) * t45 + qJDD(1) * t124 - t138;
t13 = t94 * t15;
t156 = t94 * qJD(4);
t161 = qJD(5) * t91;
t108 = t156 * t92 + t161 * t95;
t151 = t95 * qJDD(1);
t16 = qJD(1) * t108 - qJD(5) * t156 - t91 * qJDD(4) - t94 * t151;
t162 = qJD(4) * t95;
t82 = qJD(1) * qJD(2);
t83 = qJ(2) * qJDD(1);
t137 = qJDD(3) + t82 + t83;
t51 = -pkin(7) * qJDD(1) + t137;
t21 = qJDD(4) * pkin(8) + t162 * t61 + t92 * t51;
t149 = qJD(1) * qJD(4);
t134 = t95 * t149;
t152 = t92 * qJDD(1);
t44 = qJDD(5) + t134 + t152;
t158 = t91 * qJD(4);
t167 = qJD(1) * t95;
t49 = t167 * t94 + t158;
t1 = t44 * pkin(5) + t16 * qJ(6) - qJD(5) * t11 - t49 * qJD(6) - t91 * t21 + t13;
t157 = t92 * qJD(1);
t65 = qJD(5) + t157;
t10 = t94 * t26 - t91 * t32;
t7 = -t49 * qJ(6) + t10;
t6 = t65 * pkin(5) + t7;
t47 = t167 * t91 - t156;
t8 = -t47 * qJ(6) + t11;
t120 = t6 * t91 - t8 * t94;
t96 = cos(qJ(1));
t79 = g(2) * t96;
t93 = sin(qJ(1));
t172 = g(1) * t93 - t79;
t160 = qJD(5) * t94;
t145 = -t91 * t15 - t26 * t160 - t94 * t21;
t110 = -t161 * t32 - t145;
t135 = t92 * t149;
t17 = qJD(5) * t49 - t94 * qJDD(4) + (-t135 + t151) * t91;
t2 = -t17 * qJ(6) - t47 * qJD(6) + t110;
t213 = qJD(5) * t120 - t1 * t94 - t2 * t91 - t172;
t193 = t49 * t65;
t212 = t17 - t193;
t194 = t47 * t65;
t211 = -t16 + t194;
t210 = qJD(1) * t90;
t178 = t96 * t94;
t186 = t93 * t91;
t34 = t92 * t186 - t178;
t179 = t96 * t91;
t185 = t93 * t94;
t36 = -t92 * t179 - t185;
t209 = -g(1) * t36 + g(2) * t34;
t208 = -g(1) * t96 - g(2) * t93;
t198 = g(3) * t92;
t163 = qJD(4) * t92;
t130 = -qJDD(4) * pkin(4) + t61 * t163;
t181 = t95 * t51;
t20 = t130 - t181;
t207 = qJD(5) * pkin(8) * t65 - t208 * t95 - t198 + t20;
t62 = -qJD(2) + t210;
t89 = -pkin(7) + qJ(2);
t206 = qJD(4) * (qJD(2) + t62 + t210) + qJDD(4) * t89;
t205 = t49 ^ 2;
t76 = 0.2e1 * t82;
t204 = -t7 + t6;
t203 = pkin(5) * t91;
t197 = g(3) * t95;
t177 = qJ(6) + pkin(8);
t132 = qJD(5) * t177;
t168 = t94 * qJ(6);
t180 = t95 * t61;
t53 = t125 * qJD(1);
t39 = t94 * t53;
t196 = -t132 * t94 - t39 - (pkin(5) * t95 + t168 * t92) * qJD(1) + (-qJD(6) + t180) * t91;
t195 = t16 * t91;
t192 = t49 * t94;
t191 = t49 * t95;
t190 = t65 * t91;
t189 = t91 * t44;
t67 = t94 * pkin(5) + pkin(4);
t188 = t92 * t67;
t187 = t92 * t94;
t184 = t94 * t44;
t183 = t95 * t16;
t182 = t95 * t47;
t175 = t94 * t180 + t91 * t53;
t176 = t94 * qJD(6) - t175 + (-qJ(6) * t157 - t132) * t91;
t174 = t89 * t187 + t91 * t55;
t173 = t96 * pkin(1) + t93 * qJ(2);
t86 = t95 ^ 2;
t171 = t92 ^ 2 - t86;
t97 = qJD(4) ^ 2;
t98 = qJD(1) ^ 2;
t170 = -t97 - t98;
t169 = qJ(6) * t95;
t166 = qJD(2) * t92;
t165 = qJD(4) * t47;
t164 = qJD(4) * t49;
t33 = -qJD(4) * pkin(4) - t180;
t159 = t33 * qJD(5);
t154 = qJDD(4) * t92;
t153 = t90 * qJDD(1);
t148 = qJD(3) * qJD(1);
t147 = t65 * t187;
t146 = t92 * t190;
t143 = t96 * qJ(3) + t173;
t142 = t95 * t168;
t141 = t89 * t162;
t140 = t95 * t156;
t139 = qJDD(2) - t172;
t136 = -t89 + t203;
t133 = t65 * t89 + t32;
t131 = t89 * t140 + t55 * t160 + t94 * t166 + t91 * t45;
t129 = -qJD(5) * t26 - t21;
t128 = -0.2e1 * t134;
t127 = qJD(5) * t92 + qJD(1);
t126 = -t87 + t139;
t121 = t6 * t94 + t8 * t91;
t118 = t177 * t95 - t188;
t116 = -t81 + t126;
t113 = t160 * t65 + t189;
t112 = t161 * t65 - t184;
t111 = t208 + t76 + 0.2e1 * t83;
t107 = t62 * qJD(1) - t208;
t106 = t17 * pkin(5) + qJDD(6) + t130;
t105 = t107 - t51;
t104 = -pkin(8) * t44 + t33 * t65;
t103 = -qJD(6) * t95 + (qJ(6) * qJD(4) - qJD(5) * t89) * t92;
t102 = -qJD(1) * t121 + t208;
t101 = -qJD(5) * t121 - t1 * t91 + t2 * t94;
t52 = -t138 + t148;
t99 = -t89 * t97 + t148 + t153 + t172 + t52;
t75 = t96 * qJ(2);
t71 = qJDD(4) * t95;
t58 = t177 * t94;
t57 = t177 * t91;
t43 = t47 ^ 2;
t41 = t94 * t55;
t37 = t92 * t178 - t186;
t35 = -t92 * t185 - t179;
t30 = t94 * t45;
t23 = -t91 * t169 + t174;
t22 = t47 * pkin(5) + qJD(6) + t33;
t19 = -t142 + t41 + (-t89 * t91 + pkin(5)) * t92;
t5 = t106 - t181;
t4 = -qJD(5) * t142 + t103 * t91 + t131;
t3 = pkin(5) * t162 + t30 + t103 * t94 + (-t141 - t166 + (-t55 + t169) * qJD(5)) * t91;
t9 = [qJDD(1), t172, -t208, -0.2e1 * t87 + t139, t111, t150 * pkin(1) - g(1) * (-t93 * pkin(1) + t75) - g(2) * t173 + (t76 + t83) * qJ(2), qJDD(3) + t111, -t116 + 0.2e1 * t148 + t153, t52 * t90 + t62 * qJD(3) + t137 * qJ(2) + t68 * qJD(2) - g(1) * (-t90 * t93 + t75) - g(2) * t143, t86 * qJDD(1) + t128 * t92, 0.2e1 * t149 * t171 - 0.2e1 * t151 * t92, -t97 * t92 + t71, -t97 * t95 - t154, 0, t206 * t95 + t99 * t92, -t206 * t92 + t99 * t95, -t108 * t49 - t94 * t183 (t47 * t94 + t49 * t91) * t163 + (t195 - t17 * t94 + (t47 * t91 - t192) * qJD(5)) * t95 (-t156 * t65 - t16) * t92 + (-t112 + t164) * t95 (t158 * t65 - t17) * t92 + (-t113 - t165) * t95, t162 * t65 + t44 * t92, -g(1) * t35 - g(2) * t37 + t30 * t65 + t41 * t44 + (-t133 * t160 + t165 * t89 + t13) * t92 + (-qJD(2) * t47 + t10 * qJD(4) + t159 * t94 - t89 * t17) * t95 + ((-qJD(5) * t55 - t141) * t65 + t20 * t95 + (-qJD(2) * t65 - t33 * qJD(4) - t89 * t44 + t129) * t92) * t91, -t131 * t65 - t174 * t44 - g(1) * t34 - g(2) * t36 + (t133 * t161 + (-t33 * t94 + t49 * t89) * qJD(4) + t145) * t92 + (-qJD(2) * t49 - t11 * qJD(4) - t159 * t91 + t89 * t16 + t20 * t94) * t95, t121 * t163 + t19 * t16 - t23 * t17 + t213 * t95 - t3 * t49 - t4 * t47, t2 * t23 + t8 * t4 + t1 * t19 + t6 * t3 - g(1) * (-pkin(5) * t179 - t96 * pkin(7) + t75) - g(2) * (t96 * t188 + t143) - t22 * t136 * t163 + (t5 * t136 + t22 * (pkin(5) * t160 - qJD(2)) + t177 * t79) * t95 + (-g(1) * (t118 - t90) - g(2) * (-pkin(7) - t203)) * t93; 0, 0, 0, qJDD(1), -t98, -t98 * qJ(2) + t126, -t98, -qJDD(1) (-qJD(3) - t68) * qJD(1) + t116, 0, 0, 0, 0, 0, t128 - t152, 0.2e1 * t135 - t151, 0, 0, 0, 0, 0 (t146 + t182) * qJD(1) + t112 (t147 + t191) * qJD(1) + t113, t211 * t94 + t212 * t91 (t120 * t92 + t22 * t95) * qJD(1) + t213; 0, 0, 0, 0, 0, 0, qJDD(1), -t98, -t107 + t137, 0, 0, 0, 0, 0, t170 * t92 + t71, t170 * t95 - t154, 0, 0, 0, 0, 0, -t95 * t17 + (t165 - t189) * t92 + (-t127 * t94 - t158 * t95) * t65, t183 + (t164 - t184) * t92 + (t127 * t91 - t140) * t65 (t127 * t49 - t162 * t47 - t17 * t92) * t94 + (t127 * t47 - t16 * t92 + t162 * t49) * t91 (-qJD(4) * t120 - t5) * t95 + (qJD(4) * t22 + t101) * t92 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t98 * t92, -t171 * t98, t151, -t152, qJDD(4), -t105 * t95 + t198, t105 * t92 + t197, t65 * t192 - t195 (-t16 - t194) * t94 + (-t17 - t193) * t91 (t147 - t191) * qJD(1) + t113 (-t146 + t182) * qJD(1) - t112, -t65 * t167, -t10 * t167 - t47 * t54 - pkin(4) * t17 - t39 * t65 + (t65 * t180 + t104) * t91 - t207 * t94, pkin(4) * t16 + t104 * t94 + t11 * t167 + t175 * t65 + t207 * t91 - t49 * t54, t102 * t92 - t57 * t16 - t58 * t17 - t176 * t47 - t196 * t49 + t101 - t197, t2 * t58 - t1 * t57 - t5 * t67 - g(3) * t118 + t176 * t8 + t196 * t6 + (pkin(5) * t190 - t54) * t22 + t208 * (t177 * t92 + t67 * t95); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t47, -t43 + t205, t211, -t212, t44, -t32 * t160 + t11 * t65 - t33 * t49 + t13 + (t129 + t197) * t91 + t209, g(1) * t37 - g(2) * t35 + t10 * t65 + t94 * t197 + t33 * t47 - t110, pkin(5) * t16 - t204 * t47, t204 * t8 + (t91 * t197 - t22 * t49 + t1 + t209) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 - t205, -t198 + t8 * t47 + t6 * t49 + (-t208 - t51) * t95 + t106;];
tau_reg  = t9;
