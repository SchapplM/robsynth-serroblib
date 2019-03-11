% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPPRR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPPRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:19
% EndTime: 2019-03-09 01:30:22
% DurationCPUTime: 2.02s
% Computational Cost: add. (2412->317), mult. (3972->382), div. (0->0), fcn. (2190->10), ass. (0->173)
t90 = sin(qJ(5));
t93 = cos(qJ(5));
t132 = pkin(5) * t90 - pkin(8) * t93;
t183 = pkin(1) * qJDD(1);
t88 = cos(pkin(9));
t66 = t88 * t183;
t145 = -qJDD(1) * pkin(2) + qJDD(3) - t66;
t79 = qJDD(1) * qJ(4);
t135 = t79 - t145;
t133 = pkin(5) * t93 + pkin(8) * t90;
t41 = t133 * qJD(5) + qJD(4);
t12 = t41 * qJD(1) + t132 * qJDD(1) + t135;
t87 = sin(pkin(9));
t58 = pkin(1) * t87 + qJ(3);
t49 = t58 * qJD(1);
t46 = qJD(4) + t49;
t37 = -qJD(1) * pkin(7) + t46;
t25 = qJD(2) * t93 + t37 * t90;
t22 = qJD(5) * pkin(8) + t25;
t190 = pkin(2) + qJ(4);
t111 = -t132 - t190;
t209 = pkin(1) * t88;
t36 = -t111 + t209;
t23 = t36 * qJD(1) - qJD(3);
t89 = sin(qJ(6));
t92 = cos(qJ(6));
t5 = -t22 * t89 + t23 * t92;
t178 = qJD(5) * t93;
t65 = t87 * t183;
t186 = qJDD(1) * qJ(3) + t65;
t81 = qJD(3) * qJD(1);
t152 = -t81 - t186;
t136 = qJDD(4) - t152;
t33 = -qJDD(1) * pkin(7) + t136;
t155 = -t93 * qJDD(2) - t37 * t178 - t90 * t33;
t161 = qJD(2) * qJD(5);
t9 = -t90 * t161 - t155;
t7 = qJDD(5) * pkin(8) + t9;
t1 = qJD(6) * t5 + t89 * t12 + t92 * t7;
t6 = t22 * t92 + t23 * t89;
t126 = t5 * t89 - t6 * t92;
t78 = qJ(1) + pkin(9);
t69 = sin(t78);
t70 = cos(t78);
t187 = g(1) * t69 - g(2) * t70;
t11 = t92 * t12;
t2 = -qJD(6) * t6 - t89 * t7 + t11;
t218 = t126 * qJD(6) - t1 * t89 - t2 * t92 - t187;
t162 = qJD(1) * qJD(5);
t143 = t90 * t162;
t163 = t93 * qJDD(1);
t181 = qJD(1) * t93;
t44 = qJD(5) * t89 + t92 * t181;
t176 = qJD(6) * t44;
t18 = -t92 * qJDD(5) + (-t143 + t163) * t89 + t176;
t54 = qJD(1) * t90 + qJD(6);
t200 = t44 * t54;
t217 = t200 - t18;
t216 = g(1) * t70 + g(2) * t69;
t64 = -pkin(2) - t209;
t53 = qJ(4) - t64;
t182 = qJD(1) * t53;
t38 = -qJD(3) + t182;
t125 = -t38 * qJD(1) - t216;
t169 = t92 * qJD(5);
t172 = qJD(6) * t93;
t109 = t90 * t169 + t89 * t172;
t17 = t109 * qJD(1) - qJD(6) * t169 - t89 * qJDD(5) - t92 * t163;
t189 = -t17 * t90 + t44 * t178;
t42 = t89 * t181 - t169;
t117 = t42 * t54;
t215 = -t17 + t117;
t171 = t25 * qJD(5);
t27 = t93 * t33;
t10 = -t90 * qJDD(2) - t171 + t27;
t8 = -qJDD(5) * pkin(5) - t10;
t104 = g(3) * t90 - t216 * t93 - t8;
t214 = -pkin(8) * qJD(6) * t54 + t104;
t213 = qJD(5) * t126 + t8;
t142 = t93 * t162;
t164 = t90 * qJDD(1);
t40 = qJDD(6) + t142 + t164;
t193 = t92 * t40;
t212 = -t109 * t54 + t93 * t193;
t52 = -pkin(7) + t58;
t211 = qJD(5) * (qJD(3) + t38 + t182) + qJDD(5) * t52;
t24 = -qJD(2) * t90 + t37 * t93;
t100 = -(t24 * t90 - t25 * t93) * qJD(5) + t10 * t93 + t9 * t90;
t206 = g(3) * t93;
t205 = t5 * t54;
t204 = t6 * t54;
t203 = t70 * pkin(7);
t202 = t42 * t93;
t201 = t44 * t42;
t199 = t44 * t93;
t83 = t90 ^ 2;
t96 = qJD(1) ^ 2;
t198 = t83 * t96;
t197 = t89 * t40;
t196 = t89 * t90;
t195 = t90 * t18;
t194 = t90 * t92;
t192 = t93 * t17;
t191 = t93 * t18;
t84 = t93 ^ 2;
t185 = -t83 - t84;
t95 = qJD(5) ^ 2;
t184 = -t95 - t96;
t180 = qJD(5) * t42;
t179 = qJD(5) * t90;
t177 = qJD(6) * t42;
t175 = qJD(6) * t89;
t174 = qJD(6) * t90;
t173 = qJD(6) * t92;
t85 = qJDD(2) - g(3);
t168 = qJDD(1) * t53;
t166 = qJDD(5) * t90;
t165 = qJDD(5) * t93;
t160 = qJD(4) * qJD(1);
t159 = t54 * t196;
t158 = t54 * t194;
t156 = t93 * t96 * t90;
t94 = cos(qJ(1));
t153 = t94 * pkin(1) + t70 * pkin(2) + t69 * qJ(3);
t150 = t89 * t179;
t149 = t42 * t179;
t148 = t44 * t179;
t147 = t44 * t172;
t91 = sin(qJ(1));
t144 = -t91 * pkin(1) + t70 * qJ(3);
t141 = t70 * qJ(4) + t153;
t140 = -t17 + t177;
t139 = t185 * qJDD(1);
t137 = qJD(1) + t174;
t134 = t90 * t142;
t129 = g(1) * t91 - g(2) * t94;
t128 = -t52 * t174 + t41;
t127 = t5 * t92 + t6 * t89;
t123 = t24 * t93 + t25 * t90;
t120 = t145 - t187;
t34 = t135 + t160;
t119 = t38 * qJD(4) + t34 * t53;
t115 = -qJD(6) * t23 + t206 - t7;
t114 = -t79 + t120;
t113 = t54 * t173 + t197;
t112 = -t54 * t175 + t193;
t110 = -t190 * t69 + t144;
t108 = qJDD(1) * t58 + t186 - t216 + 0.2e1 * t81;
t21 = -qJD(5) * pkin(5) - t24;
t107 = -pkin(8) * t40 + t54 * t21;
t106 = -t216 * t90 - t206;
t105 = qJD(3) * t90 + qJD(6) * t36 + t52 * t178;
t103 = t161 - t125;
t102 = -t127 * qJD(6) + t1 * t92 - t2 * t89;
t99 = -t52 * t95 + t160 + t168 + t187 + t34;
t98 = qJD(5) * t21 + t102;
t75 = t84 * t96;
t48 = -t90 * t95 + t165;
t47 = -t93 * t95 - t166;
t45 = t133 * qJD(1);
t39 = t54 * t150;
t31 = t70 * t194 - t69 * t89;
t30 = -t70 * t196 - t69 * t92;
t29 = -t69 * t194 - t70 * t89;
t28 = t69 * t196 - t70 * t92;
t20 = t52 * t194 + t36 * t89;
t19 = -t52 * t196 + t36 * t92;
t15 = t24 * t92 + t45 * t89;
t14 = -t24 * t89 + t45 * t92;
t13 = t92 * t191;
t4 = -t105 * t89 + t128 * t92;
t3 = t105 * t92 + t128 * t89;
t16 = [0, 0, 0, 0, 0, qJDD(1), t129, g(1) * t94 + g(2) * t91, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t66 + t187, -0.2e1 * t65 + t216, 0 (t129 + (t87 ^ 2 + t88 ^ 2) * t183) * pkin(1), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(1) * t64 + t120, t108, -t152 * t58 + t49 * qJD(3) + t145 * t64 - g(1) * (-pkin(2) * t69 + t144) - g(2) * t153, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(4) + t108, -t114 + 0.2e1 * t160 + t168, -g(1) * t110 - g(2) * t141 + t46 * qJD(3) + t136 * t58 + t119, qJDD(1) * t84 - 0.2e1 * t134, -0.2e1 * t90 * t163 + 0.2e1 * (t83 - t84) * t162, t48, qJDD(1) * t83 + 0.2e1 * t134, t47, 0, t211 * t93 + t99 * t90, -t211 * t90 + t99 * t93, t52 * t139 + t185 * t81 - t100 + t216, -g(1) * (t110 - t203) - g(2) * (-t69 * pkin(7) + t141) + t123 * qJD(3) + t100 * t52 + t119, -t109 * t44 - t92 * t192, -t13 + (-t147 + t149) * t92 + (t148 + (t17 + t177) * t93) * t89, t189 + t212, t89 * t191 + (t92 * t172 - t150) * t42, -t195 + t39 + (-t113 - t180) * t93, t54 * t178 + t40 * t90, -g(1) * t29 - g(2) * t31 + t19 * t40 + t4 * t54 + (t2 + (-t21 * t89 + t42 * t52) * qJD(5)) * t90 + (-qJD(3) * t42 + qJD(5) * t5 + t21 * t173 - t18 * t52 + t8 * t89) * t93, -g(1) * t28 - g(2) * t30 - t20 * t40 - t3 * t54 + (-t1 + (-t21 * t92 + t44 * t52) * qJD(5)) * t90 + (-qJD(3) * t44 - qJD(5) * t6 + t17 * t52 - t21 * t175 + t8 * t92) * t93, t127 * t179 + t19 * t17 - t20 * t18 + t218 * t93 - t3 * t42 - t4 * t44, t1 * t20 + t6 * t3 + t2 * t19 + t5 * t4 - t8 * t93 * t52 - g(1) * (t144 - t203) - g(2) * (t132 * t70 + t141) + (-qJD(3) * t93 + t52 * t179) * t21 + (g(2) * pkin(7) - g(1) * t111) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, t47, -t48, 0, -qJD(5) * t123 - t10 * t90 + t9 * t93 - g(3), 0, 0, 0, 0, 0, 0, t195 + t39 + (-t113 + t180) * t93, t189 - t212, -t13 + (t147 + t149) * t92 + (t140 * t93 - t148) * t89, t213 * t90 + t98 * t93 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t96, -qJD(1) * t49 + t120, 0, 0, 0, 0, 0, 0, 0, -t96, -qJDD(1) (-qJD(4) - t46) * qJD(1) + t114, 0, 0, 0, 0, 0, 0, -0.2e1 * t142 - t164, 0.2e1 * t143 - t163, t75 + t198 (-qJD(4) - t123) * qJD(1) + t114, 0, 0, 0, 0, 0, 0 (t159 + t202) * qJD(1) - t112 (t158 + t199) * qJD(1) + t113, t215 * t92 - t217 * t89 (t126 * t90 + t21 * t93) * qJD(1) + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t96, t125 + t136, 0, 0, 0, 0, 0, 0, t184 * t90 + t165, t184 * t93 - t166, t139, t100 + t125, 0, 0, 0, 0, 0, 0, -t191 + (t180 - t197) * t90 + (-t137 * t92 - t89 * t178) * t54, t192 + (qJD(5) * t44 - t193) * t90 + (t137 * t89 - t93 * t169) * t54 (t137 * t44 - t42 * t178 - t195) * t92 + (t137 * t42 + t189) * t89, -t127 * qJD(1) - t213 * t93 + t98 * t90 - t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t75 - t198, t163, -t156, -t164, qJDD(5), t171 + t27 + (-qJD(5) * t37 - t85) * t90 - t103 * t93, t24 * qJD(5) + t103 * t90 + t155 + t206, 0, 0, -t17 * t89 + t92 * t200 (-t17 - t117) * t92 + (-t18 - t200) * t89 (t158 - t199) * qJD(1) + t113, t117 * t89 - t18 * t92 (-t159 + t202) * qJD(1) + t112, -t54 * t181, -pkin(5) * t18 + t107 * t89 - t14 * t54 - t5 * t181 + t214 * t92 - t25 * t42, pkin(5) * t17 + t107 * t92 + t15 * t54 + t6 * t181 - t214 * t89 - t25 * t44, t14 * t44 + t15 * t42 + (t1 - t205 + (-t18 + t176) * pkin(8)) * t92 + (t140 * pkin(8) - t2 - t204) * t89 + t106, -t5 * t14 - t6 * t15 - t21 * t25 + t104 * pkin(5) + (t102 + t106) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, -t42 ^ 2 + t44 ^ 2, t215, -t201, t217, t40, -g(1) * t30 + g(2) * t28 + t115 * t89 - t22 * t173 - t21 * t44 + t11 + t204, g(1) * t31 - g(2) * t29 + t21 * t42 + t205 + (qJD(6) * t22 - t12) * t89 + t115 * t92, 0, 0;];
tau_reg  = t16;
