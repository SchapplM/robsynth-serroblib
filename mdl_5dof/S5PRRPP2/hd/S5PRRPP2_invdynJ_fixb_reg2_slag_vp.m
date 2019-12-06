% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:16
% EndTime: 2019-12-05 16:10:20
% DurationCPUTime: 1.92s
% Computational Cost: add. (1693->298), mult. (3857->372), div. (0->0), fcn. (2748->10), ass. (0->159)
t101 = sin(pkin(8));
t103 = cos(pkin(8));
t106 = sin(qJ(3));
t108 = cos(qJ(3));
t198 = qJ(4) + pkin(6);
t151 = qJD(3) * t198;
t126 = -t106 * qJD(4) - t108 * t151;
t109 = cos(qJ(2));
t181 = t103 * t108;
t68 = t101 * t106 - t181;
t131 = t68 * t109;
t55 = t108 * qJD(4) - t106 * t151;
t196 = qJD(1) * t131 + t101 * t126 + t103 * t55;
t107 = sin(qJ(2));
t102 = sin(pkin(7));
t104 = cos(pkin(7));
t144 = g(1) * t104 + g(2) * t102;
t135 = t144 * t107;
t96 = g(3) * t109;
t217 = t135 - t96;
t134 = t144 * t109;
t200 = g(3) * t107;
t122 = t134 + t200;
t100 = t108 ^ 2;
t99 = t106 ^ 2;
t189 = t100 + t99;
t150 = t109 * t189;
t174 = qJD(1) * t107;
t82 = qJD(2) * pkin(6) + t174;
t169 = t109 * qJD(1);
t194 = qJD(2) * pkin(2);
t83 = -t169 - t194;
t216 = t83 * t107 + t82 * t150;
t69 = t101 * t108 + t103 * t106;
t211 = t69 * qJD(2);
t155 = t198 * t106;
t76 = t198 * t108;
t38 = -t101 * t155 + t103 * t76;
t97 = qJ(3) + pkin(8);
t94 = sin(t97);
t215 = -qJDD(3) * t38 - (qJD(1) * t211 + t144 * t94) * t107 + t94 * t96;
t206 = t211 ^ 2;
t157 = qJD(2) * t181;
t173 = qJD(2) * t106;
t59 = t101 * t173 - t157;
t56 = t59 ^ 2;
t213 = -t56 - t206;
t212 = -t56 + t206;
t61 = t69 * qJD(3);
t170 = qJD(3) * t106;
t210 = -qJD(3) * t181 + t101 * t170;
t163 = t108 * qJDD(2);
t164 = t106 * qJDD(2);
t141 = t101 * t164 - t103 * t163;
t33 = qJD(2) * t61 + t141;
t167 = qJD(2) * qJD(3);
t153 = t106 * t167;
t124 = t69 * qJDD(2) - t101 * t153;
t152 = t108 * t167;
t34 = t103 * t152 + t124;
t209 = pkin(4) * t33 - qJ(5) * t34;
t110 = qJD(3) ^ 2;
t162 = t109 * qJDD(1);
t168 = qJD(1) * qJD(2);
t147 = t107 * t168 - t162;
t186 = qJDD(2) * pkin(2);
t66 = t147 - t186;
t156 = -t66 - t96;
t207 = -pkin(6) * t110 + (t144 + t168) * t107 + t156 + t186;
t204 = pkin(3) * t106;
t203 = pkin(3) * t108;
t199 = t211 * t59;
t67 = qJDD(2) * pkin(6) + t107 * qJDD(1) + t109 * t168;
t132 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + t67;
t149 = qJ(4) * qJD(2) + t82;
t139 = qJD(3) * t149;
t15 = qJDD(3) * pkin(3) - t132 * t106 - t108 * t139;
t20 = -t106 * t139 + t132 * t108;
t4 = -t101 * t20 + t103 * t15;
t5 = t101 * t15 + t103 * t20;
t197 = t101 * t55 - t103 * t126 - t69 * t169;
t54 = t149 * t108;
t40 = t103 * t54;
t53 = t149 * t106;
t42 = qJD(3) * pkin(3) - t53;
t24 = t101 * t42 + t40;
t193 = t101 * t54;
t95 = cos(t97);
t192 = t107 * t95;
t190 = t100 - t99;
t25 = -t101 * t53 + t40;
t187 = qJD(3) * t25;
t185 = qJDD(3) * pkin(4);
t183 = t102 * t109;
t180 = t104 * t108;
t179 = t104 * t109;
t178 = t106 * t109;
t26 = -t103 * t53 - t193;
t177 = qJD(5) - t26;
t176 = qJDD(1) - g(3);
t111 = qJD(2) ^ 2;
t175 = t110 + t111;
t172 = qJD(2) * t107;
t166 = qJDD(2) * t109;
t165 = qJDD(3) * t106;
t161 = pkin(3) * t170;
t160 = t106 * t111 * t108;
t93 = pkin(2) + t203;
t154 = t189 * t67;
t148 = t189 * qJDD(2);
t146 = t106 * t152;
t145 = pkin(4) * t95 + qJ(5) * t94;
t143 = g(1) * t102 - g(2) * t104;
t142 = t33 * t68 + t59 * t61;
t23 = t103 * t42 - t193;
t138 = qJD(3) * t61 + qJDD(3) * t68;
t27 = t210 * t107 - t109 * t211;
t28 = -qJD(2) * t131 - t107 * t61;
t51 = t69 * t107;
t52 = t68 * t107;
t133 = -t211 * t27 - t28 * t59 + t52 * t33 + t34 * t51;
t47 = t104 * t95 + t94 * t183;
t49 = -t102 * t95 + t94 * t179;
t129 = g(1) * t49 + g(2) * t47 + t94 * t200 + t4;
t37 = t101 * t76 + t103 * t155;
t127 = -qJDD(3) * t37 + t144 * t192 - t95 * t96;
t123 = qJD(3) * t27 - qJDD(3) * t51 - t109 * t33 + t59 * t172;
t65 = -t93 * qJD(2) + qJD(4) - t169;
t120 = -t210 * t59 + t211 * t61 + t33 * t69 + t34 * t68;
t48 = -t104 * t94 + t95 * t183;
t50 = t102 * t94 + t95 * t179;
t119 = g(1) * t50 + g(2) * t48 + g(3) * t192 - t5;
t21 = pkin(4) * t59 - qJ(5) * t211 + t65;
t118 = -t21 * t211 - qJDD(5) + t129;
t117 = -pkin(6) * qJDD(3) + (t169 + t83 - t194) * qJD(3);
t36 = pkin(3) * t153 - t93 * qJDD(2) + qJDD(4) + t147;
t116 = qJD(3) * t28 - qJDD(3) * t52 + t109 * t34 - t172 * t211;
t115 = -t83 * qJD(2) + t122 - t67;
t114 = 0.2e1 * t211 * qJD(3) + t141;
t113 = -t196 * t59 + t197 * t211 - t38 * t33 + t34 * t37 - t122;
t112 = t36 - t217;
t98 = qJDD(3) * qJ(5);
t91 = -pkin(3) * t103 - pkin(4);
t89 = pkin(3) * t101 + qJ(5);
t87 = t102 * t203;
t79 = t109 * t93;
t35 = -qJD(3) * t210 + qJDD(3) * t69;
t32 = pkin(4) * t68 - qJ(5) * t69 - t93;
t29 = pkin(3) * t173 + pkin(4) * t211 + qJ(5) * t59;
t22 = qJD(3) * qJ(5) + t24;
t19 = -qJD(3) * pkin(4) + qJD(5) - t23;
t18 = (t59 + t157) * qJD(3) + t124;
t17 = (-t59 + t157) * qJD(3) + t124;
t12 = pkin(4) * t61 + qJ(5) * t210 - qJD(5) * t69 + t161;
t6 = -t210 * t211 + t34 * t69;
t3 = qJDD(5) - t185 - t4;
t2 = -qJD(5) * t211 + t209 + t36;
t1 = qJD(3) * qJD(5) + t5 + t98;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t176, 0, 0, 0, 0, 0, 0, -t107 * t111 + t166, -qJDD(2) * t107 - t109 * t111, 0, -g(3) + (t107 ^ 2 + t109 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, (-0.2e1 * t153 + t163) * t109 + (-t175 * t108 - t165) * t107, (-qJDD(3) * t107 - 0.2e1 * t109 * t167) * t108 + (t175 * t107 - t166) * t106, t107 * t148 + t111 * t150, t216 * qJD(2) + t107 * t154 - t109 * t66 - g(3), 0, 0, 0, 0, 0, 0, t123, -t116, t133, -t109 * t36 + t65 * t172 + t23 * t27 + t24 * t28 - t4 * t51 - t5 * t52 - g(3), 0, 0, 0, 0, 0, 0, t123, t133, t116, -t1 * t52 - t109 * t2 + t21 * t172 - t19 * t27 + t22 * t28 + t3 * t51 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t162 + t217, -t176 * t107 + t134, 0, 0, qJDD(2) * t99 + 0.2e1 * t146, 0.2e1 * t106 * t163 + 0.2e1 * t190 * t167, t108 * t110 + t165, qJDD(2) * t100 - 0.2e1 * t146, qJDD(3) * t108 - t106 * t110, 0, t117 * t106 + t207 * t108, -t207 * t106 + t117 * t108, -t200 + t154 + pkin(6) * t148 + (-t189 * t168 - t144) * t109, (t135 + t156) * pkin(2) + (t154 - t122) * pkin(6) - t216 * qJD(1), t6, -t120, t35, t142, -t138, 0, -t59 * t174 - t33 * t93 + t36 * t68 + t61 * t65 + (t59 * t204 - t197) * qJD(3) + t127, -t34 * t93 + t36 * t69 - t210 * t65 + (t204 * t211 - t196) * qJD(3) + t215, t210 * t23 - t24 * t61 - t4 * t69 - t5 * t68 + t113, t5 * t38 - t4 * t37 - t36 * t93 - g(3) * (t107 * t198 + t79) + (t161 - t174) * t65 + t196 * t24 - t197 * t23 + t144 * (t107 * t93 - t109 * t198), t6, t35, t120, 0, t138, t142, t2 * t68 + t21 * t61 + t32 * t33 + (t12 - t174) * t59 - t197 * qJD(3) + t127, -t1 * t68 - t19 * t210 - t22 * t61 + t3 * t69 + t113, t196 * qJD(3) - t12 * t211 - t2 * t69 + t21 * t210 - t32 * t34 - t215, -g(3) * t79 + t1 * t38 + t21 * t12 + t2 * t32 + t3 * t37 + t196 * t22 + t197 * t19 + (-g(3) * t145 - t144 * t198) * t109 + (-g(3) * t198 - t21 * qJD(1) + t144 * (t145 + t93)) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, -t190 * t111, t164, t160, t163, qJDD(3), t115 * t106 - t143 * t108, t143 * t106 + t115 * t108, 0, 0, t199, t212, t18, -t199, -t141, qJDD(3), t187 - t211 * t65 + (qJDD(3) * t103 - t59 * t173) * pkin(3) + t129, qJD(3) * t26 + t59 * t65 + (-qJDD(3) * t101 - t173 * t211) * pkin(3) + t119, (t24 - t25) * t211 + (-t23 + t26) * t59 + (-t101 * t33 - t103 * t34) * pkin(3), -g(1) * t87 + t23 * t25 - t24 * t26 + (g(2) * t180 + t5 * t101 + t4 * t103 + (-t65 * qJD(2) + t122) * t106) * pkin(3), t199, t18, -t212, qJDD(3), t141, -t199, t187 - t29 * t59 + (pkin(4) - t91) * qJDD(3) + t118, -t33 * t89 + t34 * t91 + (t22 - t25) * t211 + (t19 - t177) * t59, qJDD(3) * t89 - t21 * t59 + t29 * t211 + t98 + (0.2e1 * qJD(5) - t26) * qJD(3) - t119, t1 * t89 + t3 * t91 - t21 * t29 - t19 * t25 - g(1) * (-pkin(3) * t104 * t178 - pkin(4) * t49 + qJ(5) * t50 + t87) - g(2) * (-pkin(4) * t47 + qJ(5) * t48 + (-t102 * t178 - t180) * pkin(3)) + t177 * t22 - (-pkin(4) * t94 + qJ(5) * t95 - t204) * t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t17, t213, t211 * t23 + t24 * t59 + t112, 0, 0, 0, 0, 0, 0, t114, t213, -t17, t22 * t59 + (-qJD(5) - t19) * t211 + t112 + t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t199, t18, -t206 - t110, -qJD(3) * t22 - t118 - t185;];
tau_reg = t7;
