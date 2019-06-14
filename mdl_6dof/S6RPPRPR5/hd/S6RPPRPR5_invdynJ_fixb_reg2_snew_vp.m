% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:21:55
% EndTime: 2019-05-05 14:22:01
% DurationCPUTime: 2.09s
% Computational Cost: add. (8493->283), mult. (17422->367), div. (0->0), fcn. (10487->8), ass. (0->184)
t153 = sin(pkin(9));
t154 = cos(pkin(9));
t160 = cos(qJ(4));
t191 = qJD(1) * t160;
t127 = -t154 * qJD(4) + t153 * t191;
t129 = t153 * qJD(4) + t154 * t191;
t107 = t129 * t127;
t188 = qJD(1) * qJD(4);
t178 = t160 * t188;
t157 = sin(qJ(4));
t187 = t157 * qJDD(1);
t133 = t178 + t187;
t210 = -t107 + t133;
t216 = t153 * t210;
t215 = t154 * t210;
t156 = sin(qJ(6));
t130 = qJDD(6) + t133;
t159 = cos(qJ(6));
t102 = t159 * t127 + t156 * t129;
t104 = -t156 * t127 + t159 * t129;
t77 = t104 * t102;
t211 = t130 - t77;
t214 = t156 * t211;
t213 = t159 * t211;
t179 = t157 * t188;
t186 = t160 * qJDD(1);
t134 = -t179 + t186;
t116 = t153 * qJDD(4) + t154 * t134;
t174 = -t154 * qJDD(4) + t153 * t134;
t66 = -t102 * qJD(6) + t159 * t116 - t156 * t174;
t190 = t157 * qJD(1);
t140 = qJD(6) + t190;
t89 = t140 * t102;
t212 = t66 - t89;
t175 = t156 * t116 + t159 * t174;
t53 = (qJD(6) - t140) * t104 + t175;
t100 = t102 ^ 2;
t101 = t104 ^ 2;
t125 = t127 ^ 2;
t126 = t129 ^ 2;
t139 = t140 ^ 2;
t163 = qJD(1) ^ 2;
t209 = 2 * qJD(3);
t208 = pkin(5) * t157;
t155 = t163 * pkin(7);
t158 = sin(qJ(1));
t161 = cos(qJ(1));
t177 = t158 * g(1) - t161 * g(2);
t172 = -qJDD(2) + t177;
t165 = t163 * qJ(2) + t172;
t205 = pkin(1) + qJ(3);
t176 = t205 * qJDD(1);
t164 = t176 + t165;
t78 = t133 * pkin(4) - t134 * qJ(5) - t155 + (t209 + (pkin(4) * t160 + qJ(5) * t157) * qJD(4)) * qJD(1) + t164;
t162 = qJD(4) ^ 2;
t171 = pkin(4) * t157 - qJ(5) * t160;
t168 = t163 * t171;
t150 = qJDD(1) * qJ(2);
t173 = t161 * g(1) + t158 * g(2);
t170 = 0.2e1 * qJD(2) * qJD(1) - t173;
t166 = qJDD(3) + t170;
t114 = -t205 * t163 + t150 + t166;
t109 = -qJDD(1) * pkin(7) + t114;
t99 = t160 * g(3) - t157 * t109;
t83 = -t162 * pkin(4) + qJDD(4) * qJ(5) - t157 * t168 - t99;
t38 = 0.2e1 * qJD(5) * t129 + t153 * t83 - t154 * t78;
t181 = t127 * t190;
t93 = -t116 - t181;
t33 = pkin(5) * t210 + t93 * pkin(8) - t38;
t169 = pkin(5) * t190 - t129 * pkin(8);
t39 = -0.2e1 * qJD(5) * t127 + t153 * t78 + t154 * t83;
t34 = -t125 * pkin(5) - pkin(8) * t174 - t169 * t190 + t39;
t14 = t156 * t34 - t159 * t33;
t15 = t156 * t33 + t159 * t34;
t7 = -t159 * t14 + t156 * t15;
t207 = t153 * t7;
t206 = t154 * t7;
t204 = qJ(2) - pkin(7);
t95 = t107 + t133;
t203 = t153 * t95;
t202 = t154 * t95;
t98 = t157 * g(3) + t160 * t109;
t82 = qJDD(4) * pkin(4) + t162 * qJ(5) - t160 * t168 - qJDD(5) + t98;
t47 = -pkin(5) * t174 + t125 * pkin(8) - t129 * t169 + t82;
t201 = t156 * t47;
t68 = t130 + t77;
t200 = t156 * t68;
t199 = t159 * t47;
t198 = t159 * t68;
t197 = t160 * t82;
t196 = qJDD(1) * pkin(1);
t195 = t140 * t156;
t194 = t140 * t159;
t152 = t160 ^ 2;
t193 = t152 * t163;
t182 = t157 * t163 * t160;
t192 = t157 * (qJDD(4) + t182);
t185 = qJD(1) * t209;
t184 = t157 * t77;
t183 = t157 * t107;
t180 = t129 * t190;
t8 = t156 * t14 + t159 * t15;
t23 = t153 * t38 + t154 * t39;
t22 = t153 * t39 - t154 * t38;
t72 = -t157 * t99 + t160 * t98;
t167 = t171 + t205;
t91 = -t174 + t180;
t113 = t164 + t185;
t151 = t157 ^ 2;
t146 = 0.2e1 * t150;
t145 = t151 * t163;
t137 = t145 + t193;
t136 = (-t151 - t152) * qJDD(1);
t135 = -0.2e1 * t179 + t186;
t132 = 0.2e1 * t178 + t187;
t131 = t160 * (qJDD(4) - t182);
t120 = t165 + t196;
t119 = -t126 - t145;
t118 = -t126 + t145;
t117 = t125 - t145;
t112 = -t192 + t160 * (-t162 - t193);
t111 = t157 * (-t145 - t162) + t131;
t108 = -t155 + t113;
t105 = -t145 - t125;
t92 = t116 - t181;
t90 = t174 + t180;
t87 = -t125 - t126;
t86 = -t101 + t139;
t85 = t100 - t139;
t84 = -t101 - t139;
t80 = -t153 * t119 - t202;
t79 = t154 * t119 - t203;
t76 = t101 - t100;
t73 = -t139 - t100;
t71 = t154 * t105 - t216;
t70 = t153 * t105 + t215;
t65 = -t104 * qJD(6) - t175;
t64 = -t153 * t93 + t154 * t91;
t63 = t153 * t91 + t154 * t93;
t62 = (-t102 * t159 + t104 * t156) * t140;
t61 = (-t102 * t156 - t104 * t159) * t140;
t60 = -t100 - t101;
t59 = t157 * t80 - t160 * t92;
t58 = t157 * t71 - t160 * t90;
t57 = t66 + t89;
t52 = (qJD(6) + t140) * t104 + t175;
t51 = t159 * t85 - t200;
t50 = -t156 * t86 + t213;
t49 = t156 * t85 + t198;
t48 = t159 * t86 + t214;
t46 = -t104 * t195 + t159 * t66;
t45 = t104 * t194 + t156 * t66;
t44 = t102 * t194 - t156 * t65;
t43 = t102 * t195 + t159 * t65;
t42 = t157 * t64 - t160 * t87;
t41 = -t156 * t84 - t198;
t40 = t159 * t84 - t200;
t36 = t159 * t73 - t214;
t35 = t156 * t73 + t213;
t31 = t156 * t57 - t159 * t53;
t30 = -t156 * t212 - t159 * t52;
t29 = -t156 * t53 - t159 * t57;
t28 = -t156 * t52 + t159 * t212;
t27 = -pkin(8) * t40 - t199;
t26 = -t153 * t40 + t154 * t41;
t25 = t153 * t41 + t154 * t40;
t24 = -pkin(8) * t35 - t201;
t21 = -t153 * t35 + t154 * t36;
t20 = t153 * t36 + t154 * t35;
t19 = t157 * t23 + t197;
t18 = -pkin(5) * t212 + pkin(8) * t41 - t201;
t17 = t157 * t26 - t160 * t212;
t16 = -pkin(5) * t52 + pkin(8) * t36 + t199;
t12 = t157 * t21 - t160 * t52;
t11 = -t153 * t29 + t154 * t31;
t10 = t153 * t31 + t154 * t29;
t9 = t157 * t11 - t160 * t60;
t6 = pkin(5) * t47 + pkin(8) * t8;
t5 = -pkin(8) * t29 - t7;
t4 = -pkin(5) * t60 + pkin(8) * t31 + t8;
t3 = t154 * t8 - t207;
t2 = t153 * t8 + t206;
t1 = t157 * t3 + t160 * t47;
t13 = [0, 0, 0, 0, 0, qJDD(1), t177, t173, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t172 - 0.2e1 * t196, t146 + t170, qJ(2) * (-t163 * pkin(1) + t150 + t170) + pkin(1) * t120, qJDD(1), 0, 0, 0, 0, 0, 0, t146 + t166, t172 + 0.2e1 * t176 + t185, qJ(2) * t114 + t205 * t113, (t134 - t179) * t160, -t160 * t132 - t157 * t135, t131 - t157 * (t162 - t193), (t133 + t178) * t157, t160 * (t145 - t162) - t192, 0, t157 * t108 + t204 * t111 + t205 * t132, t160 * t108 + t204 * t112 + t205 * t135, t204 * t136 - t205 * t137 - t72, t205 * t108 + t204 * t72, t160 * (t154 * t116 - t153 * t180) + t183, t160 * (-t153 * t92 - t154 * t90) - t157 * (-t126 + t125), t160 * (-t153 * t118 + t215) - t157 * t93, t160 * (t153 * t174 + t154 * t181) - t183, t160 * (t154 * t117 - t203) + t157 * t91, (t133 + (-t127 * t154 + t129 * t153) * t191) * t157, -t153 * t197 - t157 * t38 + t167 * t70 + t204 * t58, -t154 * t197 - t157 * t39 + t167 * t79 + t204 * t59, -t160 * t22 + t167 * t63 + t204 * t42, t167 * t22 + t204 * t19, t160 * (-t153 * t45 + t154 * t46) + t184, t160 * (-t153 * t28 + t154 * t30) + t157 * t76, t160 * (-t153 * t48 + t154 * t50) + t157 * t57, t160 * (-t153 * t43 + t154 * t44) - t184, t160 * (-t153 * t49 + t154 * t51) - t157 * t53, t160 * (-t153 * t61 + t154 * t62) + t157 * t130, t160 * (-t153 * t16 + t154 * t24) - t157 * (-pkin(5) * t35 + t14) + t167 * t20 + t204 * t12, t160 * (-t153 * t18 + t154 * t27) - t157 * (-pkin(5) * t40 + t15) + t167 * t25 + t204 * t17, t160 * (-t153 * t4 + t154 * t5) + t29 * t208 + t204 * t9 + t167 * t10, t160 * (-pkin(8) * t206 - t153 * t6) + t7 * t208 + t167 * t2 + t204 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t163, -t120, 0, 0, 0, 0, 0, 0, 0, -t163, -qJDD(1), -t113, 0, 0, 0, 0, 0, 0, -t132, -t135, t137, -t108, 0, 0, 0, 0, 0, 0, -t70, -t79, -t63, -t22, 0, 0, 0, 0, 0, 0, -t20, -t25, -t10, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t163, t114, 0, 0, 0, 0, 0, 0, t111, t112, t136, t72, 0, 0, 0, 0, 0, 0, t58, t59, t42, t19, 0, 0, 0, 0, 0, 0, t12, t17, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, -t145 + t193, t186, -t182, -t187, qJDD(4), t98, t99, 0, 0, t153 * t116 + t154 * t180, -t153 * t90 + t154 * t92, t154 * t118 + t216, t153 * t181 - t154 * t174, t153 * t117 + t202, (-t127 * t153 - t129 * t154) * t190, -pkin(4) * t90 + qJ(5) * t71 + t154 * t82, -pkin(4) * t92 + qJ(5) * t80 - t153 * t82, -pkin(4) * t87 + qJ(5) * t64 + t23, pkin(4) * t82 + qJ(5) * t23, t153 * t46 + t154 * t45, t153 * t30 + t154 * t28, t153 * t50 + t154 * t48, t153 * t44 + t154 * t43, t153 * t51 + t154 * t49, t153 * t62 + t154 * t61, -pkin(4) * t52 + qJ(5) * t21 + t153 * t24 + t154 * t16, -pkin(4) * t212 + qJ(5) * t26 + t153 * t27 + t154 * t18, -pkin(4) * t60 + qJ(5) * t11 + t153 * t5 + t154 * t4, pkin(4) * t47 - pkin(8) * t207 + qJ(5) * t3 + t154 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t92, t87, -t82, 0, 0, 0, 0, 0, 0, t52, t212, t60, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t76, t57, -t77, -t53, t130, -t14, -t15, 0, 0;];
tauJ_reg  = t13;
