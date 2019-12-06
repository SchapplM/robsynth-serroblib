% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:58
% EndTime: 2019-12-05 16:28:10
% DurationCPUTime: 2.38s
% Computational Cost: add. (7436->306), mult. (16132->450), div. (0->0), fcn. (11698->12), ass. (0->178)
t150 = sin(pkin(10));
t152 = cos(pkin(10));
t160 = cos(qJ(3));
t189 = qJD(2) * t160;
t157 = sin(qJ(3));
t190 = qJD(2) * t157;
t116 = t150 * t190 - t152 * t189;
t118 = t150 * t189 + t152 * t190;
t96 = t118 * t116;
t210 = qJDD(3) - t96;
t216 = t150 * t210;
t215 = t152 * t210;
t156 = sin(qJ(5));
t159 = cos(qJ(5));
t101 = -t159 * qJD(3) + t118 * t156;
t103 = qJD(3) * t156 + t118 * t159;
t79 = t103 * t101;
t181 = qJD(2) * qJD(3);
t175 = t160 * t181;
t180 = t157 * qJDD(2);
t123 = t175 + t180;
t142 = t160 * qJDD(2);
t176 = t157 * t181;
t124 = t142 - t176;
t170 = t123 * t150 - t152 * t124;
t95 = qJDD(5) + t170;
t211 = -t79 + t95;
t214 = t156 * t211;
t213 = t159 * t211;
t151 = sin(pkin(5));
t153 = cos(pkin(5));
t196 = sin(pkin(9));
t197 = cos(pkin(9));
t167 = t196 * g(1) - t197 * g(2);
t191 = -g(3) + qJDD(1);
t212 = t151 * t191 + t153 * t167;
t187 = qJD(3) * t118;
t80 = t170 + t187;
t163 = qJD(2) ^ 2;
t138 = t157 * t163 * t160;
t130 = qJDD(3) + t138;
t106 = -t151 * t167 + t153 * t191;
t128 = -t197 * g(1) - t196 * g(2);
t158 = sin(qJ(2));
t161 = cos(qJ(2));
t89 = t161 * t128 + t212 * t158;
t165 = -t163 * pkin(2) + qJDD(2) * pkin(7) + t89;
t69 = -t160 * t106 + t157 * t165;
t164 = -t69 + (-t123 + t175) * qJ(4) + t130 * pkin(3);
t129 = qJD(3) * pkin(3) - qJ(4) * t190;
t147 = t160 ^ 2;
t145 = t147 * t163;
t70 = t157 * t106 + t160 * t165;
t59 = -pkin(3) * t145 + t124 * qJ(4) - qJD(3) * t129 + t70;
t30 = -0.2e1 * qJD(4) * t116 + t150 * t164 + t152 * t59;
t113 = qJD(5) + t116;
t98 = t123 * t152 + t124 * t150;
t173 = -t159 * qJDD(3) + t156 * t98;
t48 = (qJD(5) - t113) * t103 + t173;
t99 = t101 ^ 2;
t100 = t103 ^ 2;
t112 = t113 ^ 2;
t114 = t116 ^ 2;
t115 = t118 ^ 2;
t209 = 0.2e1 * qJD(4);
t208 = pkin(4) * t150;
t169 = t128 * t158 - t212 * t161;
t84 = -qJDD(2) * pkin(2) - t163 * pkin(7) + t169;
t64 = -t124 * pkin(3) - qJ(4) * t145 + t129 * t190 + qJDD(4) + t84;
t206 = t150 * t64;
t93 = qJDD(3) + t96;
t205 = t150 * t93;
t204 = t152 * t64;
t203 = t152 * t93;
t162 = qJD(3) ^ 2;
t174 = t150 * t59 - t152 * t164;
t90 = pkin(4) * t116 - pkin(8) * t118;
t22 = -qJDD(3) * pkin(4) - t162 * pkin(8) + (t209 + t90) * t118 + t174;
t202 = t156 * t22;
t62 = t79 + t95;
t201 = t156 * t62;
t29 = t118 * t209 + t174;
t16 = t150 * t30 - t152 * t29;
t200 = t157 * t16;
t199 = t159 * t22;
t198 = t159 * t62;
t195 = t113 * t156;
t194 = t113 * t159;
t193 = t157 * t130;
t131 = qJDD(3) - t138;
t192 = t160 * t131;
t188 = qJD(3) * t116;
t186 = qJD(3) * t150;
t185 = qJD(3) * t152;
t182 = qJD(5) + t113;
t179 = t150 * t79;
t178 = t152 * t79;
t177 = -pkin(4) * t152 - pkin(3);
t23 = -pkin(4) * t162 + qJDD(3) * pkin(8) - t116 * t90 + t30;
t171 = -t98 + t188;
t36 = t80 * pkin(4) + t171 * pkin(8) + t64;
t14 = t156 * t23 - t159 * t36;
t15 = t156 * t36 + t159 * t23;
t6 = t14 * t156 + t159 * t15;
t17 = t150 * t29 + t152 * t30;
t40 = t157 * t69 + t160 * t70;
t2 = t150 * t6 - t152 * t22;
t3 = t150 * t22 + t152 * t6;
t1 = -t157 * t2 + t160 * t3;
t5 = -t14 * t159 + t15 * t156;
t168 = -qJDD(3) * t156 - t159 * t98;
t81 = -t170 + t187;
t73 = -qJD(5) * t101 - t168;
t146 = t157 ^ 2;
t143 = t146 * t163;
t137 = -t145 - t162;
t136 = -t143 - t162;
t127 = t143 + t145;
t126 = (t146 + t147) * qJDD(2);
t125 = t142 - 0.2e1 * t176;
t122 = 0.2e1 * t175 + t180;
t109 = -t115 - t162;
t108 = -t115 + t162;
t107 = t114 - t162;
t105 = -t136 * t157 - t192;
t104 = t137 * t160 - t193;
t91 = -t162 - t114;
t87 = t113 * t101;
t86 = -t100 + t112;
t85 = t99 - t112;
t83 = t98 + t188;
t78 = -t114 - t115;
t77 = t100 - t99;
t76 = -t100 - t112;
t75 = -t109 * t150 - t203;
t74 = t109 * t152 - t205;
t72 = -qJD(5) * t103 - t173;
t71 = -t112 - t99;
t68 = t99 + t100;
t66 = t152 * t91 - t216;
t65 = t150 * t91 + t215;
t60 = (-t101 * t159 + t103 * t156) * t113;
t58 = t150 * t83 + t152 * t81;
t57 = t150 * t81 - t152 * t83;
t53 = t182 * t101 + t168;
t52 = t73 + t87;
t51 = t73 - t87;
t50 = -t182 * t103 - t173;
t47 = -t103 * t195 + t159 * t73;
t46 = t101 * t194 - t156 * t72;
t45 = -t157 * t74 + t160 * t75;
t44 = t159 * t85 - t201;
t43 = -t156 * t86 + t213;
t42 = -t156 * t76 - t198;
t41 = t159 * t76 - t201;
t39 = t159 * t71 - t214;
t38 = t156 * t71 + t213;
t37 = -t157 * t65 + t160 * t66;
t34 = -t157 * t57 + t160 * t58;
t33 = -t156 * t51 + t159 * t50;
t32 = t156 * t52 - t159 * t48;
t31 = -t156 * t48 - t159 * t52;
t27 = -t150 * t53 + t152 * t42;
t26 = t150 * t42 + t152 * t53;
t25 = -t150 * t50 + t152 * t39;
t24 = t150 * t39 + t152 * t50;
t21 = -t150 * t68 + t152 * t32;
t20 = t150 * t32 + t152 * t68;
t19 = -pkin(8) * t41 + t199;
t18 = -pkin(8) * t38 + t202;
t12 = -t157 * t26 + t160 * t27;
t11 = -t157 * t24 + t160 * t25;
t10 = -t157 * t20 + t160 * t21;
t9 = -pkin(4) * t41 + t15;
t8 = -pkin(4) * t38 + t14;
t7 = t160 * t17 - t200;
t4 = -pkin(8) * t31 - t5;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t191, 0, 0, 0, 0, 0, 0, (qJDD(2) * t161 - t158 * t163) * t151, (-qJDD(2) * t158 - t161 * t163) * t151, 0, t106 * t153 + (t158 * t89 - t161 * t169) * t151, 0, 0, 0, 0, 0, 0, t153 * (t130 * t160 + t157 * t137) + (t158 * t104 + t161 * t125) * t151, t153 * (-t131 * t157 + t136 * t160) + (t158 * t105 - t161 * t122) * t151, (t126 * t158 + t127 * t161) * t151, t153 * (t157 * t70 - t160 * t69) + (t158 * t40 - t161 * t84) * t151, 0, 0, 0, 0, 0, 0, t153 * (t157 * t66 + t160 * t65) + (t158 * t37 - t161 * t80) * t151, t153 * (t157 * t75 + t160 * t74) + (t158 * t45 + t161 * t171) * t151, t153 * (t157 * t58 + t160 * t57) + (t158 * t34 - t161 * t78) * t151, t153 * (t157 * t17 + t16 * t160) + (t158 * t7 - t161 * t64) * t151, 0, 0, 0, 0, 0, 0, t153 * (t157 * t25 + t160 * t24) + (t158 * t11 - t161 * t38) * t151, t153 * (t157 * t27 + t160 * t26) + (t158 * t12 - t161 * t41) * t151, t153 * (t157 * t21 + t160 * t20) + (t10 * t158 - t161 * t31) * t151, t153 * (t157 * t3 + t160 * t2) + (t1 * t158 - t161 * t5) * t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t169, -t89, 0, 0, (t123 + t175) * t157, t122 * t160 + t125 * t157, t193 + t160 * (-t143 + t162), (t124 - t176) * t160, t157 * (t145 - t162) + t192, 0, pkin(2) * t125 + pkin(7) * t104 - t160 * t84, -pkin(2) * t122 + pkin(7) * t105 + t157 * t84, pkin(2) * t127 + pkin(7) * t126 + t40, -pkin(2) * t84 + pkin(7) * t40, t157 * (-t118 * t186 + t152 * t98) + t160 * (t118 * t185 + t150 * t98), t157 * (t150 * t171 - t152 * t80) + t160 * (-t150 * t80 - t152 * t171), t157 * (-t108 * t150 + t215) + t160 * (t108 * t152 + t216), t157 * (t116 * t185 + t150 * t170) + t160 * (t116 * t186 - t152 * t170), t157 * (t107 * t152 - t205) + t160 * (t107 * t150 + t203), (t157 * (-t116 * t152 + t118 * t150) + t160 * (-t116 * t150 - t118 * t152)) * qJD(3), t157 * (-qJ(4) * t65 + t206) + t160 * (-pkin(3) * t80 + qJ(4) * t66 - t204) - pkin(2) * t80 + pkin(7) * t37, t157 * (-qJ(4) * t74 + t204) + t160 * (pkin(3) * t171 + qJ(4) * t75 + t206) + pkin(2) * t171 + pkin(7) * t45, t157 * (-qJ(4) * t57 - t16) + t160 * (-pkin(3) * t78 + qJ(4) * t58 + t17) - pkin(2) * t78 + pkin(7) * t34, -qJ(4) * t200 + t160 * (-pkin(3) * t64 + qJ(4) * t17) - pkin(2) * t64 + pkin(7) * t7, t157 * (t152 * t47 + t179) + t160 * (t150 * t47 - t178), t157 * (t150 * t77 + t152 * t33) + t160 * (t150 * t33 - t152 * t77), t157 * (t150 * t52 + t152 * t43) + t160 * (t150 * t43 - t152 * t52), t157 * (t152 * t46 - t179) + t160 * (t150 * t46 + t178), t157 * (-t150 * t48 + t152 * t44) + t160 * (t150 * t44 + t152 * t48), t157 * (t150 * t95 + t152 * t60) + t160 * (t150 * t60 - t152 * t95), t157 * (-qJ(4) * t24 - t150 * t8 + t152 * t18) + t160 * (-pkin(3) * t38 + qJ(4) * t25 + t150 * t18 + t152 * t8) - pkin(2) * t38 + pkin(7) * t11, t157 * (-qJ(4) * t26 - t150 * t9 + t152 * t19) + t160 * (-pkin(3) * t41 + qJ(4) * t27 + t150 * t19 + t152 * t9) - pkin(2) * t41 + pkin(7) * t12, t157 * (-qJ(4) * t20 + t152 * t4) + t160 * (qJ(4) * t21 + t150 * t4) + pkin(7) * t10 + (t157 * t208 + t160 * t177 - pkin(2)) * t31, (t157 * (-pkin(8) * t152 + t208) + t160 * (-pkin(8) * t150 + t177) - pkin(2)) * t5 + (pkin(7) + qJ(4)) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, t143 - t145, t180, t138, t142, qJDD(3), -t69, -t70, 0, 0, t96, t115 - t114, t83, -t96, t81, qJDD(3), pkin(3) * t65 - t29, pkin(3) * t74 - t30, pkin(3) * t57, pkin(3) * t16, t103 * t194 + t156 * t73, t156 * t50 + t159 * t51, t159 * t86 + t214, t101 * t195 + t159 * t72, t156 * t85 + t198, (-t101 * t156 - t103 * t159) * t113, pkin(3) * t24 + pkin(4) * t50 + pkin(8) * t39 - t199, pkin(3) * t26 + pkin(4) * t53 + pkin(8) * t42 + t202, pkin(3) * t20 + pkin(4) * t68 + pkin(8) * t32 + t6, pkin(3) * t2 - pkin(4) * t22 + pkin(8) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t171, t78, t64, 0, 0, 0, 0, 0, 0, t38, t41, t31, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t77, t52, -t79, -t48, t95, -t14, -t15, 0, 0;];
tauJ_reg = t13;
