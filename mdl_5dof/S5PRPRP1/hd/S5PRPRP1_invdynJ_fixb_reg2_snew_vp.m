% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRP1
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:40
% EndTime: 2019-12-05 15:28:46
% DurationCPUTime: 1.56s
% Computational Cost: add. (2015->179), mult. (4730->232), div. (0->0), fcn. (3267->8), ass. (0->119)
t107 = sin(pkin(8));
t108 = cos(pkin(8));
t111 = cos(qJ(4));
t109 = sin(qJ(4));
t140 = t107 * t109;
t81 = (-t108 * t111 + t140) * qJD(2);
t139 = t107 * t111;
t123 = t108 * t109 + t139;
t83 = t123 * qJD(2);
t154 = t83 * t81;
t161 = qJDD(4) + t154;
t152 = t109 * t161;
t112 = qJD(4) ^ 2;
t78 = t83 ^ 2;
t165 = -t78 - t112;
t17 = -t111 * t165 + t152;
t146 = t111 * t161;
t19 = t109 * t165 + t146;
t196 = qJ(3) * (t107 * t17 - t108 * t19);
t195 = pkin(6) * t17;
t194 = pkin(6) * t19;
t193 = t107 * t19 + t108 * t17;
t159 = t81 ^ 2;
t67 = t159 - t112;
t191 = t107 * (-t111 * t67 + t152) - t108 * (t109 * t67 + t146);
t75 = qJD(4) * t81;
t80 = t123 * qJDD(2);
t59 = t80 - t75;
t175 = t75 - t59;
t190 = t175 * qJ(5);
t162 = qJDD(4) - t154;
t151 = t109 * t162;
t163 = -t159 - t112;
t166 = t111 * t163 - t151;
t187 = pkin(6) * t166;
t41 = t111 * t162;
t168 = t109 * t163 + t41;
t186 = pkin(6) * t168;
t133 = t108 * qJDD(2);
t134 = t107 * qJDD(2);
t79 = t109 * t134 - t111 * t133;
t169 = -t109 * t79 - t111 * t80;
t185 = pkin(6) * t169;
t100 = t108 ^ 2;
t99 = t107 ^ 2;
t182 = t100 + t99;
t113 = qJD(2) ^ 2;
t110 = sin(qJ(2));
t143 = sin(pkin(7));
t144 = cos(pkin(7));
t122 = t143 * g(1) - t144 * g(2);
t157 = cos(qJ(2));
t86 = -t144 * g(1) - t143 * g(2);
t118 = -t110 * t122 - t157 * t86;
t173 = -t113 * pkin(2) + qJDD(2) * qJ(3) + 0.2e1 * qJD(2) * qJD(3) - t118;
t167 = t109 * t80 - t111 * t79;
t38 = t78 + t159;
t181 = pkin(3) * t38 + pkin(6) * t167;
t180 = t107 * t166 + t108 * t168;
t179 = t107 * t167 + t108 * t169;
t68 = -t78 + t112;
t178 = t107 * (-t109 * t68 + t41) + t108 * (t111 * t68 + t151);
t138 = t83 * qJD(4);
t56 = t79 + 0.2e1 * t138;
t177 = qJ(3) * (-t107 * t168 + t108 * t166) - pkin(2) * t56;
t176 = pkin(2) * t38 + qJ(3) * (-t107 * t169 + t108 * t167);
t164 = t78 - t159;
t103 = t113 * qJ(3);
t106 = qJDD(2) * pkin(2);
t127 = -t110 * t86 + t157 * t122;
t55 = qJDD(3) - t103 - t106 - t127;
t160 = t182 * t103 - t106 + t55;
t158 = 2 * qJD(5);
t156 = pkin(4) * t111;
t104 = -g(3) + qJDD(1);
t92 = t108 * t104;
t121 = t92 + (pkin(3) * t108 * t113 - pkin(6) * qJDD(2) - t173) * t107;
t141 = t100 * t113;
t36 = t107 * t104 + t173 * t108;
t25 = -pkin(3) * t141 + pkin(6) * t133 + t36;
t13 = t109 * t25 - t111 * t121;
t14 = t109 * t121 + t111 * t25;
t3 = t109 * t14 - t111 * t13;
t155 = t107 * t3;
t32 = -pkin(3) * t133 + t55 + (-t113 * t99 - t141) * pkin(6);
t153 = t109 * t32;
t150 = t109 * t56;
t147 = t111 * t32;
t145 = t111 * t56;
t142 = qJ(5) * t111;
t137 = qJD(4) * t109;
t136 = qJD(4) * t111;
t129 = -qJ(5) * t109 - pkin(3);
t35 = t173 * t107 - t92;
t128 = t107 * t35 + t108 * t36;
t4 = t109 * t13 + t111 * t14;
t45 = t81 * pkin(4) - t83 * qJ(5);
t125 = qJDD(4) * qJ(5) + qJD(4) * t158 - t81 * t45 + t14;
t10 = -qJDD(4) * pkin(4) - t112 * qJ(5) + t83 * t45 + qJDD(5) + t13;
t57 = -t79 - t138;
t120 = -t57 * pkin(4) + t190 + t32;
t119 = t83 * t158 - t120;
t117 = t107 * (-t109 * t57 + t81 * t136) + t108 * (t111 * t57 + t81 * t137);
t65 = t83 * t137;
t116 = t107 * t65 + (-t81 * t139 + t108 * (-t109 * t81 - t111 * t83)) * qJD(4);
t97 = t100 * qJDD(2);
t96 = t99 * qJDD(2);
t85 = t182 * t113;
t58 = t80 - 0.2e1 * t75;
t15 = t107 * (t111 * t59 - t65) + t108 * (t109 * t59 + t83 * t136);
t11 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t83 + t120;
t9 = -t112 * pkin(4) + t125;
t8 = t119 + (-t56 - t138) * pkin(4);
t7 = -pkin(4) * t138 + t119 - t190;
t6 = qJ(5) * t38 + t10;
t5 = (-t112 + t38) * pkin(4) + t125;
t2 = t109 * t10 + t111 * t9;
t1 = -t111 * t10 + t109 * t9;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107 * t36 - t108 * t35, 0, 0, 0, 0, 0, 0, t180, -t193, t179, t107 * t4 + t108 * t3, 0, 0, 0, 0, 0, 0, t180, t179, t193, t108 * t1 + t107 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t127, t118, 0, 0, t96, 0.2e1 * t107 * t133, 0, t97, 0, 0, -t160 * t108, t160 * t107, pkin(2) * t85 + qJ(3) * (t97 + t96) + t128, -pkin(2) * t55 + qJ(3) * t128, t15, t107 * (-t109 * t58 - t145) + t108 * (t111 * t58 - t150), t178, t117, -t191, t116, t107 * (t153 - t186) + t108 * (-pkin(3) * t56 - t147 + t187) + t177, t107 * (t147 + t195) + t108 * (-pkin(3) * t58 + t153 - t194) - pkin(2) * t58 + t196, t107 * (-t3 - t185) + t108 * (t181 + t4) + t176, -pkin(6) * t155 + t108 * (-pkin(3) * t32 + pkin(6) * t4) - pkin(2) * t32 + qJ(3) * (t108 * t4 - t155), t15, t178, t107 * (-t109 * t175 + t145) + t108 * (t111 * t175 + t150), t116, t191, t117, t107 * (-t109 * t8 - t56 * t142 - t186) + t108 * (t111 * t8 + t129 * t56 + t187) + t177, t107 * (-t109 * t5 + t111 * t6 - t185) + t108 * (t109 * t6 + t111 * t5 + t181) + t176, t107 * (t111 * t7 - t195) + t108 * (t109 * t7 + t194) - t196 - (-pkin(4) * t140 + t108 * (pkin(3) + t156) + pkin(2)) * t175, (t107 * (pkin(4) * t109 - t142) + t108 * (t129 - t156) - pkin(2)) * t11 + (qJ(3) + pkin(6)) * (-t107 * t1 + t108 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, t134, -t85, t55, 0, 0, 0, 0, 0, 0, t56, t58, -t38, t32, 0, 0, 0, 0, 0, 0, t56, -t38, t175, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t164, t80, -t154, -t79, qJDD(4), -t13, -t14, 0, 0, t154, t75 + t59, -t164, qJDD(4), t79, -t154, pkin(4) * t162 + qJ(5) * t163 - t10, -pkin(4) * t80 - qJ(5) * t79, qJ(5) * t161 + (-t112 - t165) * pkin(4) + t125, -pkin(4) * t10 + qJ(5) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, t80, t165, t10;];
tauJ_reg = t12;
