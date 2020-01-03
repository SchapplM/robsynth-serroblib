% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:24
% EndTime: 2019-12-31 17:49:28
% DurationCPUTime: 1.61s
% Computational Cost: add. (2678->202), mult. (6049->265), div. (0->0), fcn. (3913->8), ass. (0->124)
t126 = cos(qJ(4));
t125 = sin(qJ(4));
t123 = cos(pkin(8));
t121 = sin(pkin(8));
t157 = t121 * t125;
t92 = (-t123 * t126 + t157) * qJD(1);
t156 = t121 * t126;
t140 = t123 * t125 + t156;
t94 = t140 * qJD(1);
t169 = t94 * t92;
t178 = qJDD(4) + t169;
t166 = t125 * t178;
t127 = qJD(4) ^ 2;
t89 = t94 ^ 2;
t183 = -t89 - t127;
t26 = -t126 * t183 + t166;
t216 = pkin(6) * t26;
t160 = t126 * t178;
t28 = t125 * t183 + t160;
t215 = pkin(6) * t28;
t214 = t121 * t28 + t123 * t26;
t21 = t121 * t26 - t123 * t28;
t176 = t92 ^ 2;
t77 = t176 - t127;
t213 = t121 * (-t126 * t77 + t166) - t123 * (t125 * t77 + t160);
t86 = qJD(4) * t92;
t91 = t140 * qJDD(1);
t70 = t91 - t86;
t194 = t86 - t70;
t212 = t194 * qJ(5);
t122 = sin(pkin(7));
t124 = cos(pkin(7));
t179 = qJDD(4) - t169;
t165 = t125 * t179;
t180 = -t176 - t127;
t185 = t126 * t180 - t165;
t53 = t126 * t179;
t187 = t125 * t180 + t53;
t199 = -t121 * t187 + t123 * t185;
t155 = t94 * qJD(4);
t149 = t123 * qJDD(1);
t150 = t121 * qJDD(1);
t90 = t125 * t150 - t126 * t149;
t67 = t90 + 0.2e1 * t155;
t209 = pkin(1) * (t122 * t199 - t124 * t67) + qJ(3) * t199 - pkin(2) * t67;
t186 = t125 * t91 - t126 * t90;
t188 = -t125 * t90 - t126 * t91;
t197 = -t121 * t188 + t123 * t186;
t47 = t89 + t176;
t208 = pkin(2) * t47 + pkin(1) * (t122 * t197 + t124 * t47) + qJ(3) * t197;
t206 = pkin(6) * t185;
t205 = pkin(6) * t187;
t204 = pkin(6) * t188;
t172 = sin(qJ(1));
t173 = cos(qJ(1));
t139 = t173 * g(1) + t172 * g(2);
t175 = qJD(1) ^ 2;
t100 = -t175 * pkin(1) - t139;
t138 = t172 * g(1) - t173 * g(2);
t137 = qJDD(1) * pkin(1) + t138;
t168 = t124 * t100 + t122 * t137;
t193 = -t175 * pkin(2) + qJDD(1) * qJ(3) + 0.2e1 * qJD(1) * qJD(3) + t168;
t128 = t121 ^ 2;
t130 = t123 ^ 2;
t102 = (t128 + t130) * t175;
t200 = pkin(3) * t47 + pkin(6) * t186;
t198 = t121 * t185 + t123 * t187;
t196 = t121 * t186 + t123 * t188;
t78 = -t89 + t127;
t195 = t121 * (-t125 * t78 + t53) + t123 * (t126 * t78 + t165);
t146 = pkin(1) * t122 + qJ(3);
t147 = pkin(1) * t124 + pkin(2);
t144 = -t122 * t100 + t124 * t137;
t49 = -qJDD(1) * pkin(2) - t175 * qJ(3) + qJDD(3) - t144;
t184 = -t147 * qJDD(1) + t146 * t102 + t49;
t182 = t89 - t176;
t174 = 2 * qJD(5);
t171 = pkin(4) * t126;
t118 = -g(3) + qJDD(2);
t108 = t123 * t118;
t136 = t108 + (pkin(3) * t175 * t123 - pkin(6) * qJDD(1) - t193) * t121;
t152 = t130 * t175;
t42 = t121 * t118 + t193 * t123;
t34 = -pkin(3) * t152 + pkin(6) * t149 + t42;
t18 = t125 * t34 - t126 * t136;
t19 = t125 * t136 + t126 * t34;
t5 = t125 * t19 - t126 * t18;
t170 = t121 * t5;
t40 = -pkin(3) * t149 + t49 + (-t128 * t175 - t152) * pkin(6);
t167 = t125 * t40;
t164 = t125 * t67;
t161 = t126 * t40;
t159 = t126 * t67;
t158 = qJ(5) * t126;
t154 = qJD(4) * t125;
t153 = qJD(4) * t126;
t145 = -qJ(5) * t125 - pkin(3);
t41 = t193 * t121 - t108;
t23 = t121 * t41 + t123 * t42;
t6 = t125 * t18 + t126 * t19;
t57 = t92 * pkin(4) - t94 * qJ(5);
t143 = qJDD(4) * qJ(5) + qJD(4) * t174 - t92 * t57 + t19;
t12 = -qJDD(4) * pkin(4) - t127 * qJ(5) + t94 * t57 + qJDD(5) + t18;
t68 = -t90 - t155;
t135 = -t68 * pkin(4) + t212 + t40;
t134 = t94 * t174 - t135;
t133 = t121 * (-t125 * t68 + t92 * t153) + t123 * (t126 * t68 + t92 * t154);
t75 = t94 * t154;
t132 = t121 * t75 + (-t92 * t156 + t123 * (-t125 * t92 - t126 * t94)) * qJD(4);
t113 = t130 * qJDD(1);
t112 = t128 * qJDD(1);
t101 = t113 + t112;
t69 = t91 - 0.2e1 * t86;
t22 = t121 * (t126 * t70 - t75) + t123 * (t125 * t70 + t94 * t153);
t13 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t94 + t135;
t11 = -t127 * pkin(4) + t143;
t10 = (-t67 - t155) * pkin(4) + t134;
t9 = -pkin(4) * t155 + t134 - t212;
t8 = qJ(5) * t47 + t12;
t7 = (-t127 + t47) * pkin(4) + t143;
t4 = t126 * t11 + t125 * t12;
t3 = t125 * t11 - t126 * t12;
t2 = t123 * t6 - t170;
t1 = [0, 0, 0, 0, 0, qJDD(1), t138, t139, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t124 * qJDD(1) - t122 * t175) + t144, pkin(1) * (-t122 * qJDD(1) - t124 * t175) - t168, 0, pkin(1) * (t122 * t168 + t124 * t144), t112, 0.2e1 * t121 * t149, 0, t113, 0, 0, -t184 * t123, t184 * t121, pkin(2) * t102 + qJ(3) * t101 + pkin(1) * (t122 * t101 + t124 * t102) + t23, -pkin(2) * t49 + qJ(3) * t23 + pkin(1) * (t122 * t23 - t124 * t49), t22, t121 * (-t125 * t69 - t159) + t123 * (t126 * t69 - t164), t195, t133, -t213, t132, t121 * (t167 - t205) + t123 * (-pkin(3) * t67 - t161 + t206) + t209, t121 * (t161 + t216) + t123 * (-pkin(3) * t69 + t167 - t215) - pkin(2) * t69 + qJ(3) * t21 + pkin(1) * (t122 * t21 - t124 * t69), t121 * (-t5 - t204) + t123 * (t200 + t6) + t208, -pkin(6) * t170 + t123 * (-pkin(3) * t40 + pkin(6) * t6) - pkin(2) * t40 + qJ(3) * t2 + pkin(1) * (t122 * t2 - t124 * t40), t22, t195, t121 * (-t125 * t194 + t159) + t123 * (t126 * t194 + t164), t132, t213, t133, t121 * (-t125 * t10 - t67 * t158 - t205) + t123 * (t126 * t10 + t145 * t67 + t206) + t209, t121 * (-t125 * t7 + t126 * t8 - t204) + t123 * (t125 * t8 + t126 * t7 + t200) + t208, t121 * (t126 * t9 - t216) + t123 * (t125 * t9 + t215) - t146 * t21 - (-pkin(4) * t157 + t123 * (pkin(3) + t171) + t147) * t194, (t121 * (pkin(4) * t125 - t158) + t123 * (t145 - t171) - t147) * t13 + (t146 + pkin(6)) * (-t121 * t3 + t123 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t42 - t123 * t41, 0, 0, 0, 0, 0, 0, t198, -t214, t196, t121 * t6 + t123 * t5, 0, 0, 0, 0, 0, 0, t198, t196, t214, t121 * t4 + t123 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t150, -t102, t49, 0, 0, 0, 0, 0, 0, t67, t69, -t47, t40, 0, 0, 0, 0, 0, 0, t67, -t47, t194, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t182, t91, -t169, -t90, qJDD(4), -t18, -t19, 0, 0, t169, t86 + t70, -t182, qJDD(4), t90, -t169, pkin(4) * t179 + qJ(5) * t180 - t12, -pkin(4) * t91 - qJ(5) * t90, qJ(5) * t178 + (-t127 - t183) * pkin(4) + t143, -pkin(4) * t12 + qJ(5) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, t91, t183, t12;];
tauJ_reg = t1;
