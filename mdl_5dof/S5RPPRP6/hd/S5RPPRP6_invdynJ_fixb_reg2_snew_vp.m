% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:16
% EndTime: 2019-12-31 17:55:22
% DurationCPUTime: 1.84s
% Computational Cost: add. (2245->188), mult. (5176->228), div. (0->0), fcn. (3265->6), ass. (0->124)
t111 = sin(pkin(7));
t112 = cos(pkin(7));
t115 = cos(qJ(4));
t113 = sin(qJ(4));
t146 = t112 * t113;
t129 = t111 * t115 + t146;
t91 = t129 * qJD(1);
t145 = t112 * t115;
t93 = (-t111 * t113 + t145) * qJD(1);
t160 = t93 * t91;
t188 = qJDD(4) + t160;
t157 = t113 * t188;
t117 = qJD(4) ^ 2;
t88 = t93 ^ 2;
t175 = -t88 - t117;
t26 = -t115 * t175 + t157;
t151 = t115 * t188;
t28 = t113 * t175 + t151;
t15 = t111 * t28 + t112 * t26;
t159 = -qJ(3) - pkin(1);
t206 = t159 * t15;
t205 = pkin(6) * t26;
t204 = pkin(6) * t28;
t166 = t91 ^ 2;
t75 = t166 - t117;
t203 = t111 * (t113 * t75 + t151) + t112 * (-t115 * t75 + t157);
t189 = qJDD(4) - t160;
t156 = t113 * t189;
t47 = t115 * t189;
t76 = -t88 + t117;
t200 = -t111 * (t115 * t76 + t156) + t112 * (-t113 * t76 + t47);
t171 = -t166 - t117;
t176 = t115 * t171 - t156;
t179 = t113 * t171 + t47;
t186 = t111 * t176 + t112 * t179;
t144 = t93 * qJD(4);
t170 = t129 * qJDD(1);
t62 = t170 + 0.2e1 * t144;
t197 = qJ(2) * t62 + t159 * t186;
t174 = -t88 - t166;
t138 = t112 * qJDD(1);
t139 = t111 * qJDD(1);
t90 = -t113 * t139 + t115 * t138;
t177 = t113 * t90 - t115 * t170;
t180 = -t113 * t170 - t115 * t90;
t185 = t111 * t177 + t112 * t180;
t196 = qJ(2) * t174 + t159 * t185;
t195 = pkin(6) * t176;
t194 = pkin(6) * t179;
t193 = pkin(6) * t180;
t148 = qJD(4) * t91;
t65 = -t148 + t90;
t39 = t148 - t65;
t192 = qJ(5) * t39;
t187 = -pkin(3) * t174 + pkin(6) * t177;
t118 = qJD(1) ^ 2;
t114 = sin(qJ(1));
t116 = cos(qJ(1));
t135 = t114 * g(1) - t116 * g(2);
t130 = qJDD(2) - t135;
t126 = -t118 * qJ(2) + t130;
t133 = -0.2e1 * qJD(1) * qJD(3) + t159 * qJDD(1) + t126;
t105 = t111 ^ 2;
t106 = t112 ^ 2;
t143 = t105 + t106;
t178 = pkin(3) * t139 - (t143 * pkin(6) - t159) * t118;
t173 = t88 - t166;
t172 = t143 * t118;
t63 = -t170 - t144;
t169 = -t63 * pkin(4) + t192;
t165 = 2 * qJD(5);
t164 = pkin(3) * t118;
t163 = pkin(4) * t115;
t162 = t111 * g(3);
t124 = t162 + (-pkin(6) * qJDD(1) - t111 * t164 + t133) * t112;
t54 = -t112 * g(3) + t133 * t111;
t43 = -pkin(6) * t139 - t105 * t164 + t54;
t21 = t113 * t43 - t115 * t124;
t22 = t113 * t124 + t115 * t43;
t7 = t113 * t22 - t115 * t21;
t161 = t112 * t7;
t108 = qJDD(1) * qJ(2);
t131 = t116 * g(1) + t114 * g(2);
t128 = -t108 + t131;
t127 = -qJDD(3) + t128;
t140 = qJD(2) * qJD(1);
t125 = t127 - 0.2e1 * t140;
t44 = t125 - t178;
t158 = t113 * t44;
t155 = t113 * t62;
t152 = t115 * t44;
t150 = t115 * t62;
t149 = qJ(5) * t115;
t147 = qJDD(1) * pkin(1);
t142 = qJD(4) * t113;
t141 = qJD(4) * t115;
t136 = -qJ(5) * t113 - pkin(3);
t8 = t113 * t21 + t115 * t22;
t70 = -t159 * t118 + t125;
t134 = -t70 + t108;
t51 = t91 * pkin(4) - t93 * qJ(5);
t132 = qJDD(4) * qJ(5) + qJD(4) * t165 - t91 * t51 + t22;
t12 = -t117 * pkin(4) + t132;
t13 = -qJDD(4) * pkin(4) - t117 * qJ(5) + t93 * t51 + qJDD(5) + t21;
t1 = t111 * (t113 * t13 + t115 * t12) + t112 * (t113 * t12 - t115 * t13);
t24 = t111 * t54 + t112 * (t133 * t112 + t162);
t123 = t112 * (-t113 * t63 + t91 * t141) - t111 * (t115 * t63 + t91 * t142);
t72 = t93 * t142;
t122 = t112 * t72 + (-t91 * t145 - t111 * (-t113 * t91 - t115 * t93)) * qJD(4);
t121 = t93 * t165 - t169 + t44;
t103 = 0.2e1 * t140;
t96 = t143 * qJDD(1);
t95 = t111 * t172;
t94 = t112 * t172;
t85 = -t126 + t147;
t64 = -0.2e1 * t148 + t90;
t19 = t112 * (t115 * t65 - t72) - t111 * (t113 * t65 + t93 * t141);
t11 = t103 + (pkin(4) * qJD(4) - (2 * qJD(5))) * t93 - t127 + t169 + t178;
t10 = (-t62 - t144) * pkin(4) + t121;
t9 = -qJ(5) * t174 + t13;
t6 = (-t117 - t174) * pkin(4) + t132;
t5 = -pkin(4) * t144 + t121 - t192;
t2 = t111 * t8 + t161;
t3 = [0, 0, 0, 0, 0, qJDD(1), t135, t131, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t130 - 0.2e1 * t147, t103 + 0.2e1 * t108 - t131, pkin(1) * t85 + qJ(2) * (-t118 * pkin(1) + t103 - t128), t106 * qJDD(1), -0.2e1 * t111 * t138, 0, t105 * qJDD(1), 0, 0, t134 * t111 - t159 * t95, t134 * t112 - t159 * t94, -qJ(2) * t172 - t159 * t96 - t24, -qJ(2) * t70 + t159 * t24, t19, t112 * (-t113 * t64 - t150) - t111 * (t115 * t64 - t155), t200, t123, -t203, t122, t112 * (-t158 - t194) - t111 * (-pkin(3) * t62 + t152 + t195) + t197, t112 * (-t152 + t205) - t111 * (-pkin(3) * t64 - t158 - t204) + qJ(2) * t64 - t206, t112 * (-t7 - t193) - t111 * (t187 + t8) + t196, -pkin(6) * t161 - t111 * (pkin(3) * t44 + pkin(6) * t8) - qJ(2) * t44 + t159 * t2, t19, t200, t112 * (-t113 * t39 + t150) - t111 * (t115 * t39 + t155), t122, t203, t123, t112 * (-t113 * t10 - t62 * t149 - t194) - t111 * (t115 * t10 + t136 * t62 + t195) + t197, t112 * (-t113 * t6 + t115 * t9 - t193) - t111 * (t113 * t9 + t115 * t6 + t187) + t196, t112 * (t115 * t5 - t205) - t111 * (t113 * t5 + t204) + (pkin(4) * t146 - t111 * (-pkin(3) - t163) + qJ(2)) * t39 + t206, (t112 * (pkin(4) * t113 - t149) - t111 * (t136 - t163) + qJ(2)) * t11 + (t159 - pkin(6)) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t118, -t85, 0, 0, 0, 0, 0, 0, -t95, -t94, -t96, t24, 0, 0, 0, 0, 0, 0, t186, -t15, t185, t2, 0, 0, 0, 0, 0, 0, t186, t185, t15, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t138, -t172, -t70, 0, 0, 0, 0, 0, 0, t62, t64, t174, -t44, 0, 0, 0, 0, 0, 0, t62, t174, t39, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t173, t90, -t160, -t170, qJDD(4), -t21, -t22, 0, 0, t160, t65 + t148, -t173, qJDD(4), t170, -t160, pkin(4) * t189 + qJ(5) * t171 - t13, -pkin(4) * t90 - qJ(5) * t170, qJ(5) * t188 + (-t117 - t175) * pkin(4) + t132, -pkin(4) * t13 + qJ(5) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, t90, t175, t13;];
tauJ_reg = t3;
