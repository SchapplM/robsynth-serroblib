% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRP5
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:41
% EndTime: 2019-12-05 15:38:49
% DurationCPUTime: 1.73s
% Computational Cost: add. (2099->189), mult. (4901->250), div. (0->0), fcn. (3375->8), ass. (0->125)
t118 = sin(pkin(8));
t120 = cos(pkin(8));
t124 = cos(qJ(4));
t122 = sin(qJ(4));
t148 = qJD(2) * t120;
t153 = t118 * t122;
t93 = qJD(2) * t153 - t124 * t148;
t152 = t118 * t124;
t136 = t120 * t122 + t152;
t95 = t136 * qJD(2);
t164 = t95 * t93;
t170 = qJDD(4) + t164;
t162 = t122 * t170;
t126 = qJD(4) ^ 2;
t90 = t95 ^ 2;
t175 = -t90 - t126;
t26 = -t124 * t175 + t162;
t156 = t124 * t170;
t28 = t122 * t175 + t156;
t21 = t118 * t26 - t120 * t28;
t208 = qJ(3) * t21;
t123 = sin(qJ(2));
t207 = t123 * t21;
t206 = pkin(6) * t26;
t205 = pkin(6) * t28;
t169 = t93 ^ 2;
t78 = t169 - t126;
t204 = t118 * (-t124 * t78 + t162) - t120 * (t122 * t78 + t156);
t87 = t93 * qJD(4);
t92 = t136 * qJDD(2);
t66 = t92 - t87;
t185 = t87 - t66;
t203 = t185 * qJ(5);
t171 = qJDD(4) - t164;
t161 = t122 * t171;
t172 = -t169 - t126;
t176 = t124 * t172 - t161;
t49 = t124 * t171;
t178 = t122 * t172 + t49;
t188 = -t118 * t178 + t120 * t176;
t150 = t95 * qJD(4);
t143 = t120 * qJDD(2);
t144 = t118 * qJDD(2);
t91 = t122 * t144 - t124 * t143;
t63 = t91 + 0.2e1 * t150;
t200 = -pkin(2) * t63 + qJ(3) * t188;
t125 = cos(qJ(2));
t199 = t123 * t188 - t125 * t63;
t177 = t122 * t92 - t124 * t91;
t179 = -t122 * t91 - t124 * t92;
t187 = -t118 * t179 + t120 * t177;
t45 = t90 + t169;
t198 = pkin(2) * t45 + qJ(3) * t187;
t197 = t123 * t187 + t125 * t45;
t195 = pkin(6) * t176;
t194 = pkin(6) * t178;
t193 = pkin(6) * t179;
t189 = pkin(3) * t45 + pkin(6) * t177;
t79 = -t90 + t126;
t186 = t118 * (-t122 * t79 + t49) + t120 * (t124 * t79 + t161);
t168 = qJD(2) ^ 2;
t119 = sin(pkin(7));
t121 = cos(pkin(7));
t102 = -g(1) * t121 - g(2) * t119;
t149 = -g(3) + qJDD(1);
t76 = t102 * t125 + t123 * t149;
t72 = -pkin(2) * t168 + qJDD(2) * qJ(3) + t76;
t183 = 0.2e1 * qJD(2) * qJD(3) + t72;
t127 = t118 ^ 2;
t129 = t120 ^ 2;
t100 = (t127 + t129) * t168;
t174 = t90 - t169;
t167 = 2 * qJD(5);
t166 = pkin(4) * t124;
t101 = -g(1) * t119 + g(2) * t121;
t151 = t120 * t101;
t132 = t151 + (-pkin(6) * qJDD(2) - t72 + (pkin(3) * t148 - 0.2e1 * qJD(3)) * qJD(2)) * t118;
t141 = t118 * t101 + t120 * t183;
t145 = t129 * t168;
t34 = -pkin(3) * t145 + pkin(6) * t143 + t141;
t18 = t122 * t34 - t124 * t132;
t19 = t122 * t132 + t124 * t34;
t5 = t122 * t19 - t124 * t18;
t165 = t118 * t5;
t117 = qJDD(2) * pkin(2);
t75 = -t123 * t102 + t125 * t149;
t71 = -qJ(3) * t168 + qJDD(3) - t117 - t75;
t46 = -pkin(3) * t143 + t71 + (-t127 * t168 - t145) * pkin(6);
t163 = t122 * t46;
t160 = t122 * t63;
t157 = t124 * t46;
t155 = t124 * t63;
t154 = qJ(5) * t124;
t147 = qJD(4) * t122;
t146 = qJD(4) * t124;
t142 = t125 * qJDD(2);
t139 = -qJ(5) * t122 - pkin(3);
t23 = t118 * (t118 * t183 - t151) + t120 * t141;
t6 = t122 * t18 + t124 * t19;
t138 = -t71 + t117;
t53 = pkin(4) * t93 - qJ(5) * t95;
t137 = qJDD(4) * qJ(5) + qJD(4) * t167 - t53 * t93 + t19;
t11 = -pkin(4) * t126 + t137;
t12 = -qJDD(4) * pkin(4) - qJ(5) * t126 + t53 * t95 + qJDD(5) + t18;
t1 = -t118 * (t11 * t122 - t12 * t124) + t120 * (t11 * t124 + t12 * t122);
t64 = -t91 - t150;
t135 = -t64 * pkin(4) + t203 + t46;
t134 = t167 * t95 - t135;
t133 = t118 * (-t122 * t64 + t146 * t93) + t120 * (t124 * t64 + t147 * t93);
t74 = t95 * t147;
t131 = t118 * t74 + (-t93 * t152 + t120 * (-t122 * t93 - t124 * t95)) * qJD(4);
t111 = t129 * qJDD(2);
t110 = t127 * qJDD(2);
t99 = t111 + t110;
t97 = t120 * t100;
t96 = t118 * t100;
t65 = t92 - 0.2e1 * t87;
t22 = t118 * (t124 * t66 - t74) + t120 * (t122 * t66 + t146 * t95);
t13 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t95 + t135;
t10 = t134 + (-t63 - t150) * pkin(4);
t9 = -pkin(4) * t150 + t134 - t203;
t8 = qJ(5) * t45 + t12;
t7 = (-t126 + t45) * pkin(4) + t137;
t2 = t120 * t6 - t165;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t149, 0, 0, 0, 0, 0, 0, -t123 * t168 + t142, -qJDD(2) * t123 - t125 * t168, 0, t123 * t76 + t125 * t75, 0, 0, 0, 0, 0, 0, t120 * t142 - t123 * t97, -t118 * t142 + t123 * t96, t100 * t125 + t123 * t99, t123 * t23 - t125 * t71, 0, 0, 0, 0, 0, 0, t199, -t125 * t65 + t207, t197, t123 * t2 - t125 * t46, 0, 0, 0, 0, 0, 0, t199, t197, -t125 * t185 - t207, t1 * t123 - t125 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t75, -t76, 0, 0, t110, 0.2e1 * t118 * t143, 0, t111, 0, 0, -qJ(3) * t97 + t120 * t138, qJ(3) * t96 - t118 * t138, pkin(2) * t100 + qJ(3) * t99 + t23, -pkin(2) * t71 + qJ(3) * t23, t22, t118 * (-t122 * t65 - t155) + t120 * (t124 * t65 - t160), t186, t133, -t204, t131, t118 * (t163 - t194) + t120 * (-pkin(3) * t63 - t157 + t195) + t200, t118 * (t157 + t206) + t120 * (-pkin(3) * t65 + t163 - t205) - pkin(2) * t65 + t208, t118 * (-t5 - t193) + t120 * (t189 + t6) + t198, -pkin(6) * t165 + t120 * (-pkin(3) * t46 + pkin(6) * t6) - pkin(2) * t46 + qJ(3) * t2, t22, t186, t118 * (-t122 * t185 + t155) + t120 * (t124 * t185 + t160), t131, t204, t133, t118 * (-t10 * t122 - t154 * t63 - t194) + t120 * (t124 * t10 + t139 * t63 + t195) + t200, t118 * (-t122 * t7 + t124 * t8 - t193) + t120 * (t122 * t8 + t124 * t7 + t189) + t198, t118 * (t124 * t9 - t206) + t120 * (t122 * t9 + t205) - t208 - (-pkin(4) * t153 + t120 * (pkin(3) + t166) + pkin(2)) * t185, (t118 * (pkin(4) * t122 - t154) + t120 * (t139 - t166) - pkin(2)) * t13 + (qJ(3) + pkin(6)) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t144, -t100, t71, 0, 0, 0, 0, 0, 0, t63, t65, -t45, t46, 0, 0, 0, 0, 0, 0, t63, -t45, t185, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t174, t92, -t164, -t91, qJDD(4), -t18, -t19, 0, 0, t164, t87 + t66, -t174, qJDD(4), t91, -t164, pkin(4) * t171 + qJ(5) * t172 - t12, -pkin(4) * t92 - qJ(5) * t91, qJ(5) * t170 + (-t126 - t175) * pkin(4) + t137, -pkin(4) * t12 + qJ(5) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, t92, t175, t12;];
tauJ_reg = t3;
