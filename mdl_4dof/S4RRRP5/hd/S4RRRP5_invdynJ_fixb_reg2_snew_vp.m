% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:05
% EndTime: 2019-12-31 17:17:10
% DurationCPUTime: 1.74s
% Computational Cost: add. (2435->212), mult. (5260->249), div. (0->0), fcn. (3369->6), ass. (0->133)
t120 = sin(qJ(2));
t123 = cos(qJ(2));
t122 = cos(qJ(3));
t119 = sin(qJ(3));
t115 = qJDD(2) + qJDD(3);
t149 = qJD(1) * t120;
t89 = -t122 * t123 * qJD(1) + t119 * t149;
t152 = t120 * t122;
t91 = (t123 * t119 + t152) * qJD(1);
t170 = t91 * t89;
t63 = -t170 - t115;
t163 = t119 * t63;
t116 = qJD(2) + qJD(3);
t177 = t116 ^ 2;
t88 = t91 ^ 2;
t181 = -t88 - t177;
t38 = t122 * t181 + t163;
t157 = t122 * t63;
t40 = -t119 * t181 + t157;
t210 = pkin(5) * (t120 * t38 - t123 * t40);
t209 = pkin(2) * t38;
t208 = pkin(6) * t38;
t207 = pkin(6) * t40;
t178 = t89 ^ 2;
t78 = t178 - t177;
t205 = t120 * (-t122 * t78 - t163) + t123 * (-t119 * t78 + t157);
t51 = t88 + t178;
t202 = pkin(1) * t51;
t201 = pkin(2) * t51;
t179 = -t177 - t178;
t180 = -t170 + t115;
t53 = t122 * t180;
t187 = t119 * t179 + t53;
t200 = pkin(2) * t187;
t162 = t119 * t180;
t185 = t122 * t179 - t162;
t199 = pkin(6) * t185;
t198 = pkin(6) * t187;
t167 = t116 * t89;
t109 = t120 * qJDD(1);
t147 = qJD(1) * qJD(2);
t145 = t123 * t147;
t97 = t109 + t145;
t110 = t123 * qJDD(1);
t146 = t120 * t147;
t98 = t110 - t146;
t137 = t119 * t98 + t122 * t97;
t49 = -t89 * qJD(3) + t137;
t183 = t49 - t167;
t141 = t119 * t97 - t122 * t98;
t133 = t91 * qJD(3) + t141;
t84 = t116 * t91;
t30 = t133 + t84;
t195 = t120 * (t119 * t183 + t122 * t30) - t123 * (-t119 * t30 + t122 * t183);
t79 = -t88 + t177;
t194 = t120 * (-t119 * t79 + t53) + t123 * (t122 * t79 + t162);
t193 = pkin(5) * (-t120 * t187 + t123 * t185);
t188 = t183 * qJ(4);
t182 = t49 + t167;
t65 = t88 - t178;
t176 = 2 * qJD(4);
t175 = t133 * pkin(3);
t174 = cos(qJ(1));
t173 = pkin(3) * t119;
t172 = pkin(3) * t122;
t125 = qJD(1) ^ 2;
t151 = t120 * t125;
t121 = sin(qJ(1));
t136 = t174 * g(1) + t121 * g(2);
t156 = qJDD(1) * pkin(5);
t93 = -t125 * pkin(1) - t136 + t156;
t160 = t120 * t93;
t128 = qJDD(2) * pkin(2) - t97 * pkin(6) - t160 + (pkin(2) * t151 + pkin(6) * t147 - g(3)) * t123;
t118 = t123 ^ 2;
t112 = t118 * t125;
t135 = qJD(2) * pkin(2) - pkin(6) * t149;
t76 = -t120 * g(3) + t123 * t93;
t45 = -pkin(2) * t112 + t98 * pkin(6) - qJD(2) * t135 + t76;
t20 = t119 * t45 - t122 * t128;
t21 = t119 * t128 + t122 * t45;
t5 = t119 * t21 - t122 * t20;
t171 = t120 * t5;
t64 = t89 * pkin(3) - t91 * qJ(4);
t139 = t115 * qJ(4) + t116 * t176 - t89 * t64 + t21;
t11 = -pkin(3) * t177 + t139;
t13 = -t115 * pkin(3) - qJ(4) * t177 + t91 * t64 + qJDD(4) + t20;
t169 = -pkin(3) * t13 + qJ(4) * t11;
t31 = t133 - t84;
t168 = -pkin(3) * t182 - qJ(4) * t31;
t165 = t119 * t182;
t143 = t121 * g(1) - t174 * g(2);
t134 = qJDD(1) * pkin(1) + t143;
t52 = t98 * pkin(2) - t135 * t149 + (pkin(6) * t118 + pkin(5)) * t125 + t134;
t164 = t119 * t52;
t158 = t122 * t52;
t155 = t116 * t119;
t154 = t116 * t122;
t104 = t123 * t151;
t153 = t120 * (qJDD(2) + t104);
t150 = t123 * (qJDD(2) - t104);
t148 = qJD(3) + t116;
t144 = -qJ(4) * t119 - pkin(2);
t6 = t119 * t20 + t122 * t21;
t75 = t123 * g(3) + t160;
t140 = t120 * t75 + t123 * t76;
t132 = -pkin(3) * t181 - qJ(4) * t63 + t11;
t131 = pkin(3) * t180 + qJ(4) * t179 - t13;
t130 = t120 * (t119 * t133 + t89 * t154) + t123 * (-t122 * t133 + t89 * t155);
t77 = t91 * t155;
t129 = t120 * t77 + (-t89 * t152 + t123 * (-t119 * t89 - t122 * t91)) * t116;
t127 = -pkin(3) * t84 + t91 * t176 + t52;
t126 = t127 + t188;
t124 = qJD(2) ^ 2;
t117 = t120 ^ 2;
t111 = t117 * t125;
t99 = t110 - 0.2e1 * t146;
t96 = t109 + 0.2e1 * t145;
t92 = t125 * pkin(5) + t134;
t33 = -t148 * t89 + t137;
t32 = (-qJD(3) + t116) * t91 - t141;
t29 = t148 * t91 + t141;
t26 = t122 * t182;
t18 = t122 * t32 + t165;
t17 = -t122 * t31 + t165;
t16 = t119 * t32 - t26;
t15 = -t119 * t31 - t26;
t14 = t120 * (t122 * t49 - t77) + t123 * (t119 * t49 + t91 * t154);
t8 = qJ(4) * t51 + t13;
t7 = (-t177 + t51) * pkin(3) + t139;
t4 = (-t29 - t133) * pkin(3) + t126;
t3 = t127 - t175 + 0.2e1 * t188;
t1 = t119 * t11 - t122 * t13;
t2 = [0, 0, 0, 0, 0, qJDD(1), t143, t136, 0, 0, (t97 + t145) * t120, t120 * t99 + t123 * t96, t153 + t123 * (-t111 + t124), (t98 - t146) * t123, t120 * (t112 - t124) + t150, 0, t123 * t92 + pkin(1) * t99 + pkin(5) * (t123 * (-t112 - t124) - t153), -t120 * t92 - pkin(1) * t96 + pkin(5) * (-t150 - t120 * (-t111 - t124)), pkin(1) * (t111 + t112) + (t117 + t118) * t156 + t140, pkin(1) * t92 + pkin(5) * t140, t14, -t195, t194, t130, -t205, t129, t120 * (-t164 - t198) + t123 * (-pkin(2) * t30 + t158 + t199) - pkin(1) * t30 + t193, t120 * (-t158 - t208) + t123 * (-pkin(2) * t33 - t164 + t207) - pkin(1) * t33 - t210, t120 * (-pkin(6) * t16 - t5) + t123 * (pkin(6) * t18 + t201 + t6) + t202 + pkin(5) * (-t120 * t16 + t123 * t18), -pkin(6) * t171 + t123 * (pkin(2) * t52 + pkin(6) * t6) + pkin(1) * t52 + pkin(5) * (t123 * t6 - t171), t14, t194, t195, t129, t205, t130, t120 * (-t119 * t4 - t198) + t123 * (t122 * t4 + t199) + t193 + (-qJ(4) * t152 + t123 * t144 - pkin(1)) * t29, t120 * (-pkin(6) * t15 - t119 * t7 + t122 * t8) + t123 * (pkin(6) * t17 + t119 * t8 + t122 * t7 + t201) + t202 + pkin(5) * (-t120 * t15 + t123 * t17), t120 * (t122 * t3 + t208) + t123 * (t119 * t3 - t207) + t210 + (-t120 * t173 + t123 * (pkin(2) + t172) + pkin(1)) * t183, (t120 * (qJ(4) * t122 - t173) + t123 * (-t144 + t172) + pkin(1)) * (t126 - t175) + (pkin(5) + pkin(6)) * (-t120 * t1 + t123 * (t122 * t11 + t119 * t13)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t111 - t112, t109, t104, t110, qJDD(2), -t75, -t76, 0, 0, t170, t65, t182, -t170, -t31, t115, -t20 + t200, -t21 + t209, pkin(2) * t16, pkin(2) * t5, t170, t182, -t65, t115, t31, -t170, t131 + t200, pkin(2) * t15 + t168, t132 - t209, pkin(2) * t1 + t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t65, t182, -t170, -t31, t115, -t20, -t21, 0, 0, t170, t182, -t65, t115, t31, -t170, t131, t168, t132, t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180, t182, t181, t13;];
tauJ_reg = t2;
