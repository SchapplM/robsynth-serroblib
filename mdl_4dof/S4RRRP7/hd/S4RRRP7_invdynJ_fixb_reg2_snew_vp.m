% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRP7
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:04
% EndTime: 2019-12-31 17:21:09
% DurationCPUTime: 1.48s
% Computational Cost: add. (2367->210), mult. (4785->243), div. (0->0), fcn. (2992->6), ass. (0->144)
t119 = cos(qJ(3));
t116 = sin(qJ(3));
t117 = sin(qJ(2));
t157 = qJD(1) * t117;
t89 = -t119 * qJD(2) + t116 * t157;
t91 = t116 * qJD(2) + t119 * t157;
t179 = t91 * t89;
t155 = qJD(1) * qJD(2);
t105 = t117 * t155;
t120 = cos(qJ(2));
t154 = t120 * qJDD(1);
t95 = -t105 + t154;
t88 = -qJDD(3) + t95;
t127 = t88 - t179;
t174 = t116 * t127;
t103 = t120 * qJD(1) - qJD(3);
t183 = t103 ^ 2;
t87 = t91 ^ 2;
t190 = -t87 - t183;
t19 = -t119 * t190 - t174;
t221 = pkin(1) * t19;
t220 = pkin(2) * t19;
t219 = pkin(6) * t19;
t167 = t119 * t127;
t25 = -t116 * t190 + t167;
t218 = pkin(6) * t25;
t217 = t120 * t25;
t163 = t89 * t103;
t106 = t117 * qJDD(1);
t151 = t120 * t155;
t94 = t106 + t151;
t133 = -t116 * qJDD(2) - t119 * t94;
t56 = -t89 * qJD(3) - t133;
t191 = t56 + t163;
t216 = qJ(4) * t191;
t145 = t119 * qJDD(2) - t116 * t94;
t130 = t91 * qJD(3) - t145;
t77 = t91 * t103;
t34 = t130 + t77;
t184 = t89 ^ 2;
t70 = t184 - t183;
t215 = t117 * (t119 * t70 + t174) + t120 * t34;
t177 = t116 * t191;
t189 = t87 - t184;
t193 = t130 - t77;
t214 = t117 * (t119 * t193 + t177) + t120 * t189;
t48 = t88 + t179;
t166 = t119 * t48;
t187 = -t183 - t184;
t196 = t116 * t187 - t166;
t213 = pkin(2) * t196;
t173 = t116 * t48;
t195 = t119 * t187 + t173;
t212 = pkin(6) * t195;
t211 = pkin(6) * t196;
t71 = -t87 + t183;
t210 = t119 * t71 - t173;
t208 = t116 * t70 - t167;
t192 = t56 - t163;
t206 = t117 * (-t116 * t71 - t166) - t120 * t192;
t205 = pkin(5) * (t117 * t193 + t120 * t195) - pkin(1) * t196;
t188 = t87 + t184;
t204 = pkin(2) * t188;
t122 = qJD(2) ^ 2;
t123 = qJD(1) ^ 2;
t139 = -t120 * pkin(2) - t117 * pkin(6);
t118 = sin(qJ(1));
t121 = cos(qJ(1));
t137 = t121 * g(1) + t118 * g(2);
t160 = qJDD(1) * pkin(5);
t81 = -t123 * pkin(1) - t137 + t160;
t146 = t123 * t139 + t81;
t180 = t120 * g(3);
t42 = -qJDD(2) * pkin(2) - t122 * pkin(6) + t146 * t117 + t180;
t203 = pkin(3) * t130 - t216 + t42;
t201 = t117 * t188;
t194 = -t116 * t193 + t119 * t191;
t135 = t94 + t151;
t136 = -t95 + t105;
t149 = t118 * g(1) - t121 * g(2);
t80 = qJDD(1) * pkin(1) + t123 * pkin(5) + t149;
t31 = t136 * pkin(2) - t135 * pkin(6) - t80;
t181 = t117 * g(3);
t43 = -t122 * pkin(2) + qJDD(2) * pkin(6) + t146 * t120 - t181;
t18 = t116 * t31 + t119 * t43;
t59 = t89 * pkin(3) - t91 * qJ(4);
t147 = t88 * qJ(4) + t89 * t59 - t18;
t185 = -pkin(3) * (t183 + t190) - qJ(4) * t127 - t147;
t182 = pkin(3) * t119;
t176 = t116 * t192;
t175 = t116 * t42;
t102 = t120 * t123 * t117;
t171 = t117 * (qJDD(2) + t102);
t169 = t119 * t192;
t168 = t119 * t42;
t164 = t120 * (qJDD(2) - t102);
t161 = qJ(4) * t119;
t159 = t103 * t116;
t158 = t103 * t119;
t156 = qJD(4) * t103;
t153 = t120 * t179;
t152 = t89 * t158;
t150 = -qJ(4) * t116 - pkin(2);
t17 = t116 * t43 - t119 * t31;
t7 = t116 * t17 + t119 * t18;
t67 = t117 * t81 + t180;
t68 = t120 * t81 - t181;
t148 = t117 * t67 + t120 * t68;
t69 = t91 * t159;
t144 = t117 * (t119 * t56 + t69) - t153;
t99 = -0.2e1 * t156;
t143 = t99 - t147;
t142 = -t119 * t130 - t89 * t159;
t11 = t88 * pkin(3) - qJ(4) * t183 + t91 * t59 + qJDD(4) + t17;
t9 = -pkin(3) * t183 + t143;
t140 = -pkin(3) * t11 + qJ(4) * t9;
t138 = -pkin(3) * t192 - qJ(4) * t34;
t134 = t116 * t18 - t119 * t17;
t132 = -pkin(1) + t139;
t129 = (t116 * t89 + t119 * t91) * t103;
t128 = t117 * (-t69 + t152) + t120 * t88;
t126 = t117 * (t116 * t130 - t152) + t153;
t125 = 0.2e1 * qJD(4) * t91 - t203;
t124 = -pkin(3) * t48 + qJ(4) * t187 - t11;
t113 = t120 ^ 2;
t112 = t117 ^ 2;
t110 = t113 * t123;
t108 = t112 * t123;
t96 = -0.2e1 * t105 + t154;
t93 = t106 + 0.2e1 * t151;
t40 = (qJD(3) - t103) * t89 + t133;
t35 = (-qJD(3) - t103) * t91 + t145;
t28 = t116 * t56 - t91 * t158;
t15 = t119 * t35 + t176;
t14 = -t119 * t34 + t176;
t12 = -t116 * t34 - t169;
t10 = (-pkin(3) * t103 - 0.2e1 * qJD(4)) * t91 + t203;
t8 = qJ(4) * t188 + t11;
t5 = (-t183 + t188) * pkin(3) + t143;
t4 = (-t193 + t77) * pkin(3) + t125;
t3 = pkin(3) * t77 + t125 + t216;
t2 = t116 * t11 + t119 * t9;
t1 = -t119 * t11 + t116 * t9;
t6 = [0, 0, 0, 0, 0, qJDD(1), t149, t137, 0, 0, t135 * t117, t117 * t96 + t120 * t93, t171 + t120 * (-t108 + t122), -t136 * t120, t117 * (t110 - t122) + t164, 0, t120 * t80 + pkin(1) * t96 + pkin(5) * (t120 * (-t110 - t122) - t171), -t117 * t80 - pkin(1) * t93 + pkin(5) * (-t164 - t117 * (-t108 - t122)), pkin(1) * (t108 + t110) + (t112 + t113) * t160 + t148, pkin(1) * t80 + pkin(5) * t148, t144, -t214, t206, t126, t215, t128, t117 * (t175 - t211) + t120 * (t17 - t213) + t205, t117 * (t168 + t219) + t120 * (t18 + t220) + t221 + pkin(5) * (-t117 * t40 + t217), -t117 * t134 + pkin(5) * (t120 * t15 - t201) + t132 * (t116 * t35 - t169), pkin(5) * (t117 * t42 + t120 * t7) + t132 * t134, t144, t206, t214, t128, -t215, t126, t117 * (-t116 * t4 - t161 * t193 - t211) + t120 * (-t124 - t213) + t205, t117 * (-pkin(6) * t12 - t116 * t5 + t119 * t8) + t120 * (-pkin(2) * t12 - t138) - pkin(1) * t12 + pkin(5) * (t120 * t14 - t201), t117 * (-pkin(3) * t177 + t119 * t3 - t219) + t120 * (0.2e1 * t156 - t185 - t220) - t221 + pkin(5) * (-t117 * t191 - t217), t117 * (-pkin(6) * t1 + (pkin(3) * t116 - t161) * t10) + t120 * (-pkin(2) * t1 - t140) - pkin(1) * t1 + pkin(5) * (t117 * t10 + t120 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, t108 - t110, t106, t102, t154, qJDD(2), -t67, -t68, 0, 0, t28, t194, t210, t142, t208, t129, -pkin(2) * t193 - t168 + t212, pkin(2) * t40 + t175 + t218, pkin(6) * t15 + t204 + t7, -pkin(2) * t42 + pkin(6) * t7, t28, t210, -t194, t129, -t208, t142, t119 * t4 + t150 * t193 + t212, pkin(6) * t14 + t116 * t8 + t119 * t5 + t204, -t218 + t116 * t3 + (pkin(2) + t182) * t191, pkin(6) * t2 + (t150 - t182) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t189, t192, -t179, -t34, -t88, -t17, -t18, 0, 0, t179, t192, -t189, -t88, t34, -t179, t124, t138, t185 + t99, t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t192, t190, t11;];
tauJ_reg = t6;
