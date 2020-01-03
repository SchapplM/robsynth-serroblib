% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:38
% EndTime: 2019-12-31 18:24:42
% DurationCPUTime: 1.28s
% Computational Cost: add. (2720->211), mult. (5498->263), div. (0->0), fcn. (3041->8), ass. (0->127)
t116 = sin(pkin(8));
t117 = cos(pkin(8));
t124 = cos(qJ(3));
t121 = sin(qJ(3));
t127 = qJD(1) ^ 2;
t161 = t124 * t127;
t151 = t121 * t161;
t87 = qJDD(3) + t151;
t170 = t121 * t87;
t112 = t124 ^ 2;
t107 = t112 * t127;
t126 = qJD(3) ^ 2;
t93 = -t107 - t126;
t52 = -t124 * t93 + t170;
t158 = qJD(1) * qJD(3);
t102 = t121 * t158;
t155 = t124 * qJDD(1);
t82 = -0.2e1 * t102 + t155;
t191 = pkin(1) * (t116 * t52 - t117 * t82) - pkin(2) * t82 + pkin(6) * t52;
t103 = t124 * t158;
t156 = t121 * qJDD(1);
t80 = t103 + t156;
t190 = t80 + t103;
t120 = sin(qJ(5));
t123 = cos(qJ(5));
t160 = qJD(1) * t124;
t72 = t120 * qJD(3) + t123 * t160;
t74 = t123 * qJD(3) - t120 * t160;
t48 = t74 * t72;
t67 = qJDD(5) + t80;
t185 = -t48 + t67;
t189 = t120 * t185;
t188 = t123 * t185;
t81 = -t102 + t155;
t42 = -t72 * qJD(5) + t123 * qJDD(3) - t120 * t81;
t159 = t121 * qJD(1);
t97 = qJD(5) + t159;
t57 = t97 * t72;
t184 = -t57 + t42;
t163 = t121 * qJ(4);
t177 = t124 * pkin(3);
t122 = sin(qJ(1));
t125 = cos(qJ(1));
t149 = t122 * g(1) - t125 * g(2);
t75 = qJDD(1) * pkin(1) + t149;
t143 = t125 * g(1) + t122 * g(2);
t77 = -t127 * pkin(1) - t143;
t176 = t116 * t75 + t117 * t77;
t40 = -t127 * pkin(2) + qJDD(1) * pkin(6) + t176;
t147 = t127 * (-t163 - t177) + t40;
t183 = -t126 * pkin(3) + t147 * t124;
t88 = qJDD(3) - t151;
t164 = t124 * t88;
t182 = t164 + t121 * (t107 - t126);
t111 = t121 ^ 2;
t106 = t111 * t127;
t91 = -t106 - t126;
t181 = t121 * t91 + t164;
t65 = t72 ^ 2;
t66 = t74 ^ 2;
t94 = t97 ^ 2;
t180 = 2 * qJD(4);
t179 = -pkin(3) - pkin(7);
t157 = qJDD(3) * qJ(4);
t113 = -g(3) + qJDD(2);
t162 = t121 * t113;
t89 = pkin(4) * t159 - qJD(3) * pkin(7);
t15 = t157 + t162 - pkin(7) * t107 + t81 * pkin(4) + (t180 + t89) * qJD(3) + t183;
t174 = t120 * t15;
t38 = t48 + t67;
t173 = t120 * t38;
t172 = t120 * t97;
t145 = t120 * qJDD(3) + t123 * t81;
t133 = (-qJD(5) + t97) * t74 - t145;
t31 = t57 + t42;
t11 = t120 * t133 - t123 * t31;
t171 = t121 * t11;
t167 = t123 * t15;
t166 = t123 * t38;
t165 = t123 * t97;
t101 = t124 * t113;
t154 = -t66 - t94;
t84 = (t111 + t112) * qJDD(1);
t85 = t106 + t107;
t153 = pkin(2) * t85 + pkin(6) * t84 + pkin(1) * (t116 * t84 + t117 * t85);
t152 = t121 * t48;
t150 = pkin(1) * t116 + pkin(6);
t148 = -t116 * t77 + t117 * t75;
t39 = -qJDD(1) * pkin(2) - t127 * pkin(6) - t148;
t131 = -t81 * pkin(3) - t190 * qJ(4) + t39;
t146 = pkin(3) * qJD(3) - (2 * qJD(4));
t13 = -pkin(4) * t107 - t81 * pkin(7) + (t146 - t89) * t159 + t131;
t144 = -qJDD(3) * pkin(3) - t126 * qJ(4) + qJDD(4) - t101;
t16 = -qJDD(3) * pkin(7) + (t80 - t103) * pkin(4) + (-pkin(7) * t161 + t147) * t121 + t144;
t5 = t120 * t13 - t123 * t16;
t33 = t121 * t40 - t101;
t34 = t124 * t40 + t162;
t17 = t121 * t33 + t124 * t34;
t6 = t120 * t16 + t123 * t13;
t2 = t120 * t6 - t123 * t5;
t142 = t120 * t5 + t123 * t6;
t140 = t124 * (-t106 + t126) + t170;
t139 = t121 * t93 + t124 * t87;
t138 = t121 * t88 - t124 * t91;
t136 = qJD(1) * t146;
t135 = -pkin(1) * t117 - pkin(2) - t163;
t132 = -t135 + t177;
t130 = t124 * t179 + t135;
t129 = qJD(3) * t180 + t183;
t24 = t147 * t121 + t144;
t128 = t129 + t157;
t86 = t106 - t107;
t79 = 0.2e1 * t103 + t156;
t56 = -t66 + t94;
t55 = t65 - t94;
t50 = t190 * t121;
t49 = (t81 - t102) * t124;
t47 = t66 - t65;
t45 = t121 * t82 + t124 * t79;
t43 = -t94 - t65;
t41 = -t74 * qJD(5) - t145;
t36 = -t65 - t66;
t26 = (qJD(5) + t97) * t74 + t145;
t22 = t123 * t154 - t173;
t21 = t128 + t162;
t19 = t120 * t43 + t188;
t1 = [0, 0, 0, 0, 0, qJDD(1), t149, t143, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t117 * qJDD(1) - t116 * t127) + t148, pkin(1) * (-t116 * qJDD(1) - t117 * t127) - t176, 0, pkin(1) * (t116 * t176 + t117 * t148), t50, t45, t140, t49, t182, 0, -t124 * t39 - t191, t121 * t39 - pkin(2) * t79 - pkin(6) * t181 + pkin(1) * (-t116 * t181 - t117 * t79), t17 + t153, -pkin(2) * t39 + pkin(6) * t17 + pkin(1) * (t116 * t17 - t117 * t39), 0, -t140, -t182, t50, t45, t49, (pkin(3) * t85 + t128) * t124 + (qJ(4) * t85 + t101 + t24) * t121 + t153, t124 * (-pkin(3) * t82 + t131) + (-qJ(4) * t82 + t124 * t136) * t121 + t191, t121 * (-pkin(3) * t102 + t159 * t180 - t131) + t150 * t181 + t132 * t79, t150 * (t121 * t24 + t124 * t21) - t132 * (t121 * t136 + t131), t152 + t124 * (-t120 * t42 - t74 * t165), t121 * t47 + t124 * (t120 * t26 - t123 * t184), t121 * t31 + t124 * (-t123 * t56 - t189), -t152 + t124 * (-t123 * t41 - t72 * t172), t121 * t133 + t124 * (-t120 * t55 - t166), t121 * t67 + t124 * (t120 * t72 + t123 * t74) * t97, t121 * (pkin(4) * t19 - t5) + t124 * (pkin(4) * t26 + t167) + t150 * (t121 * t19 + t124 * t26) + t130 * (t123 * t43 - t189), t121 * (pkin(4) * t22 - t6) + t124 * (pkin(4) * t184 - t174) + t150 * (t121 * t22 + t124 * t184) + t130 * (-t120 * t154 - t166), pkin(4) * t171 + t124 * (pkin(4) * t36 - t142) + t150 * (t124 * t36 + t171) + t130 * (t120 * t31 + t123 * t133), t130 * t142 + (pkin(4) + t150) * (t121 * t2 + t124 * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, 0, 0, 0, 0, 0, t139, -t138, 0, t121 * t34 - t124 * t33, 0, 0, 0, 0, 0, 0, 0, -t139, t138, t121 * t21 - t124 * t24, 0, 0, 0, 0, 0, 0, t121 * t26 - t124 * t19, t121 * t184 - t124 * t22, -t124 * t11 + t121 * t36, t121 * t15 - t124 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t86, t156, t151, t155, qJDD(3), -t33, -t34, 0, 0, qJDD(3), -t156, -t155, -t151, t86, t151, (-pkin(3) * t121 + qJ(4) * t124) * qJDD(1), -pkin(3) * t87 - qJ(4) * t93 + t24, -pkin(3) * t91 + t162 + (qJDD(3) + t88) * qJ(4) + t129, -pkin(3) * t24 + qJ(4) * t21, t123 * t42 - t74 * t172, -t120 * t184 - t123 * t26, -t120 * t56 + t188, -t120 * t41 + t72 * t165, t123 * t55 - t173, (t120 * t74 - t123 * t72) * t97, qJ(4) * t26 + t179 * t19 + t174, qJ(4) * t184 + t179 * t22 + t167, qJ(4) * t36 + t179 * t11 - t2, qJ(4) * t15 + t179 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t87, t91, t24, 0, 0, 0, 0, 0, 0, t19, t22, t11, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t47, t31, -t48, t133, t67, -t5, -t6, 0, 0;];
tauJ_reg = t1;
