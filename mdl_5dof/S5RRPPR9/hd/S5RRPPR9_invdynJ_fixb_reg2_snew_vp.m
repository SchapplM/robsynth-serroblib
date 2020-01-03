% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:48
% EndTime: 2019-12-31 19:41:52
% DurationCPUTime: 1.36s
% Computational Cost: add. (2981->238), mult. (6420->276), div. (0->0), fcn. (3231->6), ass. (0->147)
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t122 = cos(qJ(2));
t164 = qJD(1) * t122;
t68 = -qJD(2) * t121 + t118 * t164;
t69 = qJD(2) * t118 + t121 * t164;
t40 = t68 * t69;
t161 = qJD(1) * qJD(2);
t100 = t122 * t161;
t119 = sin(qJ(2));
t159 = t119 * qJDD(1);
t74 = t100 + t159;
t62 = qJDD(5) + t74;
t201 = -t40 + t62;
t206 = t118 * t201;
t205 = t121 * t201;
t125 = qJD(1) ^ 2;
t197 = t122 * t125;
t91 = t119 * t197;
t83 = t91 + qJDD(2);
t84 = -t91 + qJDD(2);
t174 = t122 * t84;
t73 = 0.2e1 * t100 + t159;
t124 = qJD(2) ^ 2;
t112 = t119 ^ 2;
t171 = t112 * t125;
t87 = t124 + t171;
t204 = pkin(1) * t73 + pkin(6) * (-t119 * t87 + t174);
t180 = t119 * t83;
t113 = t122 ^ 2;
t170 = t113 * t125;
t88 = t124 + t170;
t203 = pkin(6) * (t122 * t88 + t180);
t202 = pkin(6) - qJ(4);
t158 = t122 * qJDD(1);
t99 = t119 * t161;
t75 = -t99 + t158;
t34 = qJD(5) * t68 - qJDD(2) * t118 - t121 * t75;
t165 = qJD(1) * t119;
t93 = qJD(5) + t165;
t51 = t93 * t68;
t200 = t51 + t34;
t199 = pkin(2) * t87 + qJ(3) * t84;
t198 = qJDD(1) * pkin(6);
t194 = pkin(2) + pkin(3);
t148 = -t122 * t194 - pkin(1);
t169 = t119 * qJ(3);
t196 = t148 - t169;
t46 = t174 + t119 * (-t124 + t170);
t120 = sin(qJ(1));
t123 = cos(qJ(1));
t146 = g(1) * t123 + g(2) * t120;
t56 = -pkin(1) * t125 - t146 + t198;
t142 = -pkin(2) * t122 - t169;
t71 = t142 * qJD(1);
t150 = qJD(1) * t71 + t56;
t152 = qJDD(2) * pkin(2) + qJ(3) * t124 - qJDD(3);
t191 = t122 * g(3);
t28 = t150 * t119 - t152 + t191;
t60 = t68 ^ 2;
t61 = t69 ^ 2;
t90 = t93 ^ 2;
t193 = pkin(2) + pkin(7);
t166 = t112 + t113;
t72 = t166 * t198;
t192 = t75 * pkin(2);
t190 = pkin(4) + qJ(3);
t76 = -0.2e1 * t99 + t158;
t189 = pkin(1) * t76 - t203;
t80 = t166 * t125;
t188 = pkin(1) * t80 + t72;
t187 = qJ(3) * t80;
t186 = qJ(3) * t88;
t160 = qJD(3) * qJD(2);
t105 = 0.2e1 * t160;
t48 = -g(3) * t119 + t122 * t56;
t145 = -pkin(2) * t124 + qJDD(2) * qJ(3) + t164 * t71 + t48;
t155 = 0.2e1 * qJD(1) * qJD(4);
t82 = -qJD(2) * pkin(3) - qJ(4) * t165;
t127 = pkin(3) * t170 + t75 * qJ(4) - qJD(2) * t82 + t122 * t155 - t145;
t13 = t105 - t127;
t147 = pkin(4) * t119 + pkin(7) * t122;
t10 = qJDD(2) * pkin(4) - pkin(7) * t124 - t147 * t197 + t13;
t185 = t10 * t118;
t184 = t10 * t121;
t32 = t40 + t62;
t183 = t118 * t32;
t182 = t118 * t93;
t181 = t119 * t76;
t177 = t121 * t32;
t176 = t121 * t93;
t175 = t122 * t73;
t173 = qJ(3) * t122;
t172 = qJ(4) * t119;
t168 = -0.2e1 * qJD(4) + t71;
t167 = g(1) * t120 - g(2) * t123;
t163 = qJD(2) * t122;
t162 = pkin(3) + t193;
t157 = -t61 - t90;
t156 = t119 * t40;
t154 = qJ(4) * t163;
t137 = pkin(3) * t83 + t74 * qJ(4) + t152;
t47 = t119 * t56 + t191;
t131 = -t137 + t47;
t11 = -t124 * pkin(4) - qJDD(2) * pkin(7) + (t154 + (-qJD(1) * t147 + t168) * t119) * qJD(1) + t131;
t55 = qJDD(1) * pkin(1) + t125 * pkin(6) + t167;
t135 = -pkin(2) * t99 + t55;
t130 = pkin(3) * t75 - qJ(4) * t170 + qJDD(4) + t135;
t149 = (0.2e1 * qJD(3) + t82) * t119;
t7 = t193 * t75 + t190 * t74 + (t149 + (-pkin(7) * t119 + t122 * t190) * qJD(2)) * qJD(1) + t130;
t4 = t11 * t118 - t121 * t7;
t153 = qJD(3) * t165;
t151 = t119 * t47 + t122 * t48;
t5 = t11 * t121 + t118 * t7;
t1 = t118 * t5 - t121 * t4;
t2 = t118 * t4 + t121 * t5;
t144 = t74 + t100;
t37 = t175 + t181;
t139 = t121 * qJDD(2) - t118 * t75;
t27 = t105 + t145;
t134 = t119 * t190 + t122 * t162 + pkin(1);
t133 = (qJD(5) - t93) * t69 - t139;
t132 = (t144 + t73) * t169 + t204;
t129 = t130 + t192;
t128 = t135 + 0.2e1 * t153 + t192;
t126 = t74 * qJ(3) + t129;
t14 = (t119 * t168 + t154) * qJD(1) + t131;
t81 = (t112 - t113) * t125;
t50 = -t61 + t90;
t49 = t60 - t90;
t45 = t144 * t119;
t44 = t180 + t122 * (t124 - t171);
t43 = (t75 - t99) * t122;
t38 = t61 - t60;
t35 = -t90 - t60;
t33 = qJD(5) * t69 - t139;
t30 = -t60 - t61;
t26 = -t51 + t34;
t21 = (-qJD(5) - t93) * t69 + t139;
t18 = -t118 * t157 - t177;
t17 = t121 * t157 - t183;
t16 = t121 * t35 - t206;
t15 = t118 * t35 + t205;
t12 = (qJ(3) * t163 + t149) * qJD(1) + t126;
t9 = t118 * t26 + t121 * t133;
t8 = t118 * t133 - t121 * t26;
t3 = [0, 0, 0, 0, 0, qJDD(1), t167, t146, 0, 0, t45, t37, t44, t43, t46, 0, t122 * t55 + t189, -t119 * t55 - t204, t151 + t188, pkin(1) * t55 + pkin(6) * t151, t45, t44, -t37, 0, -t46, t43, t122 * (pkin(2) * t76 + t128) + (t122 * t144 + t181) * qJ(3) + t189, t122 * (pkin(2) * t80 + t27) + (t28 + t187) * t119 + t188, pkin(2) * t175 + t119 * t128 + t132, pkin(6) * (t119 * t28 + t122 * t27) + (pkin(1) - t142) * (qJ(3) * t144 + t128), t43, t37, t46, t45, t44, 0, t122 * (-qJ(4) * t84 + t194 * t73) + (qJ(4) * t87 + qJD(1) * t149 + t129) * t119 + t132, -t83 * t172 + t122 * (-qJ(3) * t100 - qJ(4) * t88 - t165 * t82 - t126 - 0.2e1 * t153) + t203 + t196 * t76, t122 * (qJ(4) * t158 + t127 - 0.2e1 * t160) - t72 + t148 * t80 + (-t187 + (-qJ(4) * t161 - g(3)) * t122 + t137 + (qJ(4) * qJDD(1) - t150 + t155) * t119) * t119, -t196 * t12 + t202 * (t119 * t14 + t122 * t13), t156 + t122 * (-t121 * t34 - t182 * t69), t119 * t38 + t122 * (t118 * t200 + t121 * t21), t119 * t26 + t122 * (t118 * t50 - t205), -t156 + t122 * (t118 * t33 + t176 * t68), t119 * t133 + t122 * (-t121 * t49 + t183), t119 * t62 + t122 * (t118 * t69 - t121 * t68) * t93, t119 * (-qJ(4) * t16 - t4) + t122 * (-qJ(4) * t21 - t185) + pkin(6) * (t119 * t16 + t122 * t21) + t134 * t15, t119 * (-qJ(4) * t18 - t5) + t122 * (-qJ(4) * t200 - t184) + pkin(6) * (t119 * t18 + t122 * t200) + t134 * t17, -t9 * t172 + t122 * (-qJ(4) * t30 + t1) + pkin(6) * (t119 * t9 + t122 * t30) + t134 * t8, t1 * t134 + t202 * (t10 * t122 + t119 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t81, t159, t91, t158, qJDD(2), -t47, -t48, 0, 0, -t91, t159, -t81, qJDD(2), -t158, t91, pkin(2) * t83 - t186 - t28, (-pkin(2) * t119 + t173) * qJDD(1), t27 + t199, -pkin(2) * t28 + qJ(3) * t27, t91, t81, t158, -t91, t159, qJDD(2), pkin(3) * t87 + t13 + t199, -t194 * t83 + t14 + t186, (t119 * t194 - t173) * qJDD(1), qJ(3) * t13 - t14 * t194, -t118 * t34 + t176 * t69, t118 * t21 - t121 * t200, -t121 * t50 - t206, -t121 * t33 + t182 * t68, -t118 * t49 - t177, (-t118 * t68 - t121 * t69) * t93, -t16 * t162 + t190 * t21 + t184, -t162 * t18 + t190 * t200 - t185, -t162 * t9 + t190 * t30 - t2, t10 * t190 - t162 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t159, -t87, t28, 0, 0, 0, 0, 0, 0, -t87, t83, -t159, t14, 0, 0, 0, 0, 0, 0, t16, t18, t9, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t76, -t80, t12, 0, 0, 0, 0, 0, 0, t15, t17, t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t38, t26, -t40, t133, t62, -t4, -t5, 0, 0;];
tauJ_reg = t3;
