% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:54
% EndTime: 2019-12-05 16:07:02
% DurationCPUTime: 1.89s
% Computational Cost: add. (2623->211), mult. (5939->258), div. (0->0), fcn. (3933->8), ass. (0->132)
t118 = sin(qJ(3));
t120 = cos(qJ(3));
t116 = cos(pkin(8));
t115 = sin(pkin(8));
t143 = qJD(2) * t118;
t83 = -t116 * t120 * qJD(2) + t115 * t143;
t144 = t118 * t116;
t85 = (t120 * t115 + t144) * qJD(2);
t164 = t85 * t83;
t171 = qJDD(3) + t164;
t159 = t115 * t171;
t121 = qJD(3) ^ 2;
t82 = t85 ^ 2;
t175 = -t82 - t121;
t23 = -t116 * t175 + t159;
t154 = t116 * t171;
t25 = t115 * t175 + t154;
t212 = pkin(6) * (t118 * t23 - t120 * t25);
t211 = pkin(3) * t23;
t210 = qJ(4) * t23;
t209 = qJ(4) * t25;
t208 = t118 * t25 + t120 * t23;
t170 = t83 ^ 2;
t68 = t170 - t121;
t206 = t118 * (-t116 * t68 + t159) - t120 * (t115 * t68 + t154);
t103 = t118 * qJDD(2);
t140 = qJD(2) * qJD(3);
t138 = t120 * t140;
t90 = t103 + t138;
t104 = t120 * qJDD(2);
t139 = t118 * t140;
t91 = t104 - t139;
t135 = t115 * t90 - t116 * t91;
t79 = qJD(3) * t85;
t39 = t135 - t79;
t148 = qJD(3) * t83;
t65 = t115 * t91 + t116 * t90;
t43 = t65 + t148;
t190 = -t115 * t39 - t116 * t43;
t205 = pkin(3) * t190;
t204 = qJ(4) * t190;
t189 = t115 * t43 - t116 * t39;
t37 = t82 + t170;
t201 = pkin(3) * t37 + qJ(4) * t189;
t200 = t118 * t189 + t120 * t190;
t199 = pkin(2) * t37 + pkin(6) * (-t118 * t190 + t120 * t189);
t172 = qJDD(3) - t164;
t158 = t115 * t172;
t173 = -t170 - t121;
t178 = t116 * t173 - t158;
t194 = qJ(4) * t178;
t47 = t116 * t172;
t180 = t115 * t173 + t47;
t193 = qJ(4) * t180;
t176 = t65 - t148;
t192 = qJ(5) * t176;
t146 = qJD(4) * t85;
t191 = pkin(3) * t180 - 0.2e1 * t146;
t38 = t135 + t79;
t188 = t118 * (t115 * t176 + t116 * t38) - t120 * (-t115 * t38 + t116 * t176);
t187 = t118 * t178 + t120 * t180;
t69 = -t82 + t121;
t186 = t118 * (-t115 * t69 + t47) + t120 * (t116 * t69 + t158);
t185 = pkin(6) * (-t118 * t180 + t120 * t178) - pkin(2) * t38;
t174 = t82 - t170;
t112 = -g(3) + qJDD(1);
t122 = qJD(2) ^ 2;
t119 = sin(qJ(2));
t150 = sin(pkin(7));
t151 = cos(pkin(7));
t129 = t150 * g(1) - t151 * g(2);
t168 = cos(qJ(2));
t93 = -t151 * g(1) - t150 * g(2);
t128 = -t119 * t129 - t168 * t93;
t145 = qJDD(2) * pkin(6);
t59 = -t122 * pkin(2) - t128 + t145;
t45 = -t120 * t112 + t118 * t59;
t100 = t118 * t122 * t120;
t95 = qJDD(3) + t100;
t123 = (-t90 + t138) * qJ(4) + t95 * pkin(3) - t45;
t111 = t120 ^ 2;
t106 = t111 * t122;
t46 = t118 * t112 + t120 * t59;
t94 = qJD(3) * pkin(3) - qJ(4) * t143;
t27 = -pkin(3) * t106 + t91 * qJ(4) - qJD(3) * t94 + t46;
t14 = -0.2e1 * qJD(4) * t83 + t115 * t123 + t116 * t27;
t169 = 2 * qJD(5);
t167 = pkin(4) * t115;
t166 = pkin(4) * t116;
t136 = t115 * t27 - t116 * t123;
t13 = t136 + 0.2e1 * t146;
t3 = t115 * t14 - t116 * t13;
t165 = t118 * t3;
t133 = -t119 * t93 + t168 * t129;
t58 = -qJDD(2) * pkin(2) - t122 * pkin(6) - t133;
t28 = -t91 * pkin(3) - qJ(4) * t106 + t94 * t143 + qJDD(4) + t58;
t162 = t115 * t28;
t156 = t116 * t28;
t153 = t118 * t95;
t96 = qJDD(3) - t100;
t152 = t120 * t96;
t149 = qJ(5) * t116;
t142 = qJD(3) * t115;
t141 = qJD(3) * t116;
t137 = -qJ(5) * t115 - pkin(3);
t4 = t115 * t13 + t116 * t14;
t134 = t118 * t45 + t120 * t46;
t48 = t83 * pkin(4) - t85 * qJ(5);
t131 = qJDD(3) * qJ(5) + qJD(3) * t169 - t83 * t48 + t14;
t130 = -qJDD(3) * pkin(4) - t121 * qJ(5) + qJDD(5) + t136;
t10 = (0.2e1 * qJD(4) + t48) * t85 + t130;
t127 = t135 * pkin(4) - t192 + t28;
t126 = t118 * (t115 * t135 + t83 * t141) + t120 * (-t116 * t135 + t83 * t142);
t67 = t85 * t142;
t125 = t118 * t67 + (-t83 * t144 + t120 * (-t115 * t83 - t116 * t85)) * qJD(3);
t124 = t85 * t169 - t127;
t110 = t118 ^ 2;
t105 = t110 * t122;
t99 = -t106 - t121;
t98 = -t105 - t121;
t92 = t104 - 0.2e1 * t139;
t89 = t103 + 0.2e1 * t138;
t15 = t118 * (t116 * t65 - t67) + t120 * (t115 * t65 + t85 * t141);
t11 = (pkin(4) * qJD(3) - (2 * qJD(5))) * t85 + t127;
t9 = -t121 * pkin(4) + t131;
t8 = (-t38 - t79) * pkin(4) + t124;
t7 = -pkin(4) * t79 + t124 + t192;
t6 = qJ(5) * t37 + t10;
t5 = (-t121 + t37) * pkin(4) + t131;
t2 = t115 * t10 + t116 * t9;
t1 = -t116 * t10 + t115 * t9;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t112, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, 0, 0, 0, 0, 0, 0, t118 * t99 + t120 * t95, -t118 * t96 + t120 * t98, 0, t118 * t46 - t120 * t45, 0, 0, 0, 0, 0, 0, t187, -t208, t200, t118 * t4 + t120 * t3, 0, 0, 0, 0, 0, 0, t187, t200, t208, t120 * t1 + t118 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t133, t128, 0, 0, (t90 + t138) * t118, t118 * t92 + t120 * t89, t153 + t120 * (-t105 + t121), (t91 - t139) * t120, t118 * (t106 - t121) + t152, 0, -t120 * t58 + pkin(2) * t92 + pkin(6) * (t120 * t99 - t153), t118 * t58 - pkin(2) * t89 + pkin(6) * (-t118 * t98 - t152), pkin(2) * (t105 + t106) + (t110 + t111) * t145 + t134, -pkin(2) * t58 + pkin(6) * t134, t15, -t188, t186, t126, -t206, t125, t118 * (t162 - t193) + t120 * (-pkin(3) * t38 - t156 + t194) + t185, t118 * (t156 + t210) + t120 * (-pkin(3) * t176 + t162 - t209) - pkin(2) * t176 + t212, t118 * (-t3 - t204) + t120 * (t201 + t4) + t199, -qJ(4) * t165 + t120 * (-pkin(3) * t28 + qJ(4) * t4) - pkin(2) * t28 + pkin(6) * (t120 * t4 - t165), t15, t186, t188, t125, t206, t126, t118 * (-t115 * t8 - t149 * t38 - t193) + t120 * (t116 * t8 + t137 * t38 + t194) + t185, t118 * (-t115 * t5 + t116 * t6 - t204) + t120 * (t115 * t6 + t116 * t5 + t201) + t199, t118 * (t116 * t7 - t210) + t120 * (t115 * t7 + t209) - t212 + (-t118 * t167 + t120 * (pkin(3) + t166) + pkin(2)) * t176, (t118 * (-t149 + t167) + t120 * (t137 - t166) - pkin(2)) * t11 + (pkin(6) + qJ(4)) * (-t118 * t1 + t120 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t105 - t106, t103, t100, t104, qJDD(3), -t45, -t46, 0, 0, t164, t174, t43, -t164, -t39, qJDD(3), -t136 + t191, -t14 - t211, t205, pkin(3) * t3, t164, t43, -t174, qJDD(3), t39, -t164, pkin(4) * t172 + qJ(5) * t173 - t85 * t48 - t130 + t191, -pkin(4) * t43 - qJ(5) * t39 + t205, t211 + qJ(5) * t171 + (-t121 - t175) * pkin(4) + t131, pkin(3) * t1 - pkin(4) * t10 + qJ(5) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t176, -t37, t28, 0, 0, 0, 0, 0, 0, t38, -t37, -t176, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, t43, t175, t10;];
tauJ_reg = t12;
