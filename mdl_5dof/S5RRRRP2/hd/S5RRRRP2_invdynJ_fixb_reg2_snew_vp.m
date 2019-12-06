% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:48:13
% EndTime: 2019-12-05 18:48:21
% DurationCPUTime: 1.83s
% Computational Cost: add. (9373->226), mult. (12062->288), div. (0->0), fcn. (7707->8), ass. (0->157)
t168 = qJD(1) + qJD(2);
t171 = sin(qJ(4));
t175 = cos(qJ(4));
t176 = cos(qJ(3));
t172 = sin(qJ(3));
t212 = t168 * t172;
t133 = -t168 * t175 * t176 + t171 * t212;
t135 = (t171 * t176 + t172 * t175) * t168;
t105 = t135 * t133;
t165 = qJDD(3) + qJDD(4);
t233 = t105 - t165;
t236 = t233 * pkin(4);
t235 = t171 * t233;
t234 = t175 * t233;
t174 = sin(qJ(1));
t178 = cos(qJ(1));
t186 = -g(2) * t174 + g(3) * t178;
t149 = -qJD(1) ^ 2 * pkin(1) - t186;
t173 = sin(qJ(2));
t177 = cos(qJ(2));
t187 = g(2) * t178 + g(3) * t174;
t182 = qJDD(1) * pkin(1) + t187;
t108 = t177 * t149 + t173 * t182;
t164 = t168 ^ 2;
t166 = qJDD(1) + qJDD(2);
t102 = -pkin(2) * t164 + pkin(7) * t166 + t108;
t211 = t172 * t102;
t91 = t176 * g(1) + t211;
t92 = -g(1) * t172 + t102 * t176;
t58 = t172 * t91 + t176 * t92;
t167 = qJD(3) + qJD(4);
t125 = t167 * t133;
t207 = qJD(3) * t168;
t197 = t176 * t207;
t209 = t172 * t166;
t142 = t197 + t209;
t157 = t176 * t166;
t198 = t172 * t207;
t143 = t157 - t198;
t87 = -qJD(4) * t133 + t142 * t175 + t143 * t171;
t232 = -t125 + t87;
t215 = t164 * t172;
t66 = qJDD(3) * pkin(3) - t142 * pkin(8) - t211 + (pkin(3) * t215 + pkin(8) * t207 - g(1)) * t176;
t152 = qJD(3) * pkin(3) - pkin(8) * t212;
t170 = t176 ^ 2;
t159 = t170 * t164;
t67 = -pkin(3) * t159 + pkin(8) * t143 - qJD(3) * t152 + t92;
t44 = t171 * t67 - t175 * t66;
t230 = qJ(5) * t125 + 0.2e1 * qJD(5) * t135 + t236 + t44;
t120 = pkin(4) * t167 - qJ(5) * t135;
t131 = t133 ^ 2;
t45 = t171 * t66 + t175 * t67;
t190 = t171 * t142 - t143 * t175;
t86 = -qJD(4) * t135 - t190;
t30 = -pkin(4) * t131 + qJ(5) * t86 - 0.2e1 * qJD(5) * t133 - t120 * t167 + t45;
t132 = t135 ^ 2;
t163 = t167 ^ 2;
t181 = (-qJD(4) + t167) * t135 - t190;
t78 = t125 + t87;
t48 = t171 * t181 - t175 * t78;
t229 = pkin(8) * t48;
t96 = -t163 - t131;
t61 = t171 * t96 - t234;
t228 = pkin(8) * t61;
t119 = -t132 - t163;
t98 = t105 + t165;
t222 = t171 * t98;
t81 = t119 * t175 - t222;
t227 = pkin(8) * t81;
t49 = t171 * t78 + t175 * t181;
t26 = -t172 * t48 + t176 * t49;
t88 = -t131 - t132;
t226 = -pkin(2) * t88 + pkin(7) * t26;
t62 = t175 * t96 + t235;
t39 = -t172 * t61 + t176 * t62;
t73 = (qJD(4) + t167) * t135 + t190;
t225 = -pkin(2) * t73 + pkin(7) * t39;
t219 = t175 * t98;
t82 = -t119 * t171 - t219;
t52 = -t172 * t81 + t176 * t82;
t224 = -pkin(2) * t232 + pkin(7) * t52;
t107 = -t173 * t149 + t177 * t182;
t101 = -pkin(2) * t166 - pkin(7) * t164 - t107;
t72 = -pkin(3) * t143 - pkin(8) * t159 + t152 * t212 + t101;
t223 = t171 * t72;
t22 = t171 * t45 - t175 * t44;
t221 = t172 * t22;
t220 = t175 * t72;
t218 = -pkin(2) * t101 + pkin(7) * t58;
t217 = qJ(5) * t171;
t216 = qJ(5) * t175;
t214 = t167 * t171;
t213 = t167 * t175;
t154 = t176 * t215;
t210 = t172 * (qJDD(3) + t154);
t208 = t176 * (qJDD(3) - t154);
t169 = t172 ^ 2;
t158 = t169 * t164;
t179 = qJD(3) ^ 2;
t118 = -t208 - t172 * (-t158 - t179);
t141 = 0.2e1 * t197 + t209;
t204 = -pkin(2) * t141 + pkin(7) * t118 + t101 * t172;
t117 = t176 * (-t159 - t179) - t210;
t144 = t157 - 0.2e1 * t198;
t203 = pkin(2) * t144 + pkin(7) * t117 - t101 * t176;
t29 = -qJ(5) * t87 - t230;
t10 = -t171 * t29 + t175 * t30;
t37 = -pkin(4) * t86 - qJ(5) * t131 + t120 * t135 + qJDD(5) + t72;
t15 = -pkin(4) * t37 + qJ(5) * t30;
t9 = t171 * t30 + t175 * t29;
t4 = t10 * t176 - t172 * t9;
t202 = t172 * (-pkin(8) * t9 - t15 * t171 - t216 * t29) + t176 * (-pkin(3) * t37 + pkin(8) * t10 + t15 * t175 - t217 * t29) - pkin(2) * t37 + pkin(7) * t4;
t201 = -pkin(3) * t88 + pkin(8) * t49;
t200 = -pkin(3) * t73 + pkin(8) * t62;
t199 = -pkin(3) * t232 + pkin(8) * t82;
t19 = -pkin(4) * t88 + qJ(5) * t181 + t30;
t21 = (t78 + t87) * qJ(5) + t230;
t196 = t172 * (-t171 * t19 + t175 * t21 - t229) + t176 * (t171 * t21 + t175 * t19 + t201) + t226;
t23 = t171 * t44 + t175 * t45;
t195 = t172 * (-t22 - t229) + t176 * (t201 + t23) + t226;
t27 = -pkin(4) * t73 + qJ(5) * t96 - t37;
t194 = t172 * (-t171 * t27 + t216 * t233 - t228) + t176 * (t175 * t27 + t217 * t233 + t200) + t225;
t35 = -qJ(5) * t119 + t37;
t55 = -pkin(4) * t232 - qJ(5) * t98;
t193 = t172 * (-t171 * t55 + t175 * t35 - t227) + t176 * (t171 * t35 + t175 * t55 + t199) + t224;
t192 = t172 * (t223 - t228) + t176 * (t200 - t220) + t225;
t191 = t172 * (t220 - t227) + t176 * (t199 + t223) + t224;
t147 = (t169 + t170) * t166;
t148 = t158 + t159;
t189 = pkin(2) * t148 + pkin(7) * t147 + t58;
t8 = t176 * t23 - t221;
t185 = pkin(7) * t8 - pkin(8) * t221 + t176 * (-pkin(3) * t72 + pkin(8) * t23) - pkin(2) * t72;
t183 = pkin(4) * t119 - t30;
t180 = t29 - t236;
t122 = -t132 + t163;
t121 = t131 - t163;
t116 = t210 + t176 * (-t158 + t179);
t115 = t172 * (t159 - t179) + t208;
t111 = (t142 + t197) * t172;
t110 = (t143 - t198) * t176;
t106 = t141 * t176 + t144 * t172;
t103 = t132 - t131;
t80 = pkin(3) * t81;
t69 = pkin(4) * t78;
t60 = pkin(3) * t61;
t56 = (t172 * (-t133 * t175 + t135 * t171) + t176 * (-t133 * t171 - t135 * t175)) * t167;
t54 = t172 * (t121 * t175 - t222) + t176 * (t121 * t171 + t219);
t53 = t172 * (-t122 * t171 - t234) + t176 * (t122 * t175 - t235);
t47 = pkin(3) * t48;
t41 = t172 * (-t135 * t214 + t175 * t87) + t176 * (t135 * t213 + t171 * t87);
t40 = t172 * (t133 * t213 - t171 * t86) + t176 * (t133 * t214 + t175 * t86);
t34 = pkin(1) * (t173 * t52 - t177 * t232);
t31 = pkin(1) * (t173 * t39 - t177 * t73);
t28 = pkin(4) * t29;
t25 = t172 * (-t171 * t232 - t175 * t73) + t176 * (-t171 * t73 + t175 * t232);
t20 = pkin(1) * (t173 * t26 - t177 * t88);
t1 = [0, 0, 0, 0, 0, qJDD(1), t187, t186, 0, 0, 0, 0, 0, 0, 0, t166, pkin(1) * (-t164 * t173 + t166 * t177) + t107, pkin(1) * (-t164 * t177 - t166 * t173) - t108, 0, pkin(1) * (t107 * t177 + t108 * t173), t111, t106, t116, t110, t115, 0, pkin(1) * (t117 * t173 + t144 * t177) + t203, pkin(1) * (t118 * t173 - t141 * t177) + t204, pkin(1) * (t147 * t173 + t148 * t177) + t189, pkin(1) * (-t101 * t177 + t173 * t58) + t218, t41, t25, t53, t40, t54, t56, t31 + t192, t34 + t191, t20 + t195, pkin(1) * (t173 * t8 - t177 * t72) + t185, t41, t25, t53, t40, t54, t56, t31 + t194, t34 + t193, t20 + t196, pkin(1) * (t173 * t4 - t177 * t37) + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t107, -t108, 0, 0, t111, t106, t116, t110, t115, 0, t203, t204, t189, t218, t41, t25, t53, t40, t54, t56, t192, t191, t195, t185, t41, t25, t53, t40, t54, t56, t194, t193, t196, t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t158 - t159, t209, t154, t157, qJDD(3), -t91, -t92, 0, 0, t105, t103, t78, -t105, t181, t165, -t44 + t60, t80 - t45, t47, pkin(3) * t22, t105, t103, t78, -t105, t181, t165, t180 + t60, t80 + t183, -t69 + t47, pkin(3) * t9 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t103, t78, -t105, t181, t165, -t44, -t45, 0, 0, t105, t103, t78, -t105, t181, t165, t180, t183, -t69, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t232, t88, t37;];
tauJ_reg = t1;
