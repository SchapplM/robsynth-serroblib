% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP13_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:37
% EndTime: 2019-12-31 18:59:43
% DurationCPUTime: 1.78s
% Computational Cost: add. (3184->219), mult. (6182->252), div. (0->0), fcn. (3698->6), ass. (0->165)
t128 = cos(qJ(4));
t125 = sin(qJ(4));
t129 = cos(qJ(3));
t170 = qJD(1) * qJD(3);
t163 = t129 * t170;
t126 = sin(qJ(3));
t168 = t126 * qJDD(1);
t107 = -t163 - t168;
t102 = qJDD(4) - t107;
t172 = qJD(1) * t129;
t103 = -t128 * qJD(3) + t125 * t172;
t105 = t125 * qJD(3) + t128 * t172;
t181 = t105 * t103;
t204 = t102 + t181;
t192 = t125 * t204;
t101 = t105 ^ 2;
t116 = t126 * qJD(1) + qJD(4);
t200 = t116 ^ 2;
t208 = -t101 - t200;
t27 = -t128 * t208 + t192;
t244 = pkin(3) * t27;
t243 = pkin(7) * t27;
t186 = t128 * t204;
t29 = t125 * t208 + t186;
t242 = pkin(7) * t29;
t241 = qJ(2) * t27;
t240 = t126 * t29;
t118 = t129 * qJDD(1);
t164 = t126 * t170;
t108 = t118 - t164;
t158 = t128 * qJDD(3) - t125 * t108;
t141 = t105 * qJD(4) - t158;
t92 = t116 * t105;
t44 = t141 - t92;
t201 = t103 ^ 2;
t84 = t201 - t200;
t239 = t126 * t44 + t129 * (-t128 * t84 + t192);
t238 = t125 * t84 + t186;
t182 = t103 * t116;
t146 = -t125 * qJDD(3) - t128 * t108;
t65 = -t103 * qJD(4) - t146;
t211 = t65 - t182;
t195 = t125 * t211;
t207 = t101 - t201;
t212 = t141 + t92;
t236 = -t126 * t207 + t129 * (t128 * t212 + t195);
t199 = pkin(6) + pkin(1);
t203 = -t200 - t201;
t205 = t102 - t181;
t55 = t128 * t205;
t215 = t125 * t203 + t55;
t191 = t125 * t205;
t214 = t128 * t203 - t191;
t227 = t126 * t214 - t129 * t212;
t235 = qJ(2) * t215 - t199 * t227;
t234 = pkin(3) * t215;
t233 = pkin(7) * t214;
t232 = pkin(7) * t215;
t210 = t65 + t182;
t85 = -t101 + t200;
t226 = t126 * t210 + t129 * (-t125 * t85 + t55);
t206 = t101 + t201;
t225 = pkin(3) * t206;
t223 = t128 * t85 + t191;
t218 = t129 * t206;
t216 = t211 * qJ(5);
t213 = -t125 * t212 + t128 * t211;
t132 = qJD(1) ^ 2;
t209 = t199 * t132;
t169 = qJD(2) * qJD(1);
t120 = 0.2e1 * t169;
t122 = qJDD(1) * qJ(2);
t127 = sin(qJ(1));
t130 = cos(qJ(1));
t151 = t130 * g(1) + t127 * g(2);
t144 = -t122 + t151;
t140 = t120 - t144;
t148 = -t108 + t164;
t149 = -t107 + t163;
t39 = t149 * pkin(3) + t148 * pkin(7) + t140 - t209;
t131 = qJD(3) ^ 2;
t154 = pkin(3) * t126 - pkin(7) * t129;
t142 = t132 * t154;
t161 = t127 * g(1) - t130 * g(2);
t150 = qJDD(2) - t161;
t174 = t132 * qJ(2);
t138 = t150 - t174;
t88 = -t199 * qJDD(1) + t138;
t70 = t129 * g(3) - t126 * t88;
t53 = -t131 * pkin(3) + qJDD(3) * pkin(7) - t126 * t142 - t70;
t22 = t125 * t39 + t128 * t53;
t68 = t103 * pkin(4) - t105 * qJ(5);
t160 = -t102 * qJ(5) + t103 * t68 - t22;
t202 = -(t200 + t208) * pkin(4) + qJ(5) * t204 - t160;
t198 = t141 * pkin(4);
t197 = pkin(4) * t128;
t194 = t125 * t210;
t69 = t126 * g(3) + t129 * t88;
t52 = qJDD(3) * pkin(3) + t131 * pkin(7) - t129 * t142 + t69;
t193 = t125 * t52;
t187 = t128 * t52;
t184 = qJ(5) * t128;
t183 = qJDD(1) * pkin(1);
t180 = t116 * t125;
t179 = t116 * t128;
t123 = t126 ^ 2;
t178 = t123 * t132;
t124 = t129 ^ 2;
t177 = t124 * t132;
t166 = t126 * t132 * t129;
t176 = t126 * (qJDD(3) + t166);
t175 = t129 * (qJDD(3) - t166);
t173 = t123 + t124;
t171 = qJD(5) * t116;
t167 = t126 * t181;
t165 = t103 * t179;
t162 = -qJ(5) * t125 - pkin(3);
t21 = t125 * t53 - t128 * t39;
t6 = t125 * t21 + t128 * t22;
t113 = 0.2e1 * t171;
t157 = t113 - t160;
t82 = t105 * t180;
t156 = t129 * (t128 * t65 - t82) + t167;
t155 = t103 * t180 - t128 * t141;
t11 = -pkin(4) * t200 + t157;
t12 = -t102 * pkin(4) - qJ(5) * t200 + t105 * t68 + qJDD(5) + t21;
t153 = -pkin(4) * t12 + qJ(5) * t11;
t152 = -pkin(4) * t210 - qJ(5) * t44;
t147 = t125 * t22 - t128 * t21;
t38 = -t126 * t70 + t129 * t69;
t145 = qJ(2) + t154;
t139 = (-t103 * t125 - t105 * t128) * t116;
t137 = t129 * (t82 - t165) + t126 * t102;
t136 = t129 * (t125 * t141 + t165) - t167;
t135 = pkin(4) * t205 + qJ(5) * t203 - t12;
t134 = -pkin(4) * t92 + 0.2e1 * qJD(5) * t105 + t52;
t133 = t134 + t216;
t110 = t173 * qJDD(1);
t109 = t118 - 0.2e1 * t164;
t106 = 0.2e1 * t163 + t168;
t93 = -t138 + t183;
t83 = t144 - 0.2e1 * t169 + t209;
t80 = -t176 + t129 * (-t131 - t177);
t79 = t126 * (-t131 - t178) + t175;
t50 = (qJD(4) + t116) * t103 + t146;
t45 = (-qJD(4) + t116) * t105 + t158;
t41 = t128 * t210;
t40 = t105 * t179 + t125 * t65;
t26 = t128 * t45 + t194;
t25 = -t128 * t44 + t194;
t23 = -t125 * t44 - t41;
t18 = t129 * t50 - t240;
t16 = t129 * t211 + t240;
t15 = t126 * t26 + t218;
t14 = t126 * t25 + t218;
t13 = t133 - t198;
t10 = qJ(5) * t206 + t12;
t9 = (-t200 + t206) * pkin(4) + t157;
t8 = (-t212 - t141) * pkin(4) + t133;
t7 = t134 - t198 + 0.2e1 * t216;
t4 = t126 * t6 + t129 * t52;
t3 = t128 * t11 + t125 * t12;
t2 = t125 * t11 - t128 * t12;
t1 = t126 * t3 + t129 * t13;
t5 = [0, 0, 0, 0, 0, qJDD(1), t161, t151, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t150 - 0.2e1 * t183, t120 + 0.2e1 * t122 - t151, pkin(1) * t93 + qJ(2) * (-t132 * pkin(1) + t140), -t148 * t129, -t129 * t106 - t126 * t109, t175 - t126 * (t131 - t177), t149 * t126, t129 * (-t131 + t178) - t176, 0, qJ(2) * t106 - t126 * t83 - t199 * t79, qJ(2) * t109 - t129 * t83 - t199 * t80, t199 * t110 - t173 * t174 - t38, -qJ(2) * t83 - t199 * t38, t156, -t236, t226, t136, -t239, t137, t129 * (-t193 - t232) - t126 * (t21 - t234) + t235, t129 * (-t187 + t243) - t126 * (t22 + t244) - t241 - t199 * t18, -t129 * t147 + t145 * (t125 * t45 - t41) - t199 * t15, t145 * t147 - t199 * t4, t156, t226, t236, t137, t239, t136, t129 * (-t125 * t8 - t184 * t212 - t232) - t126 * (-t135 - t234) + t235, t129 * (-pkin(7) * t23 + t128 * t10 - t125 * t9) - t126 * (-pkin(3) * t23 - t152) + qJ(2) * t23 - t199 * t14, t129 * (-pkin(4) * t195 + t128 * t7 - t243) - t126 * (-0.2e1 * t171 - t202 - t244) + t241 - t199 * t16, t129 * (-pkin(7) * t2 + (-pkin(4) * t125 + t184) * t13) - t126 * (-pkin(3) * t2 - t153) + qJ(2) * t2 - t199 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t132, -t93, 0, 0, 0, 0, 0, 0, t79, t80, -t110, t38, 0, 0, 0, 0, 0, 0, t227, t18, t15, t4, 0, 0, 0, 0, 0, 0, t227, t14, t16, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, (-t123 + t124) * t132, t118, -t166, -t168, qJDD(3), t69, t70, 0, 0, t40, t213, t223, t155, t238, t139, -pkin(3) * t212 + t187 + t233, pkin(3) * t50 - t193 - t242, pkin(7) * t26 + t225 + t6, pkin(3) * t52 + pkin(7) * t6, t40, t223, -t213, t139, -t238, t155, t128 * t8 + t162 * t212 + t233, pkin(7) * t25 + t125 * t10 + t128 * t9 + t225, t242 + t125 * t7 + (pkin(3) + t197) * t211, pkin(7) * t3 + (-t162 + t197) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, t207, t210, -t181, -t44, t102, -t21, -t22, 0, 0, t181, t210, -t207, t102, t44, -t181, t135, t152, t113 + t202, t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, t210, t208, t12;];
tauJ_reg = t5;
