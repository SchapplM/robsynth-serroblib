% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRP6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:18
% EndTime: 2019-12-05 16:52:26
% DurationCPUTime: 2.20s
% Computational Cost: add. (3267->224), mult. (6735->274), div. (0->0), fcn. (4507->8), ass. (0->144)
t146 = sin(qJ(3));
t149 = cos(qJ(3));
t148 = cos(qJ(4));
t145 = sin(qJ(4));
t137 = qJDD(3) + qJDD(4);
t170 = qJD(2) * t146;
t108 = -t148 * t149 * qJD(2) + t145 * t170;
t173 = t146 * t148;
t110 = (t145 * t149 + t173) * qJD(2);
t177 = t110 * t108;
t77 = -t177 - t137;
t184 = t145 * t77;
t107 = t110 ^ 2;
t138 = qJD(3) + qJD(4);
t194 = t138 ^ 2;
t198 = -t107 - t194;
t51 = t148 * t198 + t184;
t179 = t148 * t77;
t53 = -t145 * t198 + t179;
t18 = t146 * t51 - t149 * t53;
t233 = pkin(6) * t18;
t147 = sin(qJ(2));
t150 = cos(qJ(2));
t178 = t108 * t138;
t169 = qJD(2) * qJD(3);
t165 = t149 * t169;
t168 = t146 * qJDD(2);
t115 = t165 + t168;
t166 = t146 * t169;
t167 = t149 * qJDD(2);
t160 = -t166 + t167;
t61 = -t108 * qJD(4) + t148 * t115 + t145 * t160;
t200 = t61 - t178;
t232 = t147 * t18 + t150 * t200;
t231 = pkin(3) * t51;
t230 = pkin(7) * t51;
t229 = pkin(7) * t53;
t195 = t108 ^ 2;
t93 = t195 - t194;
t228 = t146 * (-t148 * t93 - t184) + t149 * (-t145 * t93 + t179);
t197 = -t177 + t137;
t183 = t145 * t197;
t196 = -t194 - t195;
t203 = t148 * t196 - t183;
t67 = t148 * t197;
t205 = t145 * t196 + t67;
t213 = -t146 * t205 + t149 * t203;
t227 = pkin(6) * t213;
t101 = t110 * t138;
t163 = -t145 * t115 + t148 * t160;
t161 = qJD(4) * t110 - t163;
t199 = t101 + t161;
t224 = t147 * t213 - t150 * t199;
t62 = -t195 - t107;
t223 = pkin(2) * t62;
t222 = pkin(3) * t62;
t221 = pkin(3) * t205;
t220 = pkin(7) * t203;
t219 = pkin(7) * t205;
t216 = t150 * t62;
t215 = qJ(5) * t200;
t214 = t146 * (t145 * t200 + t148 * t199) - t149 * (-t145 * t199 + t148 * t200);
t94 = -t107 + t194;
t212 = t146 * (-t145 * t94 + t67) + t149 * (t148 * t94 + t183);
t201 = t178 + t61;
t79 = t107 - t195;
t193 = 2 * qJD(5);
t192 = pkin(4) * t145;
t191 = pkin(4) * t148;
t152 = qJD(2) ^ 2;
t125 = t146 * t152 * t149;
t121 = qJDD(3) + t125;
t142 = sin(pkin(8));
t143 = cos(pkin(8));
t119 = -g(1) * t142 + g(2) * t143;
t120 = -g(1) * t143 - g(2) * t142;
t171 = -g(3) + qJDD(1);
t98 = t150 * t120 + t147 * t171;
t85 = -t152 * pkin(2) + qJDD(2) * pkin(6) + t98;
t65 = -t149 * t119 + t146 * t85;
t153 = (-t115 + t165) * pkin(7) + t121 * pkin(3) - t65;
t123 = qJD(3) * pkin(3) - pkin(7) * t170;
t140 = t149 ^ 2;
t135 = t140 * t152;
t66 = t146 * t119 + t149 * t85;
t38 = -pkin(3) * t135 + t160 * pkin(7) - qJD(3) * t123 + t66;
t21 = t145 * t38 - t148 * t153;
t23 = t145 * t153 + t148 * t38;
t5 = t145 * t23 - t148 * t21;
t190 = t146 * t5;
t78 = pkin(4) * t108 - qJ(5) * t110;
t162 = t137 * qJ(5) - t108 * t78 + t138 * t193 + t23;
t14 = -pkin(4) * t194 + t162;
t16 = -t137 * pkin(4) - qJ(5) * t194 + t110 * t78 + qJDD(5) + t21;
t189 = -pkin(4) * t16 + qJ(5) * t14;
t43 = -t101 + t161;
t188 = -pkin(4) * t201 - qJ(5) * t43;
t186 = t145 * t201;
t97 = -t120 * t147 + t150 * t171;
t84 = -qJDD(2) * pkin(2) - t152 * pkin(6) - t97;
t55 = -t160 * pkin(3) - pkin(7) * t135 + t123 * t170 + t84;
t185 = t145 * t55;
t180 = t148 * t55;
t176 = t138 * t145;
t175 = t138 * t148;
t174 = t146 * t121;
t172 = t149 * (qJDD(3) - t125);
t164 = -qJ(5) * t145 - pkin(3);
t6 = t145 * t21 + t148 * t23;
t30 = t146 * t65 + t149 * t66;
t3 = t14 * t145 - t148 * t16;
t1 = -t146 * t3 + t149 * (t14 * t148 + t145 * t16);
t116 = -0.2e1 * t166 + t167;
t159 = -pkin(4) * t198 - qJ(5) * t77 + t14;
t158 = pkin(4) * t197 + qJ(5) * t196 - t16;
t157 = t161 * pkin(4) - t215 + t55;
t156 = t146 * (t108 * t175 + t145 * t161) + t149 * (t108 * t176 - t148 * t161);
t155 = t110 * t193 - t157;
t92 = t110 * t176;
t154 = t146 * t92 + (-t108 * t173 + t149 * (-t108 * t145 - t110 * t148)) * t138;
t151 = qJD(3) ^ 2;
t139 = t146 ^ 2;
t134 = t139 * t152;
t118 = t134 + t135;
t117 = (t139 + t140) * qJDD(2);
t114 = 0.2e1 * t165 + t168;
t88 = -t172 - t146 * (-t134 - t151);
t87 = t149 * (-t135 - t151) - t174;
t45 = (-qJD(4) + t138) * t110 + t163;
t39 = t148 * t201;
t27 = t148 * t45 + t186;
t26 = -t148 * t43 + t186;
t25 = t145 * t45 - t39;
t24 = -t145 * t43 - t39;
t22 = t146 * (t148 * t61 - t92) + t149 * (t110 * t175 + t145 * t61);
t17 = (pkin(4) * t138 - (2 * qJD(5))) * t110 + t157;
t12 = -qJ(5) * t62 + t16;
t11 = (-t194 - t62) * pkin(4) + t162;
t10 = (-t199 - t101) * pkin(4) + t155;
t9 = -pkin(4) * t101 + t155 + t215;
t8 = -t146 * t25 + t149 * t27;
t7 = -t146 * t24 + t149 * t26;
t2 = t149 * t6 - t190;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t171, 0, 0, 0, 0, 0, 0, qJDD(2) * t150 - t147 * t152, -qJDD(2) * t147 - t150 * t152, 0, t147 * t98 + t150 * t97, 0, 0, 0, 0, 0, 0, t116 * t150 + t147 * t87, -t114 * t150 + t147 * t88, t117 * t147 + t118 * t150, t147 * t30 - t150 * t84, 0, 0, 0, 0, 0, 0, t224, -t232, t147 * t8 - t216, t147 * t2 - t150 * t55, 0, 0, 0, 0, 0, 0, t224, t147 * t7 - t216, t232, t1 * t147 - t150 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t97, -t98, 0, 0, (t115 + t165) * t146, t114 * t149 + t116 * t146, t174 + t149 * (-t134 + t151), t116 * t149, t146 * (t135 - t151) + t172, 0, pkin(2) * t116 + pkin(6) * t87 - t149 * t84, -pkin(2) * t114 + pkin(6) * t88 + t146 * t84, pkin(2) * t118 + pkin(6) * t117 + t30, -pkin(2) * t84 + pkin(6) * t30, t22, -t214, t212, t156, -t228, t154, t146 * (t185 - t219) + t149 * (-pkin(3) * t199 - t180 + t220) - pkin(2) * t199 + t227, t146 * (t180 - t230) + t149 * (-pkin(3) * t200 + t185 + t229) - pkin(2) * t200 - t233, t146 * (-pkin(7) * t25 - t5) + t149 * (pkin(7) * t27 - t222 + t6) - t223 + pkin(6) * t8, -pkin(7) * t190 + t149 * (-pkin(3) * t55 + pkin(7) * t6) - pkin(2) * t55 + pkin(6) * t2, t22, t212, t214, t154, t228, t156, t146 * (-t10 * t145 - t219) + t149 * (t10 * t148 + t220) + t227 + (-qJ(5) * t173 + t149 * t164 - pkin(2)) * t199, t146 * (-pkin(7) * t24 - t11 * t145 + t12 * t148) + t149 * (pkin(7) * t26 + t11 * t148 + t12 * t145 - t222) - t223 + pkin(6) * t7, t146 * (t148 * t9 + t230) + t149 * (t145 * t9 - t229) + t233 + (-t146 * t192 + t149 * (pkin(3) + t191) + pkin(2)) * t200, (t146 * (-qJ(5) * t148 + t192) + t149 * (t164 - t191) - pkin(2)) * t17 + (pkin(6) + pkin(7)) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, t134 - t135, t168, t125, t167, qJDD(3), -t65, -t66, 0, 0, t177, t79, t201, -t177, -t43, t137, -t21 + t221, -t23 + t231, pkin(3) * t25, pkin(3) * t5, t177, t201, -t79, t137, t43, -t177, t158 + t221, pkin(3) * t24 + t188, t159 - t231, pkin(3) * t3 + t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t79, t201, -t177, -t43, t137, -t21, -t23, 0, 0, t177, t201, -t79, t137, t43, -t177, t158, t188, t159, t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t197, t201, t198, t16;];
tauJ_reg = t4;
