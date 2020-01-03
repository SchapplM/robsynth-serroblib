% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:07
% EndTime: 2019-12-31 17:21:12
% DurationCPUTime: 2.24s
% Computational Cost: add. (1817->333), mult. (4247->428), div. (0->0), fcn. (2628->6), ass. (0->165)
t101 = cos(qJ(3));
t102 = cos(qJ(2));
t167 = qJD(1) * t102;
t223 = -t167 + qJD(3);
t99 = sin(qJ(2));
t178 = qJD(1) * t99;
t155 = t101 * t178;
t98 = sin(qJ(3));
t177 = qJD(2) * t98;
t59 = t155 + t177;
t195 = t59 * t223;
t160 = t98 * t178;
t163 = t101 * qJD(2);
t57 = t160 - t163;
t197 = t57 * t223;
t161 = qJD(1) * qJD(2);
t152 = t102 * t161;
t162 = t99 * qJDD(1);
t174 = qJD(3) * t99;
t222 = qJD(1) * t174 - qJDD(2);
t20 = -qJD(3) * t163 + t222 * t98 + (-t152 - t162) * t101;
t21 = t222 * t101 + ((qJD(3) + t167) * qJD(2) + t162) * t98;
t227 = (t20 + t197) * t101 + (t21 + t195) * t98;
t224 = t21 - t195;
t166 = qJD(2) * t102;
t181 = t21 * t101;
t194 = t59 * t98;
t199 = t20 * t98;
t221 = ((t101 * t59 - t57 * t98) * qJD(3) + t181 - t199) * t99 + (t101 * t57 + t194) * t166;
t190 = t102 * pkin(2) + t99 * pkin(6);
t220 = -pkin(1) - t190;
t100 = sin(qJ(1));
t103 = cos(qJ(1));
t141 = g(1) * t103 + g(2) * t100;
t219 = t141 * t99;
t89 = t102 * qJDD(1);
t218 = -t99 * t161 + t89;
t164 = qJD(3) * t102;
t175 = qJD(3) * t98;
t176 = qJD(2) * t99;
t144 = pkin(2) * t99 - pkin(6) * t102;
t62 = t144 * qJD(2);
t15 = pkin(5) * (-t101 * t164 + t98 * t176) + t101 * t62 - t220 * t175;
t67 = -qJD(2) * pkin(2) + pkin(5) * t178;
t19 = pkin(3) * t57 - qJ(4) * t59 + t67;
t55 = qJDD(3) - t218;
t208 = pkin(6) * t55;
t216 = -t19 * t223 + t208;
t123 = t101 * t55 - t175 * t223;
t183 = t102 * t98;
t215 = qJD(1) * (-t183 * t223 - t57 * t99) - t123;
t165 = qJD(3) * t101;
t122 = t165 * t223 + t98 * t55;
t212 = (-t177 * t223 + t21) * t102 + (-qJD(2) * t57 - t122) * t99;
t211 = -0.2e1 * pkin(1);
t210 = t59 ^ 2;
t209 = pkin(3) * t55;
t207 = pkin(6) * t59;
t206 = g(3) * t99;
t205 = g(1) * t100;
t202 = g(2) * t103;
t201 = g(3) * t102;
t52 = t220 * qJD(1);
t86 = pkin(5) * t167;
t68 = qJD(2) * pkin(6) + t86;
t25 = t101 * t68 + t52 * t98;
t18 = qJ(4) * t223 + t25;
t200 = t18 * t223;
t198 = t25 * t223;
t196 = t59 * t57;
t193 = t98 * t99;
t139 = pkin(3) * t98 - qJ(4) * t101;
t192 = -t98 * qJD(4) + t223 * t139 - t86;
t170 = t101 * t103;
t173 = t100 * t101;
t191 = (g(1) * t170 + g(2) * t173) * t99;
t171 = t101 * t102;
t37 = pkin(5) * t171 + t220 * t98;
t189 = t103 * pkin(1) + t100 * pkin(5);
t96 = t99 ^ 2;
t97 = t102 ^ 2;
t188 = t96 - t97;
t187 = t96 + t97;
t186 = t100 * t99;
t61 = t144 * qJD(1);
t185 = t101 * t61;
t184 = t101 * t220;
t182 = t103 * t99;
t180 = t55 * qJ(4);
t179 = pkin(5) * qJDD(1);
t172 = t100 * t102;
t169 = t102 * t103;
t24 = t101 * t52 - t68 * t98;
t168 = qJD(4) - t24;
t159 = t223 * t178;
t105 = qJD(1) ^ 2;
t158 = t99 * t105 * t102;
t157 = t57 ^ 2 - t210;
t156 = pkin(5) * t98 + pkin(3);
t29 = qJD(1) * t62 + qJDD(1) * t220;
t44 = t218 * pkin(5) + qJDD(2) * pkin(6);
t151 = -t101 * t29 + t68 * t165 + t52 * t175 + t98 * t44;
t150 = pkin(2) * t169 + pkin(6) * t182 + t189;
t148 = t99 * t152;
t46 = t98 * t172 + t170;
t48 = t98 * t169 - t173;
t147 = -g(1) * t46 + g(2) * t48;
t47 = t100 * t171 - t103 * t98;
t49 = t100 * t98 + t101 * t169;
t146 = g(1) * t47 - g(2) * t49;
t84 = pkin(5) * t162;
t45 = -qJDD(2) * pkin(2) + pkin(5) * t152 + t84;
t143 = pkin(5) * t57 + t67 * t98;
t142 = pkin(5) * t59 + t101 * t67;
t140 = (qJD(3) * t57 - t20) * pkin(6);
t138 = pkin(3) * t101 + qJ(4) * t98;
t137 = qJD(3) * t67 - t208;
t16 = -pkin(3) * t223 + t168;
t136 = t101 * t16 - t18 * t98;
t134 = -t101 * t24 - t25 * t98;
t128 = pkin(2) + t138;
t127 = pkin(5) + t139;
t126 = -pkin(6) * qJD(3) * t223 - t201;
t5 = t101 * t44 + t52 * t165 - t68 * t175 + t98 * t29;
t120 = -pkin(5) * qJDD(2) + t161 * t211;
t2 = pkin(3) * t21 + qJ(4) * t20 - qJD(4) * t59 + t45;
t119 = t126 - t2;
t118 = t126 - t45;
t117 = pkin(1) * t105 + t141;
t104 = qJD(2) ^ 2;
t115 = pkin(5) * t104 + qJDD(1) * t211 + t202;
t114 = t220 * t205;
t112 = g(1) * t48 + g(2) * t46 + g(3) * t193 - t151;
t111 = t197 * t98 - t181;
t110 = -pkin(6) * t181 - t141 * t102 - t206;
t14 = t98 * t62 + t220 * t165 + (-t99 * t163 - t98 * t164) * pkin(5);
t109 = t19 * t59 + qJDD(4) - t112;
t108 = t21 * t193 + (t99 * t165 + t98 * t166) * t57;
t107 = -g(1) * t49 - g(2) * t47 - t101 * t206 + t5;
t94 = t103 * pkin(5);
t80 = g(1) * t186;
t77 = pkin(6) * t169;
t74 = pkin(6) * t172;
t50 = t98 * t61;
t42 = t127 * t99;
t36 = -pkin(5) * t183 + t184;
t34 = -pkin(5) * t155 + t50;
t33 = pkin(5) * t160 + t185;
t32 = t156 * t102 - t184;
t31 = -qJ(4) * t102 + t37;
t30 = -t102 * t55 + t176 * t223;
t28 = pkin(3) * t59 + qJ(4) * t57;
t27 = -t156 * t178 - t185;
t26 = t50 + (-pkin(5) * t101 + qJ(4)) * t178;
t13 = (t138 * qJD(3) - qJD(4) * t101) * t99 + t127 * t166;
t12 = -pkin(3) * t176 - t15;
t11 = -t20 + t197;
t10 = qJ(4) * t176 - t102 * qJD(4) + t14;
t9 = (-t171 * t223 - t59 * t99) * qJD(1) + t122;
t8 = t101 * t195 - t199;
t7 = -t174 * t194 + (t59 * t166 - t20 * t99) * t101;
t4 = qJDD(4) + t151 - t209;
t3 = (t163 * t223 + t20) * t102 + (qJD(2) * t59 + t123) * t99;
t1 = qJD(4) * t223 + t180 + t5;
t6 = [0, 0, 0, 0, 0, qJDD(1), -t202 + t205, t141, 0, 0, qJDD(1) * t96 + 0.2e1 * t148, -0.2e1 * t188 * t161 + 0.2e1 * t99 * t89, qJDD(2) * t99 + t102 * t104, qJDD(1) * t97 - 0.2e1 * t148, qJDD(2) * t102 - t104 * t99, 0, t120 * t99 + (-t115 + t205) * t102, t102 * t120 + t115 * t99 - t80, 0.2e1 * t179 * t187 - t141, -g(1) * (-pkin(1) * t100 + t94) - g(2) * t189 + (t187 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), t7, -t221, t3, t108, t212, t30, t15 * t223 + t36 * t55 + (t143 * qJD(2) + t151) * t102 + (pkin(5) * t21 + qJD(2) * t24 + t67 * t165 + t45 * t98) * t99 + t146, -t14 * t223 - t37 * t55 + (t142 * qJD(2) + t5) * t102 + (-pkin(5) * t20 - qJD(2) * t25 + t45 * t101 - t67 * t175) * t99 + t147, -t14 * t57 - t15 * t59 + t36 * t20 - t37 * t21 + t80 + t134 * t166 + (-t202 + t101 * t151 - t5 * t98 + (-t101 * t25 + t24 * t98) * qJD(3)) * t99, t5 * t37 + t25 * t14 - t151 * t36 + t24 * t15 - g(1) * t94 - g(2) * t150 - t114 + (t67 * t166 + t45 * t99) * pkin(5), t7, t3, t221, t30, -t212, t108, -t12 * t223 + t13 * t57 + t42 * t21 - t32 * t55 + (t177 * t19 + t4) * t102 + (-qJD(2) * t16 + t165 * t19 + t2 * t98) * t99 + t146, -t10 * t57 + t12 * t59 - t32 * t20 - t31 * t21 + t80 + t136 * t166 + (-t202 - t1 * t98 + t101 * t4 + (-t101 * t18 - t16 * t98) * qJD(3)) * t99, t10 * t223 - t13 * t59 + t42 * t20 + t31 * t55 + (-t163 * t19 - t1) * t102 + (qJD(2) * t18 - t2 * t101 + t175 * t19) * t99 - t147, t1 * t31 + t18 * t10 + t2 * t42 + t19 * t13 + t4 * t32 + t16 * t12 - g(1) * (-t47 * pkin(3) - t46 * qJ(4) + t94) - g(2) * (pkin(3) * t49 + qJ(4) * t48 + t150) - t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158, t188 * t105, t162, t158, t89, qJDD(2), t117 * t99 - t201 - t84, t206 + (t117 - t179) * t102, 0, 0, t8, -t227, t9, t111, -t215, -t159, -pkin(2) * t21 - t33 * t223 + t137 * t98 + t118 * t101 + (-t143 * t102 - t24 * t99) * qJD(1) + t191, pkin(2) * t20 + t34 * t223 + t137 * t101 + (-t142 * t102 + t25 * t99) * qJD(1) + (-t118 - t219) * t98, t33 * t59 + t34 * t57 + (t24 * t167 + t5 + (-t24 + t207) * qJD(3)) * t101 + (t140 + t151 - t198) * t98 + t110, -t45 * pkin(2) - t25 * t34 - t24 * t33 - t67 * t86 - g(1) * (-pkin(2) * t182 + t77) - g(2) * (-pkin(2) * t186 + t74) - g(3) * t190 + (t134 * qJD(3) + t5 * t101 + t151 * t98) * pkin(6), t8, t9, t227, -t159, t215, t111, t119 * t101 - t128 * t21 + t16 * t178 + t192 * t57 - t216 * t98 + t223 * t27 + t191, t26 * t57 - t27 * t59 + (-t16 * t167 + t1 + (t16 + t207) * qJD(3)) * t101 + (t140 + t4 - t200) * t98 + t110, -t18 * t178 - t128 * t20 - t26 * t223 - t192 * t59 + t216 * t101 + (t119 + t219) * t98, -t18 * t26 - t16 * t27 - g(1) * t77 - g(2) * t74 - g(3) * (t102 * t138 + t190) + t192 * t19 + (qJD(3) * t136 + t1 * t101 + t4 * t98) * pkin(6) + (-t2 + t219) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, -t157, t11, -t196, -t224, t55, -t59 * t67 + t112 + t198, t223 * t24 + t57 * t67 - t107, 0, 0, t196, t11, t157, t55, t224, -t196, -t28 * t57 - t109 + t198 + 0.2e1 * t209, pkin(3) * t20 - t21 * qJ(4) + (t18 - t25) * t59 + (t16 - t168) * t57, 0.2e1 * t180 - t19 * t57 + t28 * t59 - (-0.2e1 * qJD(4) + t24) * t223 + t107, t1 * qJ(4) - t4 * pkin(3) - t19 * t28 - t16 * t25 - g(1) * (-pkin(3) * t48 + qJ(4) * t49) - g(2) * (-pkin(3) * t46 + qJ(4) * t47) + t139 * t206 + t168 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 + t196, t11, -t223 ^ 2 - t210, t109 - t200 - t209;];
tau_reg = t6;
