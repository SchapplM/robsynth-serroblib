% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:16
% EndTime: 2019-12-05 17:10:22
% DurationCPUTime: 1.75s
% Computational Cost: add. (1382->212), mult. (2253->297), div. (0->0), fcn. (1716->14), ass. (0->153)
t111 = qJ(2) + qJ(3);
t101 = cos(t111);
t202 = g(3) * t101;
t112 = sin(pkin(9));
t113 = cos(pkin(9));
t147 = g(1) * t113 + g(2) * t112;
t99 = sin(t111);
t217 = t147 * t99;
t218 = -t217 + t202;
t119 = cos(qJ(4));
t174 = qJD(4) * t119;
t216 = qJD(5) * t119 + t174;
t107 = qJD(2) + qJD(3);
t114 = sin(qJ(5));
t115 = sin(qJ(4));
t118 = cos(qJ(5));
t182 = t118 * t115;
t58 = t114 * t119 + t182;
t45 = t58 * t107;
t120 = cos(qJ(3));
t199 = t120 * pkin(2);
t121 = cos(qJ(2));
t170 = t121 * qJD(1);
t116 = sin(qJ(3));
t117 = sin(qJ(2));
t177 = qJD(1) * t117;
t91 = t116 * t177;
t49 = t120 * t170 - t91;
t214 = qJD(3) * t199 - t49;
t153 = qJD(3) * t177;
t176 = qJD(3) * t116;
t169 = qJD(1) * qJD(2);
t212 = t117 * qJDD(1) + t121 * t169;
t97 = t121 * qJDD(1);
t50 = qJDD(2) * pkin(2) - t117 * t169 + t97;
t89 = qJD(2) * pkin(2) + t170;
t151 = t89 * t176 + (t153 - t50) * t120 + t212 * t116;
t105 = qJDD(2) + qJDD(3);
t201 = t105 * pkin(3);
t17 = t151 - t201;
t213 = t17 + t202;
t211 = t216 * t118;
t57 = t116 * t117 - t120 * t121;
t209 = t107 * t57;
t59 = t116 * t121 + t120 * t117;
t48 = t59 * qJD(1);
t150 = pkin(2) * t176 - t48;
t106 = qJD(4) + qJD(5);
t122 = qJD(4) ^ 2;
t94 = t116 * pkin(2) + pkin(7);
t95 = -pkin(3) - t199;
t208 = t105 * t95 + t150 * t107 + t122 * t94;
t207 = (qJD(3) * t89 + t212) * t120 + t116 * t50;
t206 = -pkin(7) - pkin(8);
t92 = g(3) * t99;
t205 = -pkin(8) - t94;
t200 = t107 * pkin(3);
t185 = t114 * t115;
t162 = t107 * t185;
t181 = t118 * t119;
t43 = -t107 * t181 + t162;
t198 = t45 * t43;
t175 = qJD(4) * t115;
t163 = pkin(4) * t175;
t197 = t163 + t150;
t196 = t106 * t209;
t42 = t116 * t89 + t120 * t177;
t36 = t107 * pkin(7) + t42;
t157 = pkin(8) * t107 + t36;
t24 = t157 * t119;
t194 = t118 * t24;
t193 = t42 * t107;
t104 = qJDD(4) + qJDD(5);
t56 = -t181 + t185;
t192 = t56 * t104;
t191 = t58 * t104;
t190 = t101 * t112;
t189 = t101 * t113;
t188 = t107 * t115;
t110 = qJ(4) + qJ(5);
t100 = cos(t110);
t187 = t112 * t100;
t186 = t113 * t100;
t183 = t115 * t105;
t180 = t119 * t105;
t179 = qJDD(1) - g(3);
t108 = t115 ^ 2;
t178 = -t119 ^ 2 + t108;
t173 = qJD(5) * t114;
t172 = qJD(5) * t118;
t41 = t120 * t89 - t91;
t35 = -t41 - t200;
t167 = t213 * t115 + t35 * t174;
t166 = t119 * t217 + t35 * t175;
t96 = -t119 * pkin(4) - pkin(3);
t161 = qJD(4) * t206;
t160 = t107 * t174;
t155 = qJD(4) * t205;
t86 = t116 * t153;
t152 = g(1) * t189 + g(2) * t190 + t86 + t92;
t149 = -t42 + t163;
t148 = t114 * t183 - t118 * t180;
t146 = g(1) * t112 - g(2) * t113;
t23 = t157 * t115;
t29 = t107 * t59;
t145 = -t57 * t105 - t29 * t107;
t22 = qJD(4) * pkin(4) - t23;
t143 = -t114 * t22 - t194;
t102 = t119 * pkin(8);
t52 = t119 * t94 + t102;
t142 = qJD(5) * t52 + t214 * t115 - t119 * t155;
t80 = t119 * pkin(7) + t102;
t141 = qJD(5) * t80 - t115 * t41 - t119 * t161;
t51 = t205 * t115;
t140 = -qJD(5) * t51 - t115 * t155 - t214 * t119;
t79 = t206 * t115;
t139 = -qJD(5) * t79 - t115 * t161 + t119 * t41;
t137 = t107 * t175 - t180;
t10 = t137 * pkin(4) + t17;
t25 = t96 * t107 - t41;
t27 = t106 * t58;
t136 = t10 * t56 - t100 * t218 + t25 * t27;
t134 = pkin(7) * t122 - t193 - t201;
t133 = t122 * t59 - t145;
t132 = -pkin(7) * qJDD(4) + (t41 - t200) * qJD(4);
t131 = 0.2e1 * t209 * qJD(4) - qJDD(4) * t59;
t130 = -t151 - t218;
t26 = t106 * t56;
t98 = sin(t110);
t128 = t10 * t58 + t218 * t98 - t25 * t26;
t12 = t105 * t182 - t106 * t162 + t211 * t107 + t114 * t180;
t127 = -qJDD(4) * t94 + (t107 * t95 - t214) * qJD(4);
t16 = t105 * pkin(7) + t207 - t86;
t126 = t147 * t101 - t35 * t107 - t16 + t92;
t4 = -t36 * t174 + qJDD(4) * pkin(4) - t115 * t16 + (-t160 - t183) * pkin(8);
t125 = t24 * t173 + t100 * t92 + t25 * t43 + (-t24 * t106 - t4) * t114 - g(1) * (-t101 * t186 - t112 * t98) - g(2) * (-t101 * t187 + t113 * t98);
t5 = -t137 * pkin(8) + t119 * t16 - t36 * t175;
t124 = -g(1) * (-t98 * t189 + t187) - g(2) * (-t98 * t190 - t186) + t143 * qJD(5) - t114 * t5 + t118 * t4 - t25 * t45 + t98 * t92;
t123 = qJD(2) ^ 2;
t103 = t107 ^ 2;
t72 = t96 - t199;
t71 = qJDD(4) * t119 - t122 * t115;
t70 = qJDD(4) * t115 + t122 * t119;
t46 = t108 * t105 + 0.2e1 * t115 * t160;
t32 = -0.2e1 * t178 * t107 * qJD(4) + 0.2e1 * t115 * t180;
t20 = -t27 * t106 - t192;
t19 = -t26 * t106 + t191;
t18 = -t43 ^ 2 + t45 ^ 2;
t13 = t27 * t107 + t148;
t8 = t43 * t106 + t12;
t2 = t12 * t58 - t45 * t26;
t1 = -t12 * t56 - t58 * t13 + t26 * t43 - t45 * t27;
t3 = [t179, 0, t121 * qJDD(2) - t123 * t117, -qJDD(2) * t117 - t123 * t121, 0, t145, -t59 * t105 + t107 * t209, 0, 0, 0, 0, 0, t131 * t115 - t133 * t119, t133 * t115 + t131 * t119, 0, 0, 0, 0, 0, t57 * t13 + t29 * t43 + t58 * t196 + ((t114 * t175 + t115 * t173 - t211) * t106 - t191) * t59, t57 * t12 + t29 * t45 - t56 * t196 + (-(-t216 * t114 - t115 * t172 - t118 * t175) * t106 + t192) * t59; 0, qJDD(2), -g(3) * t121 + t147 * t117 + t97, -t179 * t117 + t147 * t121, t105, t48 * t107 + (t105 * t120 - t107 * t176) * pkin(2) + t130, t49 * t107 + (-pkin(2) * t105 - t50) * t116 + ((-pkin(2) * t107 - t89) * qJD(3) - t212) * t120 + t152, t46, t32, t70, t71, 0, t127 * t115 + (-t213 - t208) * t119 + t166, t127 * t119 + (-t217 + t208) * t115 + t167, t2, t1, t19, t20, 0, (-t114 * t52 + t118 * t51) * t104 + t72 * t13 + t197 * t43 + (t140 * t114 - t142 * t118) * t106 + t136, -(t114 * t51 + t118 * t52) * t104 + t72 * t12 + t197 * t45 + (t142 * t114 + t140 * t118) * t106 + t128; 0, 0, 0, 0, t105, t130 + t193, t41 * t107 + t152 - t207, t46, t32, t70, t71, 0, t132 * t115 + (-t134 - t213) * t119 + t166, t132 * t119 + (-t217 + t134) * t115 + t167, t2, t1, t19, t20, 0, (-t114 * t80 + t118 * t79) * t104 + t96 * t13 + t149 * t43 + (t139 * t114 - t141 * t118) * t106 + t136, -(t114 * t79 + t118 * t80) * t104 + t96 * t12 + t149 * t45 + (t141 * t114 + t139 * t118) * t106 + t128; 0, 0, 0, 0, 0, 0, 0, -t115 * t103 * t119, t178 * t103, t183, t180, qJDD(4), t126 * t115 - t146 * t119, t146 * t115 + t126 * t119, t198, t18, t8, -t148, t104, -(t114 * t23 - t194) * t106 + (t118 * t104 - t106 * t173 - t43 * t188) * pkin(4) + t124, (-qJD(5) * t22 - t23 * t106 - t5) * t118 + (-t114 * t104 - t106 * t172 - t188 * t45) * pkin(4) + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t18, t8, -t148, t104, -t143 * t106 + t124, (-t5 + (-qJD(5) + t106) * t22) * t118 + t125;];
tau_reg = t3;
