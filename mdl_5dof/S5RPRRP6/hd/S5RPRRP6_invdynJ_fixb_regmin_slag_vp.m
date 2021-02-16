% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:27
% EndTime: 2021-01-15 18:08:36
% DurationCPUTime: 1.85s
% Computational Cost: add. (2037->308), mult. (4207->396), div. (0->0), fcn. (2637->10), ass. (0->163)
t100 = sin(pkin(8));
t85 = pkin(1) * t100 + pkin(6);
t73 = t85 * qJDD(1);
t214 = -qJD(2) * qJD(3) - t73;
t104 = sin(qJ(3));
t107 = cos(qJ(3));
t75 = t85 * qJD(1);
t42 = qJD(2) * t107 - t104 * t75;
t213 = t42 * qJD(3);
t97 = qJ(1) + pkin(8);
t90 = sin(t97);
t91 = cos(t97);
t207 = g(1) * t91 + g(2) * t90;
t170 = qJD(1) * t104;
t212 = qJD(4) * t170 - qJDD(3);
t103 = sin(qJ(4));
t106 = cos(qJ(4));
t158 = t104 * qJDD(1);
t169 = qJD(1) * t107;
t25 = ((qJD(4) + t169) * qJD(3) + t158) * t103 + t212 * t106;
t162 = t106 * qJD(3);
t62 = t103 * t170 - t162;
t81 = -qJD(4) + t169;
t194 = t62 * t81;
t160 = qJD(1) * qJD(3);
t146 = t107 * t160;
t24 = -qJD(4) * t162 + (-t146 - t158) * t106 + t212 * t103;
t211 = -t24 + t194;
t167 = qJD(3) * t103;
t64 = t106 * t170 + t167;
t193 = t64 * t81;
t210 = t25 - t193;
t164 = qJD(4) * t103;
t152 = t81 * t164;
t94 = t107 * qJDD(1);
t58 = t104 * t160 + qJDD(4) - t94;
t121 = t106 * t58 + t152;
t209 = t121 * t104;
t171 = qJDD(2) - g(3);
t208 = t171 * t107;
t175 = t103 * t104;
t174 = t103 * t107;
t37 = t106 * t91 + t90 * t174;
t39 = t106 * t90 - t91 * t174;
t206 = -g(1) * t39 + g(2) * t37 + g(3) * t175;
t204 = t64 ^ 2;
t43 = t104 * qJD(2) + t107 * t75;
t34 = qJD(3) * pkin(7) + t43;
t101 = cos(pkin(8));
t86 = -pkin(1) * t101 - pkin(2);
t54 = -pkin(3) * t107 - pkin(7) * t104 + t86;
t35 = t54 * qJD(1);
t12 = -t103 * t34 + t106 * t35;
t9 = -qJ(5) * t64 + t12;
t8 = -pkin(4) * t81 + t9;
t203 = -t9 + t8;
t202 = pkin(4) * t62;
t197 = t58 * pkin(4);
t196 = pkin(4) * t103;
t195 = g(3) * t107;
t192 = qJ(5) + pkin(7);
t148 = t107 * t162;
t173 = t104 * t106;
t191 = -t62 * t148 - t25 * t173;
t134 = pkin(3) * t104 - pkin(7) * t107;
t67 = t134 * qJD(1);
t190 = t103 * t67 + t106 * t42;
t163 = qJD(4) * t106;
t68 = t134 * qJD(3);
t189 = t103 * t68 + t54 * t163;
t143 = qJD(4) * t192;
t149 = t103 * t169;
t161 = t106 * qJD(5);
t188 = qJ(5) * t149 - t103 * t143 + t161 - t190;
t172 = t106 * t107;
t123 = pkin(4) * t104 - qJ(5) * t172;
t52 = t106 * t67;
t187 = -t123 * qJD(1) - t106 * t143 - t52 + (-qJD(5) + t42) * t103;
t166 = qJD(3) * t104;
t182 = t103 * t85;
t186 = t106 * t68 + t166 * t182;
t66 = t85 * t172;
t185 = t103 * t54 + t66;
t184 = t207 * t173;
t98 = t104 ^ 2;
t183 = -t107 ^ 2 + t98;
t181 = t106 * t81;
t180 = t24 * qJ(5);
t179 = t25 * qJ(5);
t76 = qJD(1) * t86;
t178 = qJD(3) * t62;
t177 = qJD(3) * t85;
t176 = qJD(4) * t62;
t165 = qJD(3) * t107;
t157 = t107 * qJDD(2);
t156 = t214 * t104 - t75 * t165;
t155 = t81 * t167;
t154 = t81 * t162;
t153 = t64 * t165;
t144 = -qJD(5) - t202;
t33 = -qJD(3) * pkin(3) - t42;
t23 = -t144 + t33;
t151 = t23 * t163;
t150 = pkin(6) + t196;
t147 = t104 * t163;
t17 = qJDD(3) * pkin(7) + t104 * qJDD(2) + t107 * t73 + t213;
t26 = qJD(1) * t68 + t54 * qJDD(1);
t142 = t103 * t26 + t106 * t17 + t35 * t163 - t34 * t164;
t140 = t64 * t147;
t128 = -qJDD(3) * pkin(3) - t156;
t18 = t128 - t157;
t139 = -qJD(4) * pkin(7) * t81 + t18;
t138 = -g(1) * t37 - g(2) * t39;
t38 = t103 * t91 - t90 * t172;
t40 = t103 * t90 + t91 * t172;
t137 = -g(1) * t38 - g(2) * t40;
t135 = g(1) * t90 - g(2) * t91;
t105 = sin(qJ(1));
t108 = cos(qJ(1));
t133 = g(1) * t105 - g(2) * t108;
t13 = t103 * t35 + t106 * t34;
t10 = -qJ(5) * t62 + t13;
t132 = t10 * t106 - t103 * t8;
t131 = t10 * t103 + t106 * t8;
t89 = pkin(4) * t106 + pkin(3);
t129 = t104 * t192 + t107 * t89;
t126 = -pkin(2) - t129;
t125 = t133 * pkin(1);
t124 = t207 * t104;
t122 = -t103 * t58 + t81 * t163;
t120 = -qJD(1) * t76 + t207;
t119 = -pkin(7) * t58 - t81 * t33;
t118 = 0.2e1 * t76 * qJD(3) - qJDD(3) * t85;
t117 = t25 * pkin(4) + qJDD(5) + t128;
t116 = g(1) * t40 - g(2) * t38 + g(3) * t173 - t142;
t109 = qJD(3) ^ 2;
t115 = -0.2e1 * qJDD(1) * t86 - t109 * t85 + t135;
t21 = t106 * t26;
t114 = -t13 * qJD(4) - t103 * t17 + t21;
t1 = -t64 * qJD(5) + t114 + t180 + t197;
t2 = -qJD(5) * t62 + t142 - t179;
t113 = -t131 * qJD(4) - t1 * t103 + t2 * t106;
t112 = t114 + t206;
t110 = qJD(1) ^ 2;
t84 = g(3) * t174;
t78 = t192 * t106;
t77 = t192 * t103;
t72 = qJDD(3) * t107 - t104 * t109;
t71 = qJDD(3) * t104 + t107 * t109;
t57 = t62 ^ 2;
t49 = (t85 + t196) * t104;
t46 = t64 * t166;
t45 = t106 * t54;
t28 = pkin(4) * t149 + t43;
t27 = t85 * t165 + (t103 * t165 + t147) * pkin(4);
t22 = -qJ(5) * t175 + t185;
t19 = -qJ(5) * t173 + t45 + (-pkin(4) - t182) * t107;
t7 = t117 - t157;
t6 = (-qJ(5) * qJD(4) - t177) * t173 + (-qJD(5) * t104 + (-qJ(5) * qJD(3) - qJD(4) * t85) * t107) * t103 + t189;
t5 = -t104 * t161 + t123 * qJD(3) + (-t66 + (qJ(5) * t104 - t54) * t103) * qJD(4) + t186;
t4 = (-t25 + t155) * t107 + (t122 + t178) * t104;
t3 = t46 + (t24 + t154) * t107 - t209;
t11 = [qJDD(1), t133, g(1) * t108 + g(2) * t105, (t100 ^ 2 + t101 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t125, qJDD(1) * t98 + 0.2e1 * t104 * t146, 0.2e1 * t104 * t94 - 0.2e1 * t183 * t160, t71, t72, 0, t118 * t104 + t115 * t107, -t115 * t104 + t118 * t107, t64 * t148 + (-t24 * t106 - t64 * t164) * t104, -t140 + (-t153 + (t24 + t176) * t104) * t103 + t191, t46 + (t24 - t154) * t107 + t209, (t25 + t155) * t107 + (t122 - t178) * t104, -t107 * t58 - t166 * t81, -(-t164 * t54 + t186) * t81 + t45 * t58 + (t62 * t177 - t21 + (t81 * t85 + t34) * t163 + (qJD(3) * t33 + qJD(4) * t35 - t58 * t85 + t17) * t103) * t107 + (qJD(3) * t12 + t18 * t103 + t163 * t33 + t25 * t85) * t104 + t137, t189 * t81 - t185 * t58 + (-t85 * t152 + (t106 * t33 + t64 * t85) * qJD(3) + t142) * t107 + (-t33 * t164 + t18 * t106 - t85 * t24 + (-t85 * t181 - t13) * qJD(3)) * t104 + t138, t19 * t58 + t49 * t25 + t27 * t62 - t5 * t81 + (t167 * t23 - t1) * t107 + (qJD(3) * t8 + t7 * t103 + t151) * t104 + t137, -t22 * t58 - t49 * t24 + t27 * t64 + t6 * t81 + (t162 * t23 + t2) * t107 + (-qJD(3) * t10 + t7 * t106 - t164 * t23) * t104 + t138, t19 * t24 - t22 * t25 - t5 * t64 - t6 * t62 - t131 * t165 + (-t132 * qJD(4) - t1 * t106 - t103 * t2 + t135) * t104, t1 * t19 + t10 * t6 + t2 * t22 + t23 * t27 + t7 * t49 + t8 * t5 + t125 + (-g(1) * t150 + g(2) * t126) * t91 + (-g(1) * t126 - g(2) * t150) * t90; 0, 0, 0, t171, 0, 0, 0, 0, 0, t72, -t71, 0, 0, 0, 0, 0, t4, t3, t4, t3, t140 + (t153 + (-t24 + t176) * t104) * t103 + t191, -g(3) + (qJD(3) * t132 - t7) * t107 + (qJD(3) * t23 + t113) * t104; 0, 0, 0, 0, -t104 * t110 * t107, t183 * t110, t158, t94, qJDD(3), t43 * qJD(3) + t120 * t104 + t156 + t208, t213 + (qJD(3) * t75 - t171) * t104 + (t120 + t214) * t107, -t24 * t103 - t64 * t181, -t210 * t103 + t211 * t106, (-t104 * t64 + t81 * t172) * qJD(1) - t122, (t104 * t62 - t174 * t81) * qJD(1) + t121, t81 * t170, -t12 * t170 - pkin(3) * t25 - t43 * t62 + t52 * t81 + (-t139 - t195) * t106 + (-t42 * t81 + t119) * t103 + t184, pkin(3) * t24 - t190 * t81 + t13 * t170 - t43 * t64 + t84 + t119 * t106 + (-t124 + t139) * t103, -t8 * t170 - t89 * t25 - t28 * t62 - t77 * t58 - t187 * t81 + (-t7 - t195) * t106 + (-t23 * t169 + (t23 + t202) * qJD(4)) * t103 + t184, t151 + t89 * t24 - t28 * t64 - t78 * t58 + t84 + t188 * t81 + (t10 * t104 - t172 * t23) * qJD(1) + (pkin(4) * qJD(4) * t64 - t124 + t7) * t103, -g(3) * t104 - t24 * t77 - t25 * t78 - t187 * t64 - t188 * t62 + (qJD(1) * t131 - t207) * t107 + t113, t2 * t78 - t1 * t77 - t7 * t89 - g(3) * t129 + t187 * t8 + (pkin(4) * t164 - t28) * t23 + t188 * t10 + t207 * (t104 * t89 - t107 * t192); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t62, -t57 + t204, -t24 - t194, -t193 - t25, t58, -t13 * t81 - t33 * t64 + t112, -t12 * t81 + t33 * t62 + t116, 0.2e1 * t197 + t180 - t10 * t81 + (t144 - t23) * t64 + t112, -t204 * pkin(4) + t179 - t9 * t81 + (qJD(5) + t23) * t62 + t116, t24 * pkin(4) - t203 * t62, t203 * t10 + (-t23 * t64 + t1 + t206) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, t211, -t57 - t204, t10 * t62 + t8 * t64 + t117 - t124 - t208;];
tau_reg = t11;
