% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:56
% EndTime: 2019-12-31 18:12:59
% DurationCPUTime: 1.20s
% Computational Cost: add. (1472->250), mult. (3345->270), div. (0->0), fcn. (2369->8), ass. (0->139)
t109 = sin(qJ(1));
t110 = cos(qJ(1));
t149 = g(1) * t109 - g(2) * t110;
t164 = qJDD(1) * pkin(1);
t129 = -qJDD(2) + t164 + t149;
t197 = (qJDD(3) * qJ(4) + qJD(3) * qJD(4));
t182 = g(1) * t110;
t136 = g(2) * t109 + t182;
t108 = sin(qJ(3));
t184 = cos(qJ(3));
t148 = qJD(3) * t184;
t156 = qJD(3) * t108;
t104 = sin(pkin(7));
t155 = qJD(1) * qJD(2);
t177 = pkin(6) + qJ(2);
t189 = t177 * qJDD(1) + t155;
t41 = t189 * t104;
t105 = cos(pkin(7));
t42 = t189 * t105;
t70 = t177 * t104;
t66 = qJD(1) * t70;
t71 = t177 * t105;
t67 = qJD(1) * t71;
t144 = t108 * t41 + t66 * t148 + t67 * t156 - t184 * t42;
t101 = pkin(7) + qJ(3);
t94 = sin(t101);
t95 = cos(t101);
t117 = -g(3) * t94 - t136 * t95 - t144;
t176 = -t108 * t67 - t184 * t66;
t196 = t176 * qJD(3) - t117;
t159 = -qJD(4) + t176;
t179 = pkin(3) + qJ(5);
t143 = t179 * qJDD(3);
t147 = qJDD(1) * t184;
t154 = t104 * qJDD(1);
t134 = -t105 * t147 + t108 * t154;
t65 = t184 * t104 + t108 * t105;
t62 = t65 * qJD(3);
t32 = qJD(1) * t62 + t134;
t151 = t184 * t105;
t137 = qJD(1) * t151;
t162 = t108 * t104;
t150 = qJD(1) * t162;
t57 = -t137 + t150;
t195 = t32 * qJ(5) + t57 * qJD(5);
t174 = t95 * pkin(3) + t94 * qJ(4);
t194 = 2 * t197;
t188 = t57 ^ 2;
t59 = t65 * qJD(1);
t52 = t59 ^ 2;
t193 = -t188 - t52;
t191 = qJ(2) * qJDD(1);
t190 = -t32 * pkin(4) + qJDD(5);
t22 = (qJD(2) * t104 + qJD(3) * t71) * t108 - qJD(2) * t151 + t70 * t148;
t153 = t105 * qJDD(1);
t152 = qJD(3) * t137 + t104 * t147 + t108 * t153;
t31 = qJD(3) * t150 - t152;
t187 = t31 * pkin(4);
t185 = t57 * pkin(4);
t180 = t59 * t57;
t178 = pkin(4) + t177;
t34 = -t108 * t66 + t184 * t67;
t173 = qJ(4) * t32;
t172 = t109 * t94;
t171 = t109 * t95;
t170 = t110 * t94;
t169 = t110 * t95;
t168 = t57 * qJ(4);
t167 = t95 * qJ(5);
t166 = t104 ^ 2 + t105 ^ 2;
t165 = qJ(4) * t110;
t163 = qJDD(3) * pkin(3);
t160 = t34 * qJD(3);
t19 = -t59 * pkin(4) + t176;
t158 = qJD(4) - t19;
t20 = t34 - t185;
t157 = -qJD(5) - t20;
t92 = t105 * pkin(2) + pkin(1);
t146 = t166 * qJD(1) ^ 2;
t145 = t108 * t42 + t67 * t148 - t66 * t156 + t184 * t41;
t142 = qJDD(3) - t180;
t141 = g(2) * (pkin(3) * t169 + t110 * t92 + t94 * t165);
t140 = 0.2e1 * t166;
t139 = -g(1) * t172 + g(2) * t170;
t138 = g(1) * t171 - g(2) * t169;
t27 = -qJD(3) * qJ(4) - t34;
t135 = -qJDD(4) - t145;
t61 = t104 * t156 - t105 * t148;
t133 = t61 * qJ(4) - t65 * qJD(4);
t5 = t144 - t197;
t131 = -t65 * qJ(4) - t92;
t130 = -t92 - t174;
t35 = t108 * t71 + t184 * t70;
t36 = -t108 * t70 + t184 * t71;
t127 = g(1) * t170 + g(2) * t172 - g(3) * t95 - t145;
t69 = -t92 * qJD(1) + qJD(2);
t68 = -t92 * qJDD(1) + qJDD(2);
t126 = -qJDD(4) + t127;
t125 = t129 + t164;
t124 = -t59 * qJ(4) + t69;
t123 = t22 * qJD(3) - t36 * qJDD(3) + t139;
t23 = t65 * qJD(2) + t36 * qJD(3);
t122 = -t23 * qJD(3) - t35 * qJDD(3) + t138;
t21 = t57 * pkin(3) + t124;
t121 = t21 * t59 - t126;
t120 = t32 * pkin(3) + t31 * qJ(4) + t68;
t119 = t140 * t155 - t136;
t8 = t179 * t57 + t124;
t118 = t8 * t59 - t126 - t187;
t4 = -t59 * qJD(4) + t120;
t116 = t120 - t149;
t16 = 0.2e1 * t59 * qJD(3) + t134;
t114 = -t8 * t57 + t117 + t190;
t111 = qJD(3) ^ 2;
t74 = t95 * t165;
t72 = qJ(4) * t171;
t64 = -t151 + t162;
t45 = qJD(3) * t57;
t43 = -t111 - t52;
t30 = t64 * pkin(3) + t131;
t29 = t59 * pkin(3) + t168;
t26 = -qJD(3) * pkin(3) - t159;
t25 = -t64 * pkin(4) + t36;
t24 = t65 * pkin(4) + t35;
t18 = t179 * t64 + t131;
t17 = t62 * pkin(3) + t133;
t15 = (t57 - t150) * qJD(3) + t152;
t14 = (t57 + t150) * qJD(3) - t152;
t13 = t179 * t59 + t168;
t12 = qJD(5) - t27 - t185;
t11 = -t179 * qJD(3) + t158;
t10 = -t61 * pkin(4) + t23;
t9 = -t62 * pkin(4) - t22;
t7 = t64 * qJD(5) + t179 * t62 + t133;
t6 = -t135 - t163;
t3 = -t5 + t190;
t2 = -qJD(3) * qJD(5) - t135 - t143 - t187;
t1 = t4 + t195;
t28 = [qJDD(1), t149, t136, t125 * t105, -t125 * t104, t140 * t191 + t119, t129 * pkin(1) + (t166 * t191 + t119) * qJ(2), -t31 * t65 - t59 * t61, t31 * t64 - t65 * t32 + t61 * t57 - t59 * t62, -t61 * qJD(3) + t65 * qJDD(3), -t62 * qJD(3) - t64 * qJDD(3), 0, -t92 * t32 + t69 * t62 + t68 * t64 + t122, t92 * t31 - t69 * t61 + t68 * t65 + t123, t22 * t57 + t23 * t59 - t26 * t61 + t27 * t62 - t35 * t31 - t36 * t32 + t5 * t64 + t6 * t65 - t136, -t17 * t57 - t21 * t62 - t30 * t32 - t4 * t64 - t122, -t17 * t59 + t21 * t61 + t30 * t31 - t4 * t65 - t123, t4 * t30 + t21 * t17 - t5 * t36 + t27 * t22 + t6 * t35 + t26 * t23 - t177 * t182 - t141 + (-g(1) * t130 - g(2) * t177) * t109, t10 * t59 - t11 * t61 - t12 * t62 + t2 * t65 - t24 * t31 - t25 * t32 - t3 * t64 - t9 * t57 - t136, t9 * qJD(3) + t25 * qJDD(3) - t1 * t65 + t18 * t31 - t7 * t59 + t8 * t61 - t139, -t10 * qJD(3) - t24 * qJDD(3) + t1 * t64 + t18 * t32 + t7 * t57 + t8 * t62 + t138, t1 * t18 + t8 * t7 + t2 * t24 + t11 * t10 + t3 * t25 + t12 * t9 - t141 + (-g(1) * t178 - g(2) * t167) * t110 + (-g(1) * (t130 - t167) - g(2) * t178) * t109; 0, 0, 0, -t153, t154, -t146, -qJ(2) * t146 - t129, 0, 0, 0, 0, 0, t16, -t14, t193, -t16, t31 + t45, -t27 * t57 + (-qJD(4) - t26) * t59 + t116, t193, t14, t16, t12 * t57 + (-qJD(4) - t11) * t59 + t116 + t195; 0, 0, 0, 0, 0, 0, 0, t180, -t188 + t52, t15, -t134, qJDD(3), -t69 * t59 + t127 + t160, t69 * t57 + t196, pkin(3) * t31 - t173 + (-t27 - t34) * t59 + (t26 + t159) * t57, t29 * t57 + t121 - t160 - 0.2e1 * t163, -t21 * t57 + t29 * t59 + t194 - t196, -t5 * qJ(4) - t6 * pkin(3) - t21 * t29 - t26 * t34 - g(1) * (-pkin(3) * t170 + t74) - g(2) * (-pkin(3) * t172 + t72) - g(3) * t174 + t159 * t27, -t173 + t179 * t31 + (t12 + t157) * t59 + (t11 - t158) * t57, -t19 * qJD(3) + t13 * t59 + t114 + t194, -t13 * t57 + (0.2e1 * qJD(5) + t20) * qJD(3) + 0.2e1 * t143 - t118, t3 * qJ(4) - t8 * t13 - g(1) * t74 - g(2) * t72 - g(3) * (t167 + t174) + t158 * t12 + t157 * t11 + (t136 * t94 - t2) * t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 + t45, t142, t43, t27 * qJD(3) + t121 - t163, t15, t43, -t142, -t143 + (-qJD(5) - t12) * qJD(3) + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, qJDD(3) + t180, -t188 - t111, t11 * qJD(3) + t114 + t197;];
tau_reg = t28;
