% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPP2
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
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:33
% EndTime: 2021-01-15 15:32:40
% DurationCPUTime: 1.40s
% Computational Cost: add. (1530->269), mult. (3384->350), div. (0->0), fcn. (2375->10), ass. (0->142)
t103 = qJ(4) + pkin(6);
t101 = cos(pkin(8));
t104 = sin(qJ(3));
t106 = cos(qJ(3));
t139 = qJD(3) * t103;
t122 = -t104 * qJD(4) - t106 * t139;
t107 = cos(qJ(2));
t169 = t101 * t106;
t99 = sin(pkin(8));
t133 = -t99 * t104 + t169;
t127 = t133 * t107;
t52 = t106 * qJD(4) - t104 * t139;
t182 = -qJD(1) * t127 + t101 * t52 + t99 * t122;
t100 = sin(pkin(7));
t102 = cos(pkin(7));
t136 = g(1) * t102 + g(2) * t100;
t129 = t136 * t107;
t105 = sin(qJ(2));
t186 = g(3) * t105;
t118 = t129 + t186;
t95 = qJ(3) + pkin(8);
t92 = sin(t95);
t179 = t107 * t92;
t142 = t103 * t104;
t72 = t103 * t106;
t35 = t101 * t72 - t99 * t142;
t65 = t101 * t104 + t99 * t106;
t58 = t65 * qJD(2);
t196 = g(3) * t179 - t35 * qJDD(3) - (qJD(1) * t58 + t136 * t92) * t105;
t54 = t58 ^ 2;
t143 = qJD(2) * t169;
t159 = qJD(2) * t104;
t55 = t99 * t159 - t143;
t195 = -t55 ^ 2 - t54;
t57 = t65 * qJD(3);
t156 = qJD(3) * t104;
t60 = qJD(3) * t169 - t99 * t156;
t149 = t106 * qJDD(2);
t150 = t104 * qJDD(2);
t134 = -t101 * t149 + t99 * t150;
t31 = qJD(2) * t57 + t134;
t152 = qJD(2) * qJD(3);
t141 = t104 * t152;
t120 = t65 * qJDD(2) - t99 * t141;
t140 = t106 * t152;
t32 = t101 * t140 + t120;
t194 = t31 * pkin(4) - t32 * qJ(5);
t94 = g(3) * t107;
t193 = t136 * t105 - t94;
t108 = qJD(3) ^ 2;
t148 = t107 * qJDD(1);
t153 = qJD(1) * qJD(2);
t137 = t105 * t153 - t148;
t191 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t108 + (t136 + t153) * t105 - t137 - t94;
t63 = qJDD(2) * pkin(6) + t105 * qJDD(1) + t107 * t153;
t126 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + t63;
t155 = t105 * qJD(1);
t138 = t103 * qJD(2) + t155;
t132 = qJD(3) * t138;
t14 = qJDD(3) * pkin(3) - t126 * t104 - t106 * t132;
t18 = -t104 * t132 + t126 * t106;
t5 = t101 * t18 + t99 * t14;
t189 = pkin(3) * t104;
t185 = t106 * pkin(3);
t51 = t138 * t106;
t184 = t99 * t51;
t4 = t101 * t14 - t99 * t18;
t154 = t107 * qJD(1);
t183 = -t101 * t122 - t65 * t154 + t99 * t52;
t37 = t101 * t51;
t50 = t138 * t104;
t39 = qJD(3) * pkin(3) - t50;
t22 = t99 * t39 + t37;
t97 = t104 ^ 2;
t181 = -t106 ^ 2 + t97;
t180 = qJD(2) * pkin(2);
t93 = cos(t95);
t178 = t107 * t93;
t173 = qJDD(3) * pkin(4);
t172 = t100 * t105;
t171 = t100 * t107;
t168 = t102 * t105;
t167 = t102 * t106;
t166 = t102 * t107;
t165 = t103 * t107;
t164 = t104 * t107;
t23 = -t99 * t50 + t37;
t163 = t23 * qJD(3);
t24 = -t101 * t50 - t184;
t162 = qJD(5) - t24;
t161 = qJDD(1) - g(3);
t109 = qJD(2) ^ 2;
t160 = t108 + t109;
t158 = qJD(2) * t105;
t151 = qJDD(3) * t104;
t147 = t107 * qJDD(2);
t146 = pkin(3) * t156;
t91 = pkin(2) + t185;
t135 = g(1) * t100 - g(2) * t102;
t21 = t101 * t39 - t184;
t25 = -t60 * t105 - t107 * t58;
t26 = qJD(2) * t127 - t105 * t57;
t48 = t65 * t105;
t49 = t133 * t105;
t128 = -t25 * t58 - t26 * t55 - t49 * t31 + t48 * t32;
t44 = t102 * t93 + t92 * t171;
t46 = -t100 * t93 + t92 * t166;
t125 = g(1) * t46 + g(2) * t44 + t92 * t186 + t4;
t34 = t101 * t142 + t99 * t72;
t123 = -g(3) * t178 - t34 * qJDD(3) + (g(1) * t168 + g(2) * t172) * t93;
t119 = t25 * qJD(3) - t48 * qJDD(3) - t107 * t31 + t55 * t158;
t45 = -t102 * t92 + t93 * t171;
t47 = t100 * t92 + t93 * t166;
t117 = g(1) * t47 + g(2) * t45 + t93 * t186 - t5;
t61 = -t91 * qJD(2) + qJD(4) - t154;
t19 = t55 * pkin(4) - t58 * qJ(5) + t61;
t116 = -t19 * t58 - qJDD(5) + t125;
t81 = -t154 - t180;
t115 = -pkin(6) * qJDD(3) + (t154 + t81 - t180) * qJD(3);
t33 = pkin(3) * t141 - t91 * qJDD(2) + qJDD(4) + t137;
t114 = t26 * qJD(3) + t49 * qJDD(3) + t107 * t32 - t58 * t158;
t113 = -t81 * qJD(2) + t118 - t63;
t112 = 0.2e1 * t58 * qJD(3) + t134;
t111 = -t182 * t55 + t183 * t58 - t35 * t31 + t34 * t32 - t118;
t110 = -t193 + t33;
t96 = qJDD(3) * qJ(5);
t89 = -t101 * pkin(3) - pkin(4);
t87 = t99 * pkin(3) + qJ(5);
t85 = t100 * t185;
t77 = t107 * t91;
t75 = t102 * t165;
t74 = t100 * t165;
t30 = -pkin(4) * t133 - t65 * qJ(5) - t91;
t27 = pkin(3) * t159 + t58 * pkin(4) + t55 * qJ(5);
t20 = qJD(3) * qJ(5) + t22;
t17 = -qJD(3) * pkin(4) + qJD(5) - t21;
t16 = (-t55 + t143) * qJD(3) + t120;
t11 = t57 * pkin(4) - t60 * qJ(5) - t65 * qJD(5) + t146;
t3 = qJDD(5) - t173 - t4;
t2 = -t58 * qJD(5) + t194 + t33;
t1 = qJD(3) * qJD(5) + t5 + t96;
t6 = [t161, 0, -t109 * t105 + t147, -qJDD(2) * t105 - t109 * t107, 0, 0, 0, 0, 0, (-0.2e1 * t141 + t149) * t107 + (-t160 * t106 - t151) * t105, (-qJDD(3) * t105 - 0.2e1 * t107 * t152) * t106 + (t160 * t105 - t147) * t104, t119, -t114, t128, -t33 * t107 + t61 * t158 + t21 * t25 + t22 * t26 - t4 * t48 + t5 * t49 - g(3), t119, t128, t114, t1 * t49 - t2 * t107 + t19 * t158 - t17 * t25 + t20 * t26 + t3 * t48 - g(3); 0, qJDD(2), t148 + t193, -t161 * t105 + t129, t97 * qJDD(2) + 0.2e1 * t104 * t140, 0.2e1 * t104 * t149 - 0.2e1 * t181 * t152, t108 * t106 + t151, qJDD(3) * t106 - t108 * t104, 0, t115 * t104 + t191 * t106, -t191 * t104 + t115 * t106, -t55 * t155 - t91 * t31 - t33 * t133 + t61 * t57 + (t55 * t189 - t183) * qJD(3) + t123, -t91 * t32 + t33 * t65 + t61 * t60 + (t58 * t189 - t182) * qJD(3) + t196, t133 * t5 - t21 * t60 - t22 * t57 - t4 * t65 + t111, t5 * t35 - t4 * t34 - t33 * t91 - g(1) * (-t91 * t168 + t75) - g(2) * (-t91 * t172 + t74) - g(3) * (t105 * t103 + t77) + (t146 - t155) * t61 + t182 * t22 - t183 * t21, t19 * t57 - t2 * t133 + t30 * t31 + (t11 - t155) * t55 - t183 * qJD(3) + t123, t1 * t133 + t17 * t60 - t20 * t57 + t3 * t65 + t111, t182 * qJD(3) - t11 * t58 - t19 * t60 - t2 * t65 - t30 * t32 - t196, t1 * t35 + t2 * t30 + t19 * t11 + t3 * t34 - g(1) * t75 - g(2) * t74 - g(3) * (pkin(4) * t178 + qJ(5) * t179 + t77) + t182 * t20 + t183 * t17 + (-g(3) * t103 - t19 * qJD(1) + t136 * (pkin(4) * t93 + qJ(5) * t92 + t91)) * t105; 0, 0, 0, 0, -t104 * t109 * t106, t181 * t109, t150, t149, qJDD(3), t113 * t104 - t135 * t106, t135 * t104 + t113 * t106, t163 - t61 * t58 + (qJDD(3) * t101 - t55 * t159) * pkin(3) + t125, t24 * qJD(3) + t61 * t55 + (-qJDD(3) * t99 - t58 * t159) * pkin(3) + t117, (t22 - t23) * t58 + (-t21 + t24) * t55 + (-t101 * t32 - t31 * t99) * pkin(3), -g(1) * t85 + t21 * t23 - t22 * t24 + (g(2) * t167 + t4 * t101 + t5 * t99 + (-t61 * qJD(2) + t118) * t104) * pkin(3), t163 - t27 * t55 + (pkin(4) - t89) * qJDD(3) + t116, -t87 * t31 + t89 * t32 + (t20 - t23) * t58 + (t17 - t162) * t55, t87 * qJDD(3) - t19 * t55 + t27 * t58 + t96 + (0.2e1 * qJD(5) - t24) * qJD(3) - t117, t1 * t87 + t3 * t89 - t19 * t27 - t17 * t23 - g(1) * (-t102 * pkin(3) * t164 - t46 * pkin(4) + t47 * qJ(5) + t85) - g(2) * (-t44 * pkin(4) + t45 * qJ(5) + (-t100 * t164 - t167) * pkin(3)) + t162 * t20 - (-pkin(4) * t92 + qJ(5) * t93 - t189) * t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t16, t195, t21 * t58 + t22 * t55 + t110, t112, t195, -t16, t20 * t55 + (-qJD(5) - t17) * t58 + t110 + t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t55 - qJDD(3), (t55 + t143) * qJD(3) + t120, -t54 - t108, -t20 * qJD(3) - t116 - t173;];
tau_reg = t6;
