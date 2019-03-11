% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPPRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:02
% EndTime: 2019-03-09 01:32:07
% DurationCPUTime: 1.66s
% Computational Cost: add. (3522->236), mult. (7218->315), div. (0->0), fcn. (4950->8), ass. (0->132)
t178 = 2 * qJD(3);
t166 = sin(qJ(5));
t90 = sin(pkin(10));
t126 = t166 * t90;
t92 = cos(pkin(10));
t96 = cos(qJ(5));
t67 = t96 * t92 - t126;
t77 = sin(pkin(9)) * pkin(1) + qJ(3);
t177 = qJD(1) * t77;
t66 = t166 * t92 + t96 * t90;
t102 = qJD(1) * t66;
t174 = qJD(6) + t102;
t103 = t66 * qJD(4);
t137 = qJD(1) * t92;
t76 = -cos(pkin(9)) * pkin(1) - pkin(2) - qJ(4);
t172 = t76 * qJD(1);
t65 = qJD(3) + t172;
t43 = -t90 * qJD(2) + t92 * t65;
t37 = -pkin(7) * t137 + t43;
t138 = qJD(1) * t90;
t44 = t92 * qJD(2) + t90 * t65;
t38 = -pkin(7) * t138 + t44;
t105 = t166 * t38 - t96 * t37;
t11 = -qJD(1) * t103 - qJD(5) * t105;
t63 = t66 * qJD(5);
t48 = qJD(1) * t63;
t129 = t96 * t137;
t122 = qJD(1) * t126;
t72 = qJD(5) * t122;
t49 = qJD(5) * t129 - t72;
t88 = qJD(3) * qJD(1);
t26 = t49 * pkin(5) + t48 * pkin(8) + t88;
t16 = t166 * t37 + t96 * t38;
t14 = qJD(5) * pkin(8) + t16;
t68 = qJD(4) + t177;
t55 = pkin(4) * t138 + t68;
t61 = -t122 + t129;
t22 = pkin(5) * t102 - t61 * pkin(8) + t55;
t94 = sin(qJ(6));
t95 = cos(qJ(6));
t6 = t95 * t14 + t94 * t22;
t2 = -qJD(6) * t6 - t94 * t11 + t95 * t26;
t170 = t174 * t6 + t2;
t116 = t94 * t14 - t95 * t22;
t1 = -t116 * qJD(6) + t95 * t11 + t94 * t26;
t118 = t116 * t174 + t1;
t176 = t67 * qJD(4);
t124 = t95 * t174;
t146 = t94 * t49;
t175 = -t124 * t174 - t146;
t42 = t94 * qJD(5) + t95 * t61;
t136 = qJD(6) * t42;
t24 = -t94 * t48 + t136;
t140 = t90 ^ 2 + t92 ^ 2;
t171 = t140 * qJD(4);
t135 = qJD(6) * t94;
t149 = t63 * t95;
t106 = t67 * t135 + t149;
t46 = t95 * t49;
t169 = -t106 * t174 + t67 * t46;
t168 = t61 ^ 2;
t127 = 0.2e1 * t88;
t167 = -pkin(7) + t76;
t108 = t176 * qJD(1);
t12 = t16 * qJD(5) + t108;
t57 = t167 * t90;
t58 = t167 * t92;
t30 = t166 * t57 - t96 * t58;
t165 = t12 * t30;
t164 = t12 * t66;
t163 = t12 * t67;
t162 = t12 * t94;
t133 = t95 * qJD(5);
t23 = -qJD(6) * t133 + t61 * t135 + t95 * t48;
t161 = t23 * t94;
t160 = t24 * t95;
t40 = t94 * t61 - t133;
t159 = t40 * t102;
t158 = t40 * t94;
t157 = t42 * t40;
t156 = t42 * t61;
t155 = t42 * t94;
t154 = t42 * t95;
t153 = t49 * t66;
t152 = t61 * t40;
t151 = t61 * t102;
t150 = t63 * t94;
t148 = t67 * t95;
t21 = t94 * t24;
t134 = qJD(6) * t95;
t144 = -t40 * t134 - t21;
t143 = -t24 * t148 + t40 * t149;
t64 = t67 * qJD(5);
t142 = -t23 * t66 + t42 * t64;
t141 = t102 * t63 - t67 * t49;
t51 = t64 * qJD(5);
t131 = t42 * t150;
t69 = t90 * pkin(4) + t77;
t125 = t94 * t174;
t123 = qJD(6) * t66 + qJD(1);
t120 = -t116 * t95 + t6 * t94;
t119 = -t116 * t94 - t6 * t95;
t115 = -t66 * t24 - t64 * t40;
t29 = t66 * pkin(5) - t67 * pkin(8) + t69;
t31 = t166 * t58 + t96 * t57;
t9 = t95 * t29 - t94 * t31;
t10 = t94 * t29 + t95 * t31;
t114 = t154 + t158;
t113 = t43 * t92 + t44 * t90;
t112 = -t67 * t48 - t61 * t63;
t111 = -t66 * t48 + t61 * t64;
t110 = t102 * t64 + t153;
t109 = t46 + (-t102 * t94 - t135) * t174;
t107 = t67 * t134 - t150;
t13 = -qJD(5) * pkin(5) + t105;
t104 = -pkin(8) * t49 + t13 * t174;
t101 = t114 * qJD(6) - t161;
t100 = t105 * t63 + t11 * t66 + t16 * t64 - t163;
t99 = -t107 * t174 - t67 * t146;
t98 = -t120 * qJD(6) + t1 * t95 - t2 * t94;
t97 = qJD(1) ^ 2;
t56 = t102 ^ 2;
t50 = t63 * qJD(5);
t34 = t61 * pkin(5) + pkin(8) * t102;
t32 = t64 * pkin(5) + t63 * pkin(8) + qJD(3);
t20 = t31 * qJD(5) + t176;
t19 = -qJD(5) * t30 - t103;
t8 = -t105 * t95 + t94 * t34;
t7 = t105 * t94 + t95 * t34;
t4 = -t10 * qJD(6) - t94 * t19 + t95 * t32;
t3 = t9 * qJD(6) + t95 * t19 + t94 * t32;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, t177 * t178, 0, 0, 0, 0, 0, 0, t90 * t127, t92 * t127, 0.2e1 * qJD(1) * t171 (t68 + t177) * qJD(3) + (-t140 * t172 - t113) * qJD(4), t112, -t111 + t141, -t50, t110, -t51, 0, -t20 * qJD(5) + t102 * t178 + t69 * t49 + t55 * t64, -t19 * qJD(5) - t69 * t48 - t55 * t63 + (qJD(1) * t67 + t61) * qJD(3), -t102 * t19 + t20 * t61 - t30 * t48 - t31 * t49 - t100, t11 * t31 + t165 + t105 * t20 + t16 * t19 + (qJD(1) * t69 + t55) * qJD(3), -t106 * t42 - t23 * t148, t131 + (t161 + (-t154 + t158) * qJD(6)) * t67 + t143, t142 + t169, t107 * t40 + t67 * t21, t115 + t99, t174 * t64 + t153, t107 * t13 - t116 * t64 + t67 * t162 + t174 * t4 + t2 * t66 + t20 * t40 + t30 * t24 + t9 * t49, -t1 * t66 - t10 * t49 - t106 * t13 + t12 * t148 - t174 * t3 + t20 * t42 - t30 * t23 - t6 * t64, -t10 * t24 + t9 * t23 - t3 * t40 - t4 * t42 + t120 * t63 + (t119 * qJD(6) - t1 * t94 - t2 * t95) * t67, t1 * t10 - t116 * t4 + t13 * t20 + t2 * t9 + t6 * t3 + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50, t111 + t141, t105 * t64 + t11 * t67 - t16 * t63 + t164, 0, 0, 0, 0, 0, 0, -t115 + t99, t142 - t169, t101 * t67 - t131 + t143, t119 * t63 + t13 * t64 + t98 * t67 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t177 * qJD(1), 0, 0, 0, 0, 0, 0, -t97 * t90, -t97 * t92, 0 (-t68 - t171) * qJD(1), 0, 0, 0, 0, 0, 0, -qJD(1) * t102 - t50, -qJD(1) * t61 - t51, -t110 - t112, -t55 * qJD(1) + t100, 0, 0, 0, 0, 0, 0, -t66 * t146 - t67 * t24 + t63 * t40 + (-t123 * t95 - t64 * t94) * t174, -t66 * t46 + t67 * t23 + t63 * t42 + (t123 * t94 - t64 * t95) * t174 (-t40 * t95 + t155) * t64 + t114 * qJD(1) + (t101 - t160) * t66, -t120 * qJD(1) - t119 * t64 + t13 * t63 + t98 * t66 - t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140 * t97, t113 * qJD(1) + t88, 0, 0, 0, 0, 0, 0, -t72 + (t61 + t129) * qJD(5), -0.2e1 * t102 * qJD(5), -t56 - t168, t102 * t16 - t105 * t61 + t88, 0, 0, 0, 0, 0, 0, t109 - t152, -t156 + t175 (t23 - t159) * t95 + t42 * t125 + t144, t118 * t94 - t13 * t61 + t170 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, -t56 + t168, 0, -t151, t72 + (t61 - t129) * qJD(5), 0, -t55 * t61 - t108 (qJD(4) + t55) * t102, 0, 0, t42 * t124 - t161 (-t23 - t159) * t95 - t174 * t155 + t144, -t156 - t175, t40 * t125 - t160, t109 + t152, -t174 * t61, -pkin(5) * t24 - t12 * t95 - t16 * t40 + t116 * t61 + (-pkin(8) * t134 - t7) * t174 + t104 * t94, pkin(5) * t23 + t162 - t16 * t42 + t6 * t61 + (pkin(8) * t135 + t8) * t174 + t104 * t95, t8 * t40 + t7 * t42 + ((-t24 + t136) * pkin(8) + t118) * t95 + ((qJD(6) * t40 - t23) * pkin(8) - t170) * t94, -t12 * pkin(5) + pkin(8) * t98 + t116 * t7 - t13 * t16 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, -t40 ^ 2 + t42 ^ 2, t174 * t40 - t23, -t157, t174 * t42 - t24, t49, -t13 * t42 + t170, t13 * t40 - t118, 0, 0;];
tauc_reg  = t5;
