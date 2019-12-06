% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:38
% EndTime: 2019-12-05 15:51:42
% DurationCPUTime: 1.54s
% Computational Cost: add. (1035->227), mult. (2544->345), div. (0->0), fcn. (2207->12), ass. (0->135)
t97 = cos(qJ(4));
t140 = t97 * qJD(2);
t74 = -qJD(5) + t140;
t179 = t74 + qJD(5);
t89 = sin(pkin(5));
t151 = qJD(2) * t89;
t129 = qJD(1) * t151;
t138 = t89 * qJDD(1);
t95 = sin(qJ(2));
t98 = cos(qJ(2));
t178 = t98 * t129 + t95 * t138;
t94 = sin(qJ(4));
t144 = qJD(5) * t94;
t177 = -qJD(2) * t144 + qJDD(4);
t137 = t94 * qJDD(2);
t93 = sin(qJ(5));
t96 = cos(qJ(5));
t31 = ((qJD(5) + t140) * qJD(4) + t137) * t93 - t177 * t96;
t87 = sin(pkin(10));
t90 = cos(pkin(10));
t52 = t95 * t87 - t98 * t90;
t115 = t98 * t87 + t95 * t90;
t92 = cos(pkin(5));
t47 = t115 * t92;
t88 = sin(pkin(9));
t91 = cos(pkin(9));
t116 = -t88 * t47 - t91 * t52;
t117 = t91 * t47 - t88 * t52;
t162 = t89 * t97;
t46 = t115 * t89;
t32 = t46 * t94 - t92 * t97;
t110 = g(1) * (-t116 * t94 + t88 * t162) + g(2) * (-t117 * t94 - t91 * t162) - g(3) * t32;
t120 = pkin(4) * t94 - pkin(8) * t97;
t71 = t92 * qJDD(1) + qJDD(3);
t156 = t97 * t71;
t70 = t98 * t138;
t44 = qJDD(2) * pkin(2) - t95 * t129 + t70;
t18 = t178 * t90 + t87 * t44;
t16 = qJDD(2) * pkin(7) + t18;
t152 = qJD(1) * t89;
t132 = t98 * t152;
t61 = qJD(2) * pkin(2) + t132;
t134 = t95 * t152;
t65 = t90 * t134;
t37 = t87 * t61 + t65;
t35 = qJD(2) * pkin(7) + t37;
t73 = t92 * qJD(1) + qJD(3);
t22 = t97 * t35 + t94 * t73;
t2 = -qJDD(4) * pkin(4) + t22 * qJD(4) + t94 * t16 - t156;
t176 = (pkin(8) * qJD(5) + t120 * qJD(2)) * t74 - t110 - t2;
t141 = t96 * qJD(4);
t131 = t97 * t141;
t108 = -t93 * t144 + t131;
t136 = qJD(2) * qJD(4);
t127 = t94 * t136;
t83 = t97 * qJDD(2);
t51 = qJDD(5) - t83 + t127;
t157 = t96 * t51;
t175 = -t108 * t74 + t94 * t157;
t40 = t87 * t132 + t65;
t114 = -t97 * pkin(4) - t94 * pkin(8) - pkin(3);
t171 = t90 * pkin(2);
t50 = t114 - t171;
t60 = t120 * qJD(4);
t174 = t50 * t51 + (t40 - t60) * t74;
t111 = t52 * t92;
t25 = -t91 * t111 - t115 * t88;
t28 = t88 * t111 - t115 * t91;
t45 = t52 * t89;
t109 = g(1) * t28 + g(2) * t25 - g(3) * t45;
t77 = t87 * pkin(2) + pkin(7);
t166 = t74 * t77;
t20 = qJD(4) * pkin(8) + t22;
t172 = (t20 + t166) * qJD(5) - t109;
t30 = qJD(2) * t131 + qJD(5) * t141 + t96 * t137 + t177 * t93;
t170 = t30 * t93;
t150 = qJD(2) * t94;
t54 = t93 * t150 - t141;
t168 = t54 * t74;
t142 = t93 * qJD(4);
t56 = t96 * t150 + t142;
t167 = t56 * t74;
t165 = t74 * t96;
t164 = t74 * t97;
t163 = t89 * t94;
t161 = t92 * t95;
t160 = t92 * t98;
t159 = t94 * t71;
t85 = t94 ^ 2;
t153 = -t97 ^ 2 + t85;
t149 = qJD(4) * t54;
t148 = qJD(4) * t77;
t147 = qJD(4) * t94;
t146 = qJD(4) * t97;
t145 = qJD(5) * t74;
t118 = t94 * t35 - t97 * t73;
t19 = -qJD(4) * pkin(4) + t118;
t143 = t19 * qJD(5);
t139 = qJDD(1) - g(3);
t133 = t74 * t142;
t125 = t56 * t147 - t30 * t97;
t64 = t87 * t134;
t36 = t90 * t61 - t64;
t34 = -qJD(2) * pkin(3) - t36;
t123 = -qJD(2) * t34 - t16;
t17 = -t178 * t87 + t90 * t44;
t23 = t114 * qJD(2) - t36;
t4 = t96 * t20 + t93 * t23;
t119 = t93 * t20 - t96 * t23;
t33 = t46 * t97 + t92 * t94;
t10 = t33 * t96 + t45 * t93;
t9 = -t33 * t93 + t45 * t96;
t112 = t96 * t145 - t93 * t51;
t43 = t90 * t132 - t64;
t78 = -pkin(3) - t171;
t106 = -qJDD(4) * t77 + (qJD(2) * t78 + t34 + t43) * qJD(4);
t105 = -g(1) * t116 - g(2) * t117 - g(3) * t46 + t50 * t145;
t104 = -pkin(8) * t51 + (-t19 + t118) * t74;
t1 = qJDD(4) * pkin(8) - t118 * qJD(4) + t97 * t16 + t159;
t103 = qJD(4) * t19 + qJD(5) * t23 - t43 * t74 - t51 * t77 + t1;
t99 = qJD(4) ^ 2;
t102 = -qJD(2) * t40 + t77 * t99 + t109 - t17 + (-pkin(3) + t78) * qJDD(2);
t101 = -g(1) * (-t88 * t160 - t91 * t95) - g(2) * (t91 * t160 - t88 * t95) - g(3) * t89 * t98;
t100 = qJD(2) ^ 2;
t67 = qJDD(4) * t97 - t99 * t94;
t66 = qJDD(4) * t94 + t99 * t97;
t42 = t52 * t151;
t41 = qJD(2) * t46;
t14 = t116 * t97 + t88 * t163;
t12 = t117 * t97 - t91 * t163;
t8 = -t32 * qJD(4) - t42 * t97;
t7 = t33 * qJD(4) - t42 * t94;
t6 = qJD(2) * t60 + t114 * qJDD(2) - t17;
t5 = t96 * t6;
t3 = [t139, 0, (qJDD(2) * t98 - t100 * t95) * t89, (-qJDD(2) * t95 - t100 * t98) * t89, -t17 * t45 + t18 * t46 - t36 * t41 - t37 * t42 + t71 * t92 - g(3), 0, 0, 0, 0, 0, -t45 * t83 - t7 * qJD(4) - t32 * qJDD(4) + (t45 * t147 - t41 * t97) * qJD(2), t45 * t137 - t8 * qJD(4) - t33 * qJDD(4) + (t45 * t146 + t41 * t94) * qJD(2), 0, 0, 0, 0, 0, -(-t10 * qJD(5) + t41 * t96 - t8 * t93) * t74 + t9 * t51 + t7 * t54 + t32 * t31, (qJD(5) * t9 + t41 * t93 + t8 * t96) * t74 - t10 * t51 + t7 * t56 + t32 * t30; 0, qJDD(2), t70 + t101, -g(1) * (t88 * t161 - t91 * t98) - g(2) * (-t91 * t161 - t88 * t98) - t139 * t95 * t89, t36 * t40 - t37 * t43 + (t17 * t90 + t18 * t87 + t101) * pkin(2), t85 * qJDD(2) + 0.2e1 * t97 * t127, -0.2e1 * t153 * t136 + 0.2e1 * t94 * t83, t66, t67, 0, -t102 * t97 + t106 * t94, t102 * t94 + t106 * t97, t30 * t96 * t94 + t108 * t56, (-t54 * t96 - t56 * t93) * t146 + (-t170 - t31 * t96 + (t54 * t93 - t56 * t96) * qJD(5)) * t94, t125 + t175, (t31 + t133) * t97 + (t112 - t149) * t94, -t74 * t147 - t51 * t97, t174 * t96 + t105 * t93 + (t103 * t93 + t54 * t148 + t172 * t96 - t5) * t97 + (t96 * t143 + t2 * t93 + t77 * t31 - t43 * t54 + (-t93 * t166 - t119) * qJD(4)) * t94, -t174 * t93 + t105 * t96 + (t56 * t148 + t103 * t96 + (-t172 + t6) * t93) * t97 + (-t93 * t143 + t2 * t96 + t77 * t30 - t43 * t56 + (-t77 * t165 - t4) * qJD(4)) * t94; 0, 0, 0, 0, -g(3) * t92 + (-g(1) * t88 + g(2) * t91) * t89 + t71, 0, 0, 0, 0, 0, t67, -t66, 0, 0, 0, 0, 0, (-t31 + t133) * t97 + (t112 + t149) * t94, t125 - t175; 0, 0, 0, 0, 0, -t94 * t100 * t97, t153 * t100, t137, t83, qJDD(4), t123 * t94 - t110 + t156, g(1) * t14 + g(2) * t12 + g(3) * t33 + t123 * t97 - t159, -t56 * t165 + t170, (t30 + t168) * t96 + (-t31 + t167) * t93, (t96 * t164 - t94 * t56) * qJD(2) - t112, t93 * t145 + t157 + (-t93 * t164 + t94 * t54) * qJD(2), t74 * t150, -pkin(4) * t31 + t104 * t93 + t119 * t150 + t176 * t96 - t22 * t54, -pkin(4) * t30 + t104 * t96 + t4 * t150 - t176 * t93 - t22 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, t30 - t168, -t167 - t31, t51, -t93 * t1 + t5 - t19 * t56 - g(1) * (-t14 * t93 - t28 * t96) - g(2) * (-t12 * t93 - t25 * t96) - g(3) * t9 - t179 * t4, -t96 * t1 - t93 * t6 + t19 * t54 - g(1) * (-t14 * t96 + t28 * t93) - g(2) * (-t12 * t96 + t25 * t93) + g(3) * t10 + t179 * t119;];
tau_reg = t3;
