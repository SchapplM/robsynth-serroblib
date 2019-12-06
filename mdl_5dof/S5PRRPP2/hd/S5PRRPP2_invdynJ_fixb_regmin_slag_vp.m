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
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 16:10:12
% EndTime: 2019-12-05 16:10:17
% DurationCPUTime: 1.33s
% Computational Cost: add. (1342->239), mult. (2970->315), div. (0->0), fcn. (2065->10), ass. (0->129)
t162 = qJ(4) + pkin(6);
t128 = qJD(3) * t162;
t96 = sin(qJ(3));
t98 = cos(qJ(3));
t111 = -t96 * qJD(4) - t98 * t128;
t93 = cos(pkin(8));
t166 = t93 * t98;
t91 = sin(pkin(8));
t62 = t91 * t96 - t166;
t99 = cos(qJ(2));
t115 = t62 * t99;
t50 = t98 * qJD(4) - t96 * t128;
t160 = qJD(1) * t115 + t91 * t111 + t93 * t50;
t92 = sin(pkin(7));
t94 = cos(pkin(7));
t124 = g(1) * t94 + g(2) * t92;
t86 = g(3) * t99;
t97 = sin(qJ(2));
t107 = t124 * t97 - t86;
t116 = t124 * t99;
t173 = g(3) * t97;
t108 = t116 + t173;
t63 = t91 * t98 + t93 * t96;
t56 = t63 * qJD(2);
t52 = t56 ^ 2;
t133 = qJD(2) * t166;
t156 = qJD(2) * t96;
t53 = t91 * t156 - t133;
t178 = -t53 ^ 2 - t52;
t55 = t63 * qJD(3);
t153 = qJD(3) * t96;
t58 = qJD(3) * t166 - t91 * t153;
t142 = t98 * qJDD(2);
t143 = t96 * qJDD(2);
t122 = -t93 * t142 + t91 * t143;
t30 = qJD(2) * t55 + t122;
t138 = qJD(2) * qJD(3);
t130 = t96 * t138;
t110 = t63 * qJDD(2) - t91 * t130;
t129 = t98 * t138;
t31 = t93 * t129 + t110;
t177 = t30 * pkin(4) - t31 * qJ(5);
t100 = qJD(3) ^ 2;
t139 = qJD(1) * qJD(2);
t141 = t99 * qJDD(1);
t125 = t97 * t139 - t141;
t176 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t100 + (t124 + t139) * t97 - t125 - t86;
t171 = t98 * pkin(3);
t61 = qJDD(2) * pkin(6) + t97 * qJDD(1) + t99 * t139;
t113 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + t61;
t149 = t97 * qJD(1);
t126 = t162 * qJD(2) + t149;
t119 = qJD(3) * t126;
t14 = qJDD(3) * pkin(3) - t113 * t96 - t98 * t119;
t17 = t113 * t98 - t96 * t119;
t5 = t91 * t14 + t93 * t17;
t49 = t126 * t98;
t170 = t91 * t49;
t168 = t92 * t99;
t36 = t93 * t49;
t165 = t94 * t98;
t164 = t94 * t99;
t163 = t96 * t99;
t4 = t93 * t14 - t91 * t17;
t148 = t99 * qJD(1);
t161 = -t93 * t111 - t63 * t148 + t91 * t50;
t48 = t126 * t96;
t38 = qJD(3) * pkin(3) - t48;
t21 = t91 * t38 + t36;
t89 = t96 ^ 2;
t159 = -t98 ^ 2 + t89;
t158 = qJD(2) * pkin(2);
t155 = qJD(2) * t97;
t151 = qJDD(3) * pkin(4);
t83 = pkin(2) + t171;
t59 = -t83 * qJD(2) + qJD(4) - t148;
t150 = t59 * qJD(2);
t23 = -t93 * t48 - t170;
t147 = qJD(5) - t23;
t146 = qJDD(1) - g(3);
t101 = qJD(2) ^ 2;
t145 = t100 + t101;
t144 = qJDD(3) * t96;
t140 = t99 * qJDD(2);
t137 = qJDD(3) * qJ(5) + t5;
t136 = pkin(3) * t153;
t132 = -qJDD(5) + t4;
t131 = t162 * t96;
t11 = t55 * pkin(4) - t58 * qJ(5) - t63 * qJD(5) + t136;
t127 = -t11 + t149;
t123 = g(1) * t92 - g(2) * t94;
t87 = qJ(3) + pkin(8);
t84 = sin(t87);
t85 = cos(t87);
t121 = pkin(4) * t85 + qJ(5) * t84;
t20 = t93 * t38 - t170;
t24 = -t99 * t56 - t58 * t97;
t25 = -qJD(2) * t115 - t97 * t55;
t46 = t63 * t97;
t47 = t62 * t97;
t114 = -t24 * t56 - t25 * t53 + t47 * t30 + t46 * t31;
t73 = -t148 - t158;
t106 = -pkin(6) * qJDD(3) + (t148 + t73 - t158) * qJD(3);
t32 = pkin(3) * t130 - t83 * qJDD(2) + qJDD(4) + t125;
t105 = -t73 * qJD(2) + t108 - t61;
t18 = t53 * pkin(4) - t56 * qJ(5) + t59;
t42 = t84 * t168 + t94 * t85;
t44 = t84 * t164 - t92 * t85;
t104 = g(1) * t44 + g(2) * t42 + t84 * t173 - t18 * t56 + t132;
t103 = t32 - t107;
t68 = t162 * t98;
t33 = t93 * t131 + t91 * t68;
t34 = -t91 * t131 + t93 * t68;
t102 = -t160 * t53 + t161 * t56 - t34 * t30 + t33 * t31 - t108;
t81 = -t93 * pkin(3) - pkin(4);
t79 = t91 * pkin(3) + qJ(5);
t77 = t92 * t171;
t70 = t99 * t83;
t45 = t85 * t164 + t92 * t84;
t43 = t85 * t168 - t94 * t84;
t29 = t62 * pkin(4) - t63 * qJ(5) - t83;
t26 = pkin(3) * t156 + t56 * pkin(4) + t53 * qJ(5);
t22 = -t91 * t48 + t36;
t19 = qJD(3) * qJ(5) + t21;
t16 = -qJD(3) * pkin(4) + qJD(5) - t20;
t3 = -t132 - t151;
t2 = -t56 * qJD(5) + t177 + t32;
t1 = qJD(3) * qJD(5) + t137;
t6 = [t146, 0, -t101 * t97 + t140, -qJDD(2) * t97 - t101 * t99, 0, 0, 0, 0, 0, (-0.2e1 * t130 + t142) * t99 + (-t145 * t98 - t144) * t97, (-qJDD(3) * t97 - 0.2e1 * t99 * t138) * t98 + (t145 * t97 - t140) * t96, t114, t97 * t150 + t20 * t24 + t21 * t25 - t32 * t99 - t4 * t46 - t5 * t47 - g(3), t24 * qJD(3) - t46 * qJDD(3) + t53 * t155 - t99 * t30, t114, t25 * qJD(3) - t47 * qJDD(3) - t155 * t56 + t99 * t31, -t1 * t47 + t155 * t18 - t16 * t24 + t19 * t25 - t2 * t99 + t3 * t46 - g(3); 0, qJDD(2), t141 + t107, -t146 * t97 + t116, t89 * qJDD(2) + 0.2e1 * t96 * t129, -0.2e1 * t159 * t138 + 0.2e1 * t96 * t142, t100 * t98 + t144, qJDD(3) * t98 - t100 * t96, 0, t106 * t96 + t176 * t98, t106 * t98 - t176 * t96, -t20 * t58 - t21 * t55 - t4 * t63 - t5 * t62 + t102, t5 * t34 - t4 * t33 - t32 * t83 - g(3) * (t162 * t97 + t70) + (t136 - t149) * t59 + t160 * t21 - t161 * t20 + t124 * (-t162 * t99 + t83 * t97), -t161 * qJD(3) - t33 * qJDD(3) + t107 * t85 - t127 * t53 + t18 * t55 + t2 * t62 + t29 * t30, -t1 * t62 + t16 * t58 - t19 * t55 + t3 * t63 + t102, qJD(3) * t160 + t34 * qJDD(3) + t107 * t84 + t127 * t56 - t18 * t58 - t2 * t63 - t29 * t31, -g(3) * t70 + t1 * t34 + t18 * t11 + t2 * t29 + t3 * t33 + t160 * t19 + t161 * t16 + (-g(3) * t121 - t124 * t162) * t99 + (-g(3) * t162 - t18 * qJD(1) + t124 * (t121 + t83)) * t97; 0, 0, 0, 0, -t96 * t101 * t98, t159 * t101, t143, t142, qJDD(3), t105 * t96 - t123 * t98, t105 * t98 + t123 * t96, (t21 - t22) * t56 + (-t20 + t23) * t53 + (-t30 * t91 - t31 * t93) * pkin(3), -g(1) * t77 + t20 * t22 - t21 * t23 + (g(2) * t165 + t4 * t93 + t5 * t91 + (t108 - t150) * t96) * pkin(3), t22 * qJD(3) - t26 * t53 + (pkin(4) - t81) * qJDD(3) + t104, -t79 * t30 + t81 * t31 + (t19 - t22) * t56 + (t16 - t147) * t53, -t85 * t173 - g(1) * t45 - g(2) * t43 + t79 * qJDD(3) - t18 * t53 + t26 * t56 + (0.2e1 * qJD(5) - t23) * qJD(3) + t137, t1 * t79 + t3 * t81 - t18 * t26 - t16 * t22 - g(1) * (-t94 * pkin(3) * t163 - t44 * pkin(4) + t45 * qJ(5) + t77) - g(2) * (-t42 * pkin(4) + t43 * qJ(5) + (-t163 * t92 - t165) * pkin(3)) + t147 * t19 - (-pkin(3) * t96 - pkin(4) * t84 + qJ(5) * t85) * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, t20 * t56 + t21 * t53 + t103, 0.2e1 * t56 * qJD(3) + t122, t178, (t53 - t133) * qJD(3) - t110, t19 * t53 + (-qJD(5) - t16) * t56 + t103 + t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t53 - qJDD(3), (t53 + t133) * qJD(3) + t110, -t52 - t100, -t19 * qJD(3) - t104 - t151;];
tau_reg = t6;
