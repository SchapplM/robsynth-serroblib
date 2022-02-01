% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP2
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
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:12
% EndTime: 2022-01-23 09:28:15
% DurationCPUTime: 0.86s
% Computational Cost: add. (1303->208), mult. (2181->241), div. (0->0), fcn. (1229->12), ass. (0->135)
t82 = sin(pkin(8));
t160 = pkin(1) * t82;
t126 = qJD(1) * t160;
t83 = cos(pkin(8));
t68 = t83 * pkin(1) + pkin(2);
t165 = -qJD(3) * t126 + t68 * qJDD(1);
t164 = qJDD(2) - g(3);
t53 = t68 * qJD(1);
t163 = qJD(3) * t53 + qJDD(1) * t160;
t79 = qJ(1) + pkin(8);
t71 = qJ(3) + t79;
t66 = sin(t71);
t60 = g(2) * t66;
t67 = cos(t71);
t62 = g(1) * t67;
t142 = -t60 - t62;
t86 = sin(qJ(3));
t89 = cos(qJ(3));
t141 = t89 * t160 + t86 * t68;
t78 = qJD(1) + qJD(3);
t138 = qJ(5) * t78;
t31 = t89 * t126 + t86 * t53;
t23 = t78 * pkin(7) + t31;
t117 = t23 + t138;
t88 = cos(qJ(4));
t104 = t117 * t88;
t85 = sin(qJ(4));
t132 = t85 * qJD(4);
t123 = t78 * t132;
t155 = t88 * pkin(4);
t69 = pkin(3) + t155;
t77 = qJDD(1) + qJDD(3);
t149 = t69 * t77;
t162 = -pkin(4) * t123 + t149;
t118 = -t86 * t160 + t89 * t68;
t161 = -t163 * t89 - t165 * t86;
t80 = t85 ^ 2;
t159 = pkin(4) * t80;
t61 = g(1) * t66;
t158 = g(2) * t67;
t157 = g(3) * t88;
t156 = t77 * pkin(3);
t30 = -t86 * t126 + t89 * t53;
t22 = -t78 * pkin(3) - t30;
t154 = t22 * t78;
t153 = t31 * t78;
t33 = t141 * qJD(3);
t152 = t33 * t78;
t151 = t66 * t88;
t150 = t67 * t85;
t76 = t78 ^ 2;
t148 = t76 * t88;
t147 = t78 * t85;
t146 = t78 * t88;
t65 = t85 * t77;
t145 = t88 * t77;
t84 = -qJ(5) - pkin(7);
t73 = t88 * qJD(2);
t12 = -t117 * t85 + t73;
t137 = qJD(4) * pkin(4);
t11 = t12 + t137;
t143 = t11 - t12;
t81 = t88 ^ 2;
t140 = -t80 - t81;
t139 = t80 - t81;
t36 = pkin(7) + t141;
t136 = -qJ(5) - t36;
t134 = qJD(4) * t78;
t133 = qJDD(4) * pkin(4);
t131 = t88 * qJD(4);
t19 = t23 * t132;
t9 = t77 * pkin(7) - t161;
t112 = -qJD(4) * qJD(2) - t9;
t100 = -qJ(5) * t77 + t112;
t95 = qJD(5) * t78 - t100;
t3 = -t19 + (-qJ(5) * t134 + qJDD(2)) * t85 + t95 * t88;
t130 = t3 * t88 + t142;
t16 = -t69 * t78 + qJD(5) - t30;
t46 = g(2) * t150;
t114 = -t163 * t86 + t165 * t89;
t5 = qJDD(5) - t114 - t162;
t129 = t16 * t131 + t5 * t85 + t46;
t10 = -t114 - t156;
t128 = t10 * t85 + t22 * t131 + t46;
t47 = g(1) * t151;
t127 = t30 * t132 + t31 * t146 + t47;
t125 = pkin(4) * t132;
t122 = -t5 - t158;
t121 = -t10 - t158;
t120 = -t66 * t84 + t67 * t69;
t116 = qJD(4) * t84;
t115 = 0.2e1 * t78 * t131;
t113 = qJD(4) * t136;
t35 = -pkin(3) - t118;
t91 = qJD(4) ^ 2;
t110 = -pkin(7) * t91 + t156;
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t109 = g(1) * t87 - g(2) * t90;
t108 = -t153 - t61;
t13 = t85 * qJD(2) + t104;
t107 = t11 * t85 - t13 * t88;
t28 = t33 + t125;
t29 = t35 - t155;
t106 = t28 * t78 + t29 * t77;
t105 = -t66 * t69 - t67 * t84;
t70 = t88 * qJDD(2);
t103 = g(1) * t150 + t85 * t60 - t157 + t70;
t102 = -pkin(3) * t134 - pkin(7) * qJDD(4);
t99 = g(2) * t151 - t164 * t85 + t88 * t62 + t19;
t98 = -t114 - t61 + t158;
t97 = t35 * t77 + t36 * t91 + t152;
t32 = t118 * qJD(3);
t96 = -qJDD(4) * t36 + (t35 * t78 - t32) * qJD(4);
t94 = (-qJD(5) - t16) * t78 + t100;
t93 = -t142 + t161;
t74 = t88 * qJ(5);
t72 = t88 * qJD(5);
t55 = t88 * pkin(7) + t74;
t54 = t84 * t85;
t43 = qJDD(4) * t88 - t91 * t85;
t42 = qJDD(4) * t85 + t91 * t88;
t38 = -t85 * qJD(5) + t88 * t116;
t37 = t85 * t116 + t72;
t34 = t85 * t115 + t80 * t77;
t27 = t88 * t36 + t74;
t26 = t136 * t85;
t25 = t30 * t131;
t20 = -0.2e1 * t139 * t134 + 0.2e1 * t85 * t145;
t17 = t22 * t132;
t14 = t16 * t132;
t7 = (-qJD(5) - t32) * t85 + t88 * t113;
t6 = t85 * t113 + t88 * t32 + t72;
t2 = -qJD(4) * t104 - t95 * t85 + t133 + t70;
t1 = [qJDD(1), t109, g(1) * t90 + g(2) * t87, (t109 + (t82 ^ 2 + t83 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t77, t118 * t77 - t152 - t98, -t141 * t77 - t32 * t78 + t93, t34, t20, t42, t43, 0, t17 + t47 + t96 * t85 + (t121 - t97) * t88, t96 * t88 + (t97 - t61) * t85 + t128, t26 * qJDD(4) + t14 + t47 + (t29 * t147 + t7) * qJD(4) + (-t106 + t122) * t88, -t27 * qJDD(4) + (t29 * t146 - t6) * qJD(4) + (t106 - t61) * t85 + t129, (t27 * t77 + t6 * t78 + (-t26 * t78 - t11) * qJD(4)) * t88 + (-t26 * t77 - t7 * t78 - t2 + (-t27 * t78 - t13) * qJD(4)) * t85 + t130, t3 * t27 + t13 * t6 + t2 * t26 + t11 * t7 + t5 * t29 + t16 * t28 - g(1) * (-pkin(2) * sin(t79) - t87 * pkin(1) + t105) - g(2) * (pkin(2) * cos(t79) + t90 * pkin(1) + t120); 0, 0, 0, t164, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t42, t43, -t42, 0, -t107 * qJD(4) + t2 * t88 + t3 * t85 - g(3); 0, 0, 0, 0, t77, -t98 + t153, t30 * t78 + t93, t34, t20, t42, t43, 0, t17 + t102 * t85 + (t110 + t121) * t88 + t127, t25 + t102 * t88 + (t108 - t110) * t85 + t128, t54 * qJDD(4) + t14 + (-t69 * t147 + t38) * qJD(4) + (t122 + t162) * t88 + t127, -t55 * qJDD(4) + t25 + (t108 - t149) * t85 + (-t37 + (-t69 * t88 + t159) * t78) * qJD(4) + t129, (-qJD(4) * t11 + t55 * t77) * t88 + (-t13 * qJD(4) - t54 * t77 - t2) * t85 + (t37 * t88 - t38 * t85 + t140 * t30 + (-t54 * t88 - t55 * t85) * qJD(4)) * t78 + t130, t3 * t55 + t2 * t54 - t5 * t69 - g(1) * t105 - g(2) * t120 + (-t31 + t125) * t16 + (-t88 * t30 + t37) * t13 + (t85 * t30 + t38) * t11; 0, 0, 0, 0, 0, 0, 0, -t85 * t148, t139 * t76, t65, t145, qJDD(4), (-t9 - t154) * t85 + t103, (-t85 * t23 + t73) * qJD(4) + (t112 - t154) * t88 + t99, 0.2e1 * t133 + (t13 - t104) * qJD(4) + (pkin(4) * t148 + t94) * t85 + t103, -t76 * t159 + (t85 * t138 + t12) * qJD(4) + t94 * t88 + t99, -pkin(4) * t65 + (-t137 + t143) * t146, t143 * t13 + (-t157 + t2 + (-t16 * t78 - t142) * t85) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t123 - t145, t65 + t115, t140 * t76, t107 * t78 - t122 - t61;];
tau_reg = t1;
