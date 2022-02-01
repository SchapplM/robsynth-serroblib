% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:07
% EndTime: 2022-01-23 09:19:08
% DurationCPUTime: 0.63s
% Computational Cost: add. (1006->150), mult. (1667->192), div. (0->0), fcn. (1095->16), ass. (0->107)
t93 = sin(pkin(8));
t150 = pkin(1) * t93;
t128 = qJD(3) * t150;
t95 = cos(pkin(8));
t75 = t95 * pkin(1) + pkin(2);
t157 = -qJD(1) * t128 + t75 * qJDD(1);
t91 = qJ(1) + pkin(8);
t83 = qJ(3) + t91;
t74 = cos(t83);
t149 = g(2) * t74;
t73 = sin(t83);
t68 = g(1) * t73;
t156 = t149 - t68;
t60 = t75 * qJD(1);
t155 = qJD(3) * t60 + qJDD(1) * t150;
t94 = cos(pkin(9));
t99 = cos(qJ(5));
t139 = t99 * t94;
t92 = sin(pkin(9));
t96 = sin(qJ(5));
t142 = t96 * t92;
t43 = -t139 + t142;
t44 = t99 * t92 + t96 * t94;
t90 = qJD(1) + qJD(3);
t36 = t44 * t90;
t154 = g(1) * t74 + g(2) * t73;
t100 = cos(qJ(3));
t97 = sin(qJ(3));
t151 = -t155 * t100 - t157 * t97;
t86 = qJDD(1) + qJDD(3);
t12 = t86 * qJ(4) + t90 * qJD(4) - t151;
t78 = t94 * qJDD(2);
t8 = -t92 * t12 + t78;
t9 = t92 * qJDD(2) + t94 * t12;
t153 = -t8 * t92 + t9 * t94;
t136 = t100 * t150 + t97 * t75;
t135 = t92 ^ 2 + t94 ^ 2;
t152 = t135 * t90;
t129 = qJD(1) * t150;
t31 = t100 * t60 - t97 * t129;
t116 = qJD(4) - t31;
t147 = t86 * pkin(3);
t146 = t94 * pkin(4);
t32 = t100 * t129 + t97 * t60;
t145 = t32 * t90;
t37 = t136 * qJD(3);
t144 = t37 * t90;
t143 = t94 * t86;
t138 = t74 * pkin(3) + t73 * qJ(4);
t134 = t100 * t75;
t132 = t90 * t142;
t131 = t90 * t139;
t130 = qJD(5) * t131 + t44 * t86;
t76 = -pkin(3) - t146;
t126 = -t73 * pkin(3) + t74 * qJ(4);
t121 = t157 * t100 - t155 * t97;
t109 = qJDD(4) - t121;
t13 = t109 - t147;
t125 = -t13 - t149;
t124 = t135 * t86;
t122 = -t97 * t150 + t134;
t119 = t43 * t86;
t39 = -pkin(3) - t122;
t101 = cos(qJ(1));
t98 = sin(qJ(1));
t117 = g(1) * t98 - g(2) * t101;
t27 = t90 * qJ(4) + t32;
t20 = t94 * qJD(2) - t92 * t27;
t21 = t92 * qJD(2) + t94 * t27;
t115 = t20 * t92 - t21 * t94;
t38 = qJ(4) + t136;
t28 = (-pkin(7) - t38) * t92;
t84 = t94 * pkin(7);
t29 = t94 * t38 + t84;
t114 = t99 * t28 - t96 * t29;
t113 = t96 * t28 + t99 * t29;
t51 = (-pkin(7) - qJ(4)) * t92;
t52 = t94 * qJ(4) + t84;
t112 = t99 * t51 - t96 * t52;
t111 = t96 * t51 + t99 * t52;
t110 = -t154 + t153;
t108 = qJD(3) * t134 - t97 * t128;
t10 = t76 * t86 + t109;
t22 = t76 * t90 + t116;
t40 = t43 * qJD(5);
t89 = pkin(9) + qJ(5);
t81 = sin(t89);
t107 = t10 * t44 + t156 * t81 - t22 * t40;
t41 = t44 * qJD(5);
t82 = cos(t89);
t106 = t10 * t43 - t156 * t82 + t22 * t41;
t104 = -t121 + t156;
t103 = t154 + t151;
t50 = t94 * t68;
t34 = -t131 + t132;
t33 = qJD(4) + t108;
t30 = t39 - t146;
t26 = -t90 * pkin(3) + t116;
t24 = -t41 * qJD(5) - t43 * qJDD(5);
t23 = -t40 * qJD(5) + t44 * qJDD(5);
t19 = t90 * t41 + t119;
t18 = -qJD(5) * t132 + t130;
t4 = pkin(7) * t143 + t9;
t3 = t78 + (-pkin(7) * t86 - t12) * t92;
t2 = t18 * t44 - t36 * t40;
t1 = -t18 * t43 - t44 * t19 + t40 * t34 - t36 * t41;
t5 = [qJDD(1), t117, g(1) * t101 + g(2) * t98, (t117 + (t93 ^ 2 + t95 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t86, t122 * t86 - t104 - t144, -t108 * t90 - t136 * t86 + t103, t50 + (-t39 * t86 + t125 - t144) * t94, t38 * t124 + t33 * t152 + t110, t13 * t39 + t26 * t37 - g(1) * (-pkin(2) * sin(t91) - t98 * pkin(1) + t126) - g(2) * (pkin(2) * cos(t91) + t101 * pkin(1) + t138) + t153 * t38 - t115 * t33, t2, t1, t23, t24, 0, t37 * t34 + t30 * t19 + t114 * qJDD(5) + (-t113 * qJD(5) - t44 * t33) * qJD(5) + t106, t37 * t36 + t30 * t18 - t113 * qJDD(5) + (-t114 * qJD(5) + t43 * t33) * qJD(5) + t107; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, t8 * t94 + t9 * t92 - g(3), 0, 0, 0, 0, 0, t24, -t23; 0, 0, 0, 0, t86, -t104 + t145, t31 * t90 + t103, t50 + (t125 + t145 + t147) * t94, qJ(4) * t124 + t116 * t152 + t110, -t13 * pkin(3) - t26 * t32 - g(1) * t126 - g(2) * t138 + (t9 * qJ(4) + t116 * t21) * t94 + (-t8 * qJ(4) - t116 * t20) * t92, t2, t1, t23, t24, 0, t76 * t19 + t112 * qJDD(5) - t32 * t34 + (-t111 * qJD(5) - t116 * t44) * qJD(5) + t106, t76 * t18 - t111 * qJDD(5) - t32 * t36 + (-t112 * qJD(5) + t116 * t43) * qJD(5) + t107; 0, 0, 0, 0, 0, 0, 0, -t143, -t135 * t90 ^ 2, t115 * t90 - t125 - t68, 0, 0, 0, 0, 0, 0.2e1 * t36 * qJD(5) + t119, (-t34 - t132) * qJD(5) + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t34, -t34 ^ 2 + t36 ^ 2, (t34 - t132) * qJD(5) + t130, -t119, qJDD(5), -g(3) * t82 + t154 * t81 - t22 * t36 + t99 * t3 - t96 * t4, g(3) * t81 + t154 * t82 + t22 * t34 - t96 * t3 - t99 * t4;];
tau_reg = t5;
