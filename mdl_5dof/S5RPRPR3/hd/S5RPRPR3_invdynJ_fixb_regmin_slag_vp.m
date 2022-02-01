% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR3
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:55
% EndTime: 2022-01-23 09:20:57
% DurationCPUTime: 0.81s
% Computational Cost: add. (1123->174), mult. (1882->249), div. (0->0), fcn. (1143->14), ass. (0->119)
t71 = sin(pkin(8));
t145 = pkin(1) * t71;
t110 = qJD(3) * t145;
t73 = cos(pkin(8));
t58 = t73 * pkin(1) + pkin(2);
t150 = -qJD(1) * t110 + t58 * qJDD(1);
t43 = t58 * qJD(1);
t149 = qJD(3) * t43 + qJDD(1) * t145;
t75 = sin(qJ(3));
t78 = cos(qJ(3));
t122 = t78 * t145 + t75 * t58;
t70 = sin(pkin(9));
t63 = qJDD(1) + qJDD(3);
t118 = t63 * qJ(4);
t147 = -t149 * t78 - t150 * t75;
t66 = qJD(1) + qJD(3);
t12 = t66 * qJD(4) + t118 - t147;
t72 = cos(pkin(9));
t8 = -t72 * qJDD(2) + t70 * t12;
t142 = t8 * t70;
t9 = t70 * qJDD(2) + t72 * t12;
t148 = t9 * t72 + t142;
t111 = qJD(1) * t145;
t27 = -t75 * t111 + t78 * t43;
t94 = qJD(4) - t27;
t131 = t72 * t66;
t45 = -qJD(5) + t131;
t114 = qJD(5) + t45;
t28 = t78 * t111 + t75 * t43;
t20 = t66 * qJ(4) + t28;
t17 = -t72 * qJD(2) + t70 * t20;
t141 = t17 * t70;
t18 = t70 * qJD(2) + t72 * t20;
t146 = -t114 * t18 - t66 * t141;
t67 = qJ(1) + pkin(8);
t61 = qJ(3) + t67;
t56 = sin(t61);
t52 = g(1) * t56;
t57 = cos(t61);
t144 = g(2) * t57;
t143 = t63 * pkin(3);
t140 = t28 * t66;
t125 = t78 * t58;
t90 = qJD(3) * t125 - t75 * t110;
t29 = qJD(4) + t90;
t139 = t29 * t66;
t74 = sin(qJ(5));
t138 = t29 * t74;
t30 = t122 * qJD(3);
t137 = t30 * t66;
t132 = t72 * t63;
t42 = -qJDD(5) + t132;
t136 = t42 * t72;
t135 = t63 * t74;
t62 = t66 ^ 2;
t64 = t70 ^ 2;
t134 = t64 * t62;
t133 = t70 * t74;
t130 = t72 * t74;
t77 = cos(qJ(5));
t129 = t72 * t77;
t128 = t74 * t42;
t127 = t74 * t77;
t126 = t77 * t42;
t124 = t57 * pkin(3) + t56 * qJ(4);
t123 = g(1) * t57 + g(2) * t56;
t121 = t72 ^ 2 + t64;
t69 = t77 ^ 2;
t120 = t74 ^ 2 - t69;
t119 = qJ(4) * t72;
t116 = qJD(5) * t74;
t115 = qJD(5) * t77;
t113 = qJ(4) * qJD(5);
t109 = t66 * t115;
t108 = t45 * t116;
t106 = -t56 * pkin(3) + t57 * qJ(4);
t97 = -t149 * t75 + t150 * t78;
t91 = qJDD(4) - t97;
t13 = t91 - t143;
t105 = -t13 - t144;
t104 = t121 * t63;
t102 = t42 - t132;
t101 = t42 + t132;
t100 = -t75 * t145 + t125;
t99 = t66 * t114;
t98 = t94 * t77;
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t95 = g(1) * t76 - g(2) * t79;
t93 = t18 * t72 + t141;
t92 = -t123 + t148;
t37 = -t72 * pkin(4) - t70 * pkin(7) - pkin(3);
t89 = -t74 * t119 + t77 * t37;
t14 = t37 * t66 + t94;
t23 = t56 * t130 + t57 * t77;
t25 = -t57 * t130 + t56 * t77;
t7 = t37 * t63 + t91;
t88 = -g(1) * t23 - g(2) * t25 + (t74 * t7 + t77 * t9 + (t77 * t14 - t74 * t18) * qJD(5)) * t72 + t77 * t142;
t24 = -t56 * t129 + t57 * t74;
t26 = t57 * t129 + t56 * t74;
t87 = -g(1) * t24 - g(2) * t26 + t115 * t141 + t8 * t133;
t85 = -t52 - t97 + t144;
t84 = g(3) * t70 - t114 * t14 - t9;
t83 = t77 * t113 + t94 * t74;
t82 = -t45 ^ 2 - t134;
t81 = t123 + t147;
t38 = t72 * t52;
t34 = t70 * t116 * t131;
t32 = -pkin(3) - t100;
t31 = qJ(4) + t122;
t22 = (-0.2e1 * t74 * t109 + t63 * t69) * t64;
t21 = -t100 + t37;
t19 = -t66 * pkin(3) + t94;
t16 = 0.2e1 * (t120 * t66 * qJD(5) - t63 * t127) * t64;
t11 = (t101 * t74 + (t45 + t131) * t115) * t70;
t10 = t34 + (-t101 * t77 + t108) * t70;
t5 = t77 * t7;
t2 = -t74 * t9 + t5 + (-t74 * t14 - t77 * t18) * qJD(5);
t1 = [qJDD(1), t95, g(1) * t79 + g(2) * t76, (t95 + (t71 ^ 2 + t73 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t63, t100 * t63 - t137 - t85, -t122 * t63 - t90 * t66 + t81, t38 + (-t32 * t63 + t105 - t137) * t72, t31 * t104 + t121 * t139 + t92, t13 * t32 + t19 * t30 - g(1) * (-pkin(2) * sin(t67) - t76 * pkin(1) + t106) - g(2) * (pkin(2) * cos(t67) + t79 * pkin(1) + t124) + t148 * t31 + t93 * t29, t22, t16, t10, t11, t136, -(-t21 * t116 + t77 * t30) * t45 - t21 * t126 + (-(-t31 * t115 - t138) * t45 + t31 * t128 - t2) * t72 + (t66 * t138 + (t109 + t135) * t31) * t64 + t87, (t29 * t129 + t74 * t30) * t45 + (t31 * t129 + t74 * t21) * t42 + (t31 * t63 + t139) * t77 * t64 + (t77 * t21 * t45 + (-t141 + (-t45 * t72 - t64 * t66) * t31) * t74) * qJD(5) + t88; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, t9 * t70 - t8 * t72 - g(3), 0, 0, 0, 0, 0, (t102 * t74 + (t45 - t131) * t115) * t70, t34 + (t102 * t77 - t108) * t70; 0, 0, 0, 0, t63, -t85 + t140, t27 * t66 + t81, t38 + (t105 + t140 + t143) * t72, t94 * t66 * t121 + qJ(4) * t104 + t92, -t13 * pkin(3) - t19 * t28 - g(1) * t106 - g(2) * t124 + (t9 * qJ(4) + t94 * t18) * t72 + (t8 * qJ(4) + t94 * t17) * t70, t22, t16, t10, t11, t136, -t89 * t42 - t2 * t72 + (t37 * t116 + t77 * t28 + t83 * t72) * t45 + (t74 * t118 + t83 * t66) * t64 + t87, (t77 * t119 + t74 * t37) * t42 + (-t74 * t28 + t72 * t98) * t45 + (-t17 * t133 + t89 * t45) * qJD(5) + (t77 * t118 + (-t74 * t113 + t98) * t66) * t64 + t88; 0, 0, 0, 0, 0, 0, 0, -t132, -t121 * t62, -t93 * t66 - t105 - t52, 0, 0, 0, 0, 0, t82 * t74 - t126, t82 * t77 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127 * t134, -t120 * t134, (t63 * t77 - t74 * t99) * t70, (-t77 * t99 - t135) * t70, -t42, -g(1) * t25 + g(2) * t23 + t146 * t77 + t84 * t74 + t5, g(1) * t26 - g(2) * t24 + t84 * t77 + (-t146 - t7) * t74;];
tau_reg = t1;
