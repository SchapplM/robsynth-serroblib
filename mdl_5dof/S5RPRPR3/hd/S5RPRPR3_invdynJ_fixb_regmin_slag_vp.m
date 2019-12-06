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
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:51:35
% EndTime: 2019-12-05 17:51:39
% DurationCPUTime: 0.85s
% Computational Cost: add. (1167->178), mult. (1952->253), div. (0->0), fcn. (1186->14), ass. (0->120)
t69 = qJ(1) + pkin(8);
t63 = qJ(3) + t69;
t58 = sin(t63);
t59 = cos(t63);
t153 = g(2) * t59 + g(3) * t58;
t73 = sin(pkin(8));
t149 = pkin(1) * t73;
t115 = qJD(3) * t149;
t75 = cos(pkin(8));
t60 = t75 * pkin(1) + pkin(2);
t155 = -qJD(1) * t115 + t60 * qJDD(1);
t45 = t60 * qJD(1);
t154 = qJD(3) * t45 + qJDD(1) * t149;
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t127 = t80 * t149 + t77 * t60;
t72 = sin(pkin(9));
t65 = qJDD(1) + qJDD(3);
t123 = t65 * qJ(4);
t151 = -t80 * t154 - t155 * t77;
t68 = qJD(1) + qJD(3);
t12 = t68 * qJD(4) + t123 - t151;
t74 = cos(pkin(9));
t8 = -t74 * qJDD(2) + t72 * t12;
t147 = t8 * t72;
t9 = t72 * qJDD(2) + t74 * t12;
t152 = t9 * t74 + t147;
t116 = qJD(1) * t149;
t28 = -t77 * t116 + t80 * t45;
t99 = qJD(4) - t28;
t136 = t74 * t68;
t48 = -qJD(5) + t136;
t119 = qJD(5) + t48;
t29 = t80 * t116 + t77 * t45;
t21 = t68 * qJ(4) + t29;
t18 = -t74 * qJD(2) + t72 * t21;
t146 = t18 * t72;
t19 = t72 * qJD(2) + t74 * t21;
t150 = -t119 * t19 - t68 * t146;
t148 = t65 * pkin(3);
t145 = t29 * t68;
t130 = t80 * t60;
t92 = qJD(3) * t130 - t77 * t115;
t30 = qJD(4) + t92;
t144 = t30 * t68;
t76 = sin(qJ(5));
t143 = t30 * t76;
t31 = t127 * qJD(3);
t142 = t31 * t68;
t137 = t74 * t65;
t44 = -qJDD(5) + t137;
t141 = t44 * t74;
t140 = t65 * t76;
t64 = t68 ^ 2;
t66 = t72 ^ 2;
t139 = t66 * t64;
t138 = t72 * t76;
t135 = t74 * t76;
t79 = cos(qJ(5));
t134 = t74 * t79;
t133 = t76 * t44;
t132 = t76 * t79;
t131 = t79 * t44;
t129 = t153 * t74;
t128 = g(2) * t58 - g(3) * t59;
t126 = t74 ^ 2 + t66;
t71 = t79 ^ 2;
t125 = t76 ^ 2 - t71;
t124 = qJ(4) * t74;
t121 = qJD(5) * t76;
t120 = qJD(5) * t79;
t118 = qJ(4) * qJD(5);
t114 = t68 * t120;
t113 = t48 * t121;
t111 = -t58 * pkin(3) + t59 * qJ(4);
t110 = t126 * t65;
t108 = t44 - t137;
t107 = t44 + t137;
t106 = -t77 * t149 + t130;
t105 = t68 * t119;
t104 = t99 * t79;
t103 = -t154 * t77 + t155 * t80;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t100 = g(2) * t81 + g(3) * t78;
t98 = t145 + t148;
t97 = -t59 * pkin(3) - t58 * qJ(4);
t96 = t19 * t74 + t146;
t33 = -pkin(3) - t106;
t95 = -t33 * t65 - t142;
t94 = t128 + t152;
t93 = qJDD(4) - t103;
t38 = -t74 * pkin(4) - t72 * pkin(7) - pkin(3);
t91 = t103 + t153;
t90 = -t76 * t124 + t79 * t38;
t15 = t38 * t68 + t99;
t24 = t58 * t135 + t59 * t79;
t26 = t59 * t135 - t58 * t79;
t7 = t38 * t65 + t93;
t89 = -g(2) * t26 - g(3) * t24 + (t76 * t7 + t79 * t9 + (t79 * t15 - t76 * t19) * qJD(5)) * t74 + t79 * t147;
t25 = t58 * t134 - t59 * t76;
t27 = -t59 * t134 - t58 * t76;
t88 = -g(2) * t27 + g(3) * t25 + t120 * t146 + t8 * t138;
t86 = g(1) * t72 - t119 * t15 - t9;
t14 = t93 - t148;
t85 = t79 * t118 + t99 * t76;
t84 = -t48 ^ 2 - t139;
t83 = -t128 + t151;
t35 = t72 * t121 * t136;
t32 = qJ(4) + t127;
t23 = (-0.2e1 * t76 * t114 + t65 * t71) * t66;
t22 = -t106 + t38;
t20 = -t68 * pkin(3) + t99;
t17 = 0.2e1 * (t125 * t68 * qJD(5) - t65 * t132) * t66;
t13 = t14 * t72;
t11 = (t107 * t76 + (t48 + t136) * t120) * t72;
t10 = t35 + (-t107 * t79 + t113) * t72;
t5 = t79 * t7;
t2 = -t76 * t9 + t5 + (-t76 * t15 - t79 * t19) * qJD(5);
t1 = [qJDD(1), t100, -g(2) * t78 + g(3) * t81, (t100 + (t73 ^ 2 + t75 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t65, t106 * t65 - t142 + t91, -t127 * t65 - t92 * t68 + t83, (-t14 + t95) * t74 + t129, t13 + (-t153 - t95) * t72, t32 * t110 + t126 * t144 + t94, t14 * t33 + t20 * t31 - g(2) * (-pkin(2) * cos(t69) - t81 * pkin(1) + t97) - g(3) * (-pkin(2) * sin(t69) - t78 * pkin(1) + t111) + t152 * t32 + t96 * t30, t23, t17, t10, t11, t141, -(-t22 * t121 + t79 * t31) * t48 - t22 * t131 + (-(-t32 * t120 - t143) * t48 + t32 * t133 - t2) * t74 + (t68 * t143 + (t114 + t140) * t32) * t66 + t88, (t30 * t134 + t76 * t31) * t48 + (t32 * t134 + t76 * t22) * t44 + (t32 * t65 + t144) * t79 * t66 + (t79 * t22 * t48 + (-t146 + (-t48 * t74 - t66 * t68) * t32) * t76) * qJD(5) + t89; 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, t9 * t72 - t8 * t74 - g(1), 0, 0, 0, 0, 0, (t108 * t76 + (t48 - t136) * t120) * t72, t35 + (t108 * t79 - t113) * t72; 0, 0, 0, 0, t65, t91 + t145, t28 * t68 + t83, (-t14 + t98) * t74 + t129, t13 + (-t153 - t98) * t72, t99 * t68 * t126 + qJ(4) * t110 + t94, -t14 * pkin(3) - t20 * t29 - g(2) * t97 - g(3) * t111 + (t9 * qJ(4) + t99 * t19) * t74 + (t8 * qJ(4) + t99 * t18) * t72, t23, t17, t10, t11, t141, -t90 * t44 - t2 * t74 + (t38 * t121 + t79 * t29 + t74 * t85) * t48 + (t76 * t123 + t68 * t85) * t66 + t88, (t79 * t124 + t76 * t38) * t44 + (t104 * t74 - t76 * t29) * t48 + (-t18 * t138 + t48 * t90) * qJD(5) + (t79 * t123 + (-t76 * t118 + t104) * t68) * t66 + t89; 0, 0, 0, 0, 0, 0, 0, -t137, t72 * t65, -t126 * t64, -t68 * t96 + t14 - t153, 0, 0, 0, 0, 0, t76 * t84 - t131, t84 * t79 + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132 * t139, -t125 * t139, (-t105 * t76 + t65 * t79) * t72, (-t105 * t79 - t140) * t72, -t44, -g(2) * t24 + g(3) * t26 + t150 * t79 + t86 * t76 + t5, -g(2) * t25 - g(3) * t27 + t86 * t79 + (-t150 - t7) * t76;];
tau_reg = t1;
