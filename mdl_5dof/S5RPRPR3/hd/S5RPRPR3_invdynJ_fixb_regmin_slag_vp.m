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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:36:26
% EndTime: 2020-01-03 11:36:31
% DurationCPUTime: 0.88s
% Computational Cost: add. (1167->176), mult. (1952->252), div. (0->0), fcn. (1186->14), ass. (0->120)
t69 = qJ(1) + pkin(8);
t63 = qJ(3) + t69;
t58 = sin(t63);
t59 = cos(t63);
t101 = g(2) * t59 + g(3) * t58;
t65 = qJDD(1) + qJDD(3);
t149 = t65 * pkin(3);
t73 = sin(pkin(8));
t152 = pkin(1) * t73;
t75 = cos(pkin(8));
t60 = t75 * pkin(1) + pkin(2);
t45 = t60 * qJD(1);
t156 = qJD(3) * t45 + qJDD(1) * t152;
t115 = qJD(3) * t152;
t157 = -qJD(1) * t115 + t60 * qJDD(1);
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t103 = -t156 * t77 + t157 * t80;
t94 = qJDD(4) - t103;
t14 = t94 - t149;
t92 = -t101 - t14;
t128 = t80 * t152 + t77 * t60;
t72 = sin(pkin(9));
t124 = t65 * qJ(4);
t154 = -t80 * t156 - t157 * t77;
t68 = qJD(1) + qJD(3);
t12 = t68 * qJD(4) + t124 - t154;
t74 = cos(pkin(9));
t8 = -t74 * qJDD(2) + t72 * t12;
t148 = t8 * t72;
t9 = t72 * qJDD(2) + t74 * t12;
t155 = t9 * t74 + t148;
t116 = qJD(1) * t152;
t28 = -t77 * t116 + t80 * t45;
t99 = qJD(4) - t28;
t137 = t74 * t68;
t48 = -qJD(5) + t137;
t120 = qJD(5) + t48;
t29 = t80 * t116 + t77 * t45;
t21 = t68 * qJ(4) + t29;
t18 = -t74 * qJD(2) + t72 * t21;
t147 = t18 * t72;
t19 = t72 * qJD(2) + t74 * t21;
t153 = -t120 * t19 - t68 * t147;
t146 = t29 * t68;
t131 = t80 * t60;
t93 = qJD(3) * t131 - t77 * t115;
t30 = qJD(4) + t93;
t145 = t30 * t68;
t76 = sin(qJ(5));
t144 = t30 * t76;
t31 = t128 * qJD(3);
t143 = t31 * t68;
t138 = t74 * t65;
t44 = -qJDD(5) + t138;
t142 = t44 * t74;
t141 = t65 * t76;
t64 = t68 ^ 2;
t66 = t72 ^ 2;
t140 = t66 * t64;
t139 = t72 * t76;
t136 = t74 * t76;
t79 = cos(qJ(5));
t135 = t74 * t79;
t134 = t76 * t44;
t133 = t76 * t79;
t132 = t79 * t44;
t130 = t59 * pkin(3) + t58 * qJ(4);
t129 = g(2) * t58 - g(3) * t59;
t127 = t74 ^ 2 + t66;
t71 = t79 ^ 2;
t126 = t76 ^ 2 - t71;
t125 = qJ(4) * t74;
t122 = qJD(5) * t76;
t121 = qJD(5) * t79;
t119 = qJ(4) * qJD(5);
t117 = t92 * t72;
t114 = t68 * t121;
t113 = t48 * t122;
t111 = t127 * t65;
t109 = t44 - t138;
t108 = t44 + t138;
t107 = -t77 * t152 + t131;
t106 = t58 * pkin(3) - t59 * qJ(4);
t105 = t68 * t120;
t104 = t99 * t79;
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t100 = -g(2) * t81 - g(3) * t78;
t98 = -t146 - t149;
t97 = t19 * t74 + t147;
t33 = -pkin(3) - t107;
t96 = t33 * t65 + t143;
t95 = -t129 + t155;
t38 = -t74 * pkin(4) - t72 * pkin(7) - pkin(3);
t91 = -t76 * t125 + t79 * t38;
t15 = t38 * t68 + t99;
t24 = -t58 * t136 - t59 * t79;
t26 = t59 * t136 - t58 * t79;
t7 = t38 * t65 + t94;
t90 = g(2) * t26 - g(3) * t24 + (t76 * t7 + t79 * t9 + (t79 * t15 - t76 * t19) * qJD(5)) * t74 + t79 * t148;
t25 = t58 * t135 - t59 * t76;
t27 = t59 * t135 + t58 * t76;
t89 = -g(2) * t27 - g(3) * t25 + t121 * t147 + t8 * t139;
t87 = g(1) * t72 - t120 * t15 - t9;
t86 = t79 * t119 + t99 * t76;
t85 = -t48 ^ 2 - t140;
t84 = t101 - t103;
t83 = t129 + t154;
t35 = t72 * t122 * t137;
t32 = qJ(4) + t128;
t23 = (-0.2e1 * t76 * t114 + t65 * t71) * t66;
t22 = -t107 + t38;
t20 = -t68 * pkin(3) + t99;
t17 = 0.2e1 * (t126 * t68 * qJD(5) - t65 * t133) * t66;
t11 = (t108 * t76 + (t48 + t137) * t121) * t72;
t10 = t35 + (-t108 * t79 + t113) * t72;
t5 = t79 * t7;
t2 = -t76 * t9 + t5 + (-t76 * t15 - t79 * t19) * qJD(5);
t1 = [qJDD(1), t100, g(2) * t78 - g(3) * t81, (t100 + (t73 ^ 2 + t75 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t65, t107 * t65 - t143 - t84, -t128 * t65 - t93 * t68 + t83, (t92 - t96) * t74, t72 * t96 - t117, t32 * t111 + t127 * t145 + t95, t14 * t33 + t20 * t31 - g(2) * (pkin(2) * cos(t69) + t81 * pkin(1) + t130) - g(3) * (pkin(2) * sin(t69) + t78 * pkin(1) + t106) + t155 * t32 + t97 * t30, t23, t17, t10, t11, t142, -(-t22 * t122 + t79 * t31) * t48 - t22 * t132 + (-(-t32 * t121 - t144) * t48 + t32 * t134 - t2) * t74 + (t68 * t144 + (t114 + t141) * t32) * t66 + t89, (t30 * t135 + t76 * t31) * t48 + (t32 * t135 + t76 * t22) * t44 + (t32 * t65 + t145) * t79 * t66 + (t79 * t22 * t48 + (-t147 + (-t48 * t74 - t66 * t68) * t32) * t76) * qJD(5) + t90; 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, t9 * t72 - t8 * t74 - g(1), 0, 0, 0, 0, 0, (t109 * t76 + (t48 - t137) * t121) * t72, t35 + (t109 * t79 - t113) * t72; 0, 0, 0, 0, t65, -t84 + t146, t28 * t68 + t83, (t92 - t98) * t74, t72 * t98 - t117, t99 * t68 * t127 + qJ(4) * t111 + t95, -t14 * pkin(3) - t20 * t29 - g(2) * t130 - g(3) * t106 + (t9 * qJ(4) + t99 * t19) * t74 + (t8 * qJ(4) + t99 * t18) * t72, t23, t17, t10, t11, t142, -t91 * t44 - t2 * t74 + (t38 * t122 + t79 * t29 + t74 * t86) * t48 + (t76 * t124 + t68 * t86) * t66 + t89, (t79 * t125 + t76 * t38) * t44 + (t104 * t74 - t76 * t29) * t48 + (-t18 * t139 + t48 * t91) * qJD(5) + (t79 * t124 + (-t76 * t119 + t104) * t68) * t66 + t90; 0, 0, 0, 0, 0, 0, 0, -t138, t72 * t65, -t127 * t64, -t68 * t97 + qJDD(4) - t149 + t84, 0, 0, 0, 0, 0, t76 * t85 - t132, t79 * t85 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 * t140, -t126 * t140, (-t105 * t76 + t65 * t79) * t72, (-t105 * t79 - t141) * t72, -t44, -g(2) * t24 - g(3) * t26 + t153 * t79 + t87 * t76 + t5, g(2) * t25 - g(3) * t27 + t87 * t79 + (-t153 - t7) * t76;];
tau_reg = t1;
