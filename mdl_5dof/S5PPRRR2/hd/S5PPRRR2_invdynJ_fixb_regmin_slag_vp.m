% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRR2
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:52
% EndTime: 2019-12-05 15:14:56
% DurationCPUTime: 0.92s
% Computational Cost: add. (696->150), mult. (1518->215), div. (0->0), fcn. (1225->14), ass. (0->108)
t71 = sin(pkin(9));
t73 = cos(pkin(9));
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t150 = -t77 * t71 + t80 * t73;
t26 = t150 * qJD(1);
t79 = cos(qJ(4));
t58 = -t79 * pkin(4) - pkin(3);
t19 = t58 * qJD(3) - t26;
t152 = t19 + t26;
t119 = qJDD(3) * pkin(3);
t34 = t80 * t71 + t77 * t73;
t29 = t34 * qJD(3);
t93 = -qJD(1) * t29 + qJDD(1) * t150;
t13 = -t93 - t119;
t81 = qJD(4) ^ 2;
t27 = t34 * qJD(1);
t72 = sin(pkin(8));
t74 = cos(pkin(8));
t105 = g(1) * t74 + g(2) * t72;
t66 = pkin(9) + qJ(3);
t59 = sin(t66);
t60 = cos(t66);
t90 = -g(3) * t60 + t105 * t59;
t86 = t27 * qJD(3) + t90;
t151 = -pkin(6) * t81 + t119 - t13 + t86;
t67 = qJD(4) + qJD(5);
t145 = pkin(6) + pkin(7);
t115 = t79 * qJDD(3);
t76 = sin(qJ(4));
t116 = t76 * qJDD(3);
t75 = sin(qJ(5));
t78 = cos(qJ(5));
t36 = t75 * t79 + t78 * t76;
t149 = t36 * t67;
t8 = qJD(3) * t149 - t78 * t115 + t75 * t116;
t35 = t75 * t76 - t78 * t79;
t94 = t35 * t67;
t108 = t145 * qJD(3) + t27;
t15 = t79 * qJD(2) - t108 * t76;
t148 = t67 * t79;
t147 = t34 * qJDD(1);
t16 = t76 * qJD(2) + t108 * t79;
t143 = g(3) * t59;
t125 = qJD(3) * t79;
t111 = t78 * t125;
t126 = qJD(3) * t76;
t112 = t75 * t126;
t30 = -t111 + t112;
t32 = -t75 * t125 - t78 * t126;
t141 = t32 * t30;
t65 = qJDD(4) + qJDD(5);
t140 = t35 * t65;
t139 = t36 * t65;
t70 = qJ(4) + qJ(5);
t63 = sin(t70);
t138 = t72 * t63;
t64 = cos(t70);
t137 = t72 * t64;
t136 = t74 * t63;
t135 = t74 * t64;
t132 = t78 * t16;
t68 = t76 ^ 2;
t129 = -t79 ^ 2 + t68;
t128 = qJD(3) * pkin(3);
t127 = qJD(3) * t150;
t124 = qJD(4) * t76;
t122 = qJD(5) * t75;
t121 = qJD(5) * t78;
t113 = qJD(3) * qJD(4);
t110 = qJD(4) * t145;
t109 = t79 * t113;
t12 = qJDD(3) * pkin(6) + qJD(1) * t127 + t147;
t107 = pkin(7) * qJDD(3) + t12;
t106 = pkin(4) * t124 - t27;
t104 = -g(1) * t72 + g(2) * t74;
t14 = qJD(4) * pkin(4) + t15;
t102 = -t75 * t14 - t132;
t101 = -t67 * t94 + t139;
t43 = t145 * t76;
t44 = t145 * t79;
t100 = -t78 * t43 - t75 * t44;
t99 = -t75 * t43 + t78 * t44;
t97 = -t29 * qJD(3) + qJDD(3) * t150;
t96 = -qJDD(2) - t104;
t92 = t34 * t81 - t97;
t91 = t105 * t60 + t143;
t88 = -0.2e1 * t127 * qJD(4) - qJDD(4) * t34;
t24 = -t26 - t128;
t87 = -pkin(6) * qJDD(4) + (t24 + t26 - t128) * qJD(4);
t7 = qJD(5) * t111 - t67 * t112 + t75 * t115 + (t109 + t116) * t78;
t85 = -t24 * qJD(3) - t12 + t91;
t61 = t79 * qJDD(2);
t2 = qJDD(4) * pkin(4) - t16 * qJD(4) - t107 * t76 + t61;
t84 = -g(1) * (-t60 * t135 - t138) - g(2) * (-t60 * t137 + t136) + t19 * t30 + t16 * t122 + t64 * t143 + (-t16 * t67 - t2) * t75;
t3 = t15 * qJD(4) + t76 * qJDD(2) + t107 * t79;
t83 = -g(1) * (-t60 * t136 + t137) - g(2) * (-t60 * t138 - t135) + t102 * qJD(5) + t19 * t32 + t78 * t2 - t75 * t3 + t63 * t143;
t82 = qJD(3) ^ 2;
t42 = qJDD(4) * t79 - t81 * t76;
t41 = qJDD(4) * t76 + t81 * t79;
t38 = t79 * t110;
t37 = t76 * t110;
t10 = -t30 ^ 2 + t32 ^ 2;
t9 = (t76 * t113 - t115) * pkin(4) + t13;
t6 = -t149 * t67 - t140;
t5 = -t32 * t67 - t8;
t4 = t30 * t67 + t7;
t1 = [qJDD(1) - g(3), -g(3) + (t71 ^ 2 + t73 ^ 2) * qJDD(1), 0, t97, -qJD(3) * t127 - t34 * qJDD(3), 0, 0, 0, 0, 0, t88 * t76 - t92 * t79, t92 * t76 + t88 * t79, 0, 0, 0, 0, 0, t29 * t30 - t150 * t8 - t127 * t149 + ((t76 * t122 + t75 * t124 - t148 * t78) * t67 - t139) * t34, -t29 * t32 - t150 * t7 + t127 * t94 + (-(-t76 * t121 - t78 * t124 - t148 * t75) * t67 + t140) * t34; 0, -t96, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t41, 0, 0, 0, 0, 0, t6, -t101; 0, 0, qJDD(3), t86 + t93, -t147 + t91, t68 * qJDD(3) + 0.2e1 * t76 * t109, -0.2e1 * t129 * t113 + 0.2e1 * t76 * t115, t41, t42, 0, t151 * t79 + t87 * t76, -t151 * t76 + t87 * t79, t32 * t94 + t7 * t36, t149 * t32 + t30 * t94 - t7 * t35 - t36 * t8, t101, t6, 0, (-t99 * qJD(5) + t75 * t37 - t78 * t38) * t67 + t100 * t65 + t58 * t8 + t9 * t35 + t106 * t30 + t90 * t64 + t152 * t149, -(t100 * qJD(5) - t78 * t37 - t75 * t38) * t67 - t99 * t65 + t58 * t7 + t9 * t36 - t106 * t32 - t90 * t63 - t152 * t94; 0, 0, 0, 0, 0, -t76 * t82 * t79, t129 * t82, t116, t115, qJDD(4), t104 * t79 + t85 * t76 + t61, t96 * t76 + t85 * t79, -t141, t10, t4, t5, t65, -(-t75 * t15 - t132) * t67 + (-t67 * t122 - t30 * t126 + t78 * t65) * pkin(4) + t83, (-qJD(5) * t14 + t15 * t67 - t3) * t78 + (-t67 * t121 + t32 * t126 - t75 * t65) * pkin(4) + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t10, t4, t5, t65, -t102 * t67 + t83, (-t3 + (-qJD(5) + t67) * t14) * t78 + t84;];
tau_reg = t1;
