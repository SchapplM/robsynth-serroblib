% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tau_reg [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:03
% EndTime: 2019-12-31 16:45:04
% DurationCPUTime: 0.52s
% Computational Cost: add. (697->147), mult. (1622->179), div. (0->0), fcn. (1135->8), ass. (0->88)
t107 = qJDD(1) * pkin(1);
t66 = sin(qJ(1));
t67 = cos(qJ(1));
t123 = g(1) * t66 - g(2) * t67;
t80 = -qJDD(2) + t107 + t123;
t87 = g(1) * t67 + g(2) * t66;
t114 = cos(qJ(3));
t62 = sin(pkin(6));
t63 = cos(pkin(6));
t65 = sin(qJ(3));
t37 = t114 * t62 + t65 * t63;
t32 = t37 * qJD(1);
t122 = qJ(2) * qJDD(1);
t110 = pkin(5) + qJ(2);
t42 = t110 * t62;
t43 = t110 * t63;
t19 = t114 * t43 - t65 * t42;
t61 = pkin(6) + qJ(3);
t56 = sin(t61);
t111 = t65 * t62;
t96 = t114 * t63;
t77 = t96 - t111;
t78 = -t114 * t42 - t65 * t43;
t8 = t77 * qJD(2) + t78 * qJD(3);
t121 = -t8 * qJD(3) - t19 * qJDD(3) - t123 * t56;
t57 = cos(t61);
t101 = qJD(1) * qJD(2);
t119 = t110 * qJDD(1) + t101;
t24 = t119 * t62;
t25 = t119 * t63;
t38 = qJD(1) * t42;
t92 = qJD(3) * t114;
t99 = t114 * t25 - t65 * t24 - t38 * t92;
t120 = -g(3) * t56 - t87 * t57 + t99;
t118 = t32 ^ 2;
t88 = qJD(1) * t96;
t97 = qJD(1) * t111;
t30 = -t88 + t97;
t113 = t32 * t30;
t39 = qJD(1) * t43;
t112 = t65 * t39;
t109 = t62 ^ 2 + t63 ^ 2;
t108 = qJD(3) * t65;
t106 = qJDD(3) * pkin(3);
t17 = t114 * t39 - t65 * t38;
t105 = t17 * qJD(3);
t16 = -t114 * t38 - t112;
t104 = qJD(4) - t16;
t103 = t62 * qJDD(1);
t102 = t63 * qJDD(1);
t100 = qJDD(3) * qJ(4);
t91 = qJDD(1) * t114;
t98 = qJD(3) * t88 + t65 * t102 + t62 * t91;
t54 = t63 * pkin(2) + pkin(1);
t94 = t109 * qJD(1) ^ 2;
t93 = t16 + t112;
t90 = -t38 * t108 + t114 * t24 + t65 * t25 + t39 * t92;
t89 = 0.2e1 * t109;
t85 = t65 * t103 - t63 * t91;
t84 = t57 * pkin(3) + t56 * qJ(4);
t79 = t54 + t84;
t41 = -t54 * qJD(1) + qJD(2);
t40 = -t54 * qJDD(1) + qJDD(2);
t76 = -g(3) * t57 + t87 * t56 - t90;
t35 = t37 * qJD(3);
t75 = t80 + t107;
t9 = t37 * qJD(2) + t19 * qJD(3);
t74 = -t9 * qJD(3) + t78 * qJDD(3) + t123 * t57;
t7 = t30 * pkin(3) - t32 * qJ(4) + t41;
t73 = t7 * t32 + qJDD(4) - t76;
t72 = t89 * t101 - t87;
t14 = qJD(3) * t97 - t98;
t15 = qJD(1) * t35 + t85;
t71 = t15 * pkin(3) + t14 * qJ(4) + t40;
t70 = 0.2e1 * t32 * qJD(3) + t85;
t34 = t62 * t108 - t63 * t92;
t29 = t30 ^ 2;
t13 = -pkin(3) * t77 - t37 * qJ(4) - t54;
t12 = t32 * pkin(3) + t30 * qJ(4);
t11 = qJD(3) * qJ(4) + t17;
t10 = -qJD(3) * pkin(3) + t104;
t6 = t35 * pkin(3) + t34 * qJ(4) - t37 * qJD(4);
t5 = (t30 - t97) * qJD(3) + t98;
t4 = (t30 + t97) * qJD(3) - t98;
t3 = qJDD(4) + t90 - t106;
t2 = t100 + (qJD(4) - t112) * qJD(3) + t99;
t1 = -t32 * qJD(4) + t71;
t18 = [qJDD(1), t123, t87, t75 * t63, -t75 * t62, t89 * t122 + t72, t80 * pkin(1) + (t109 * t122 + t72) * qJ(2), -t14 * t37 - t32 * t34, -t14 * t77 - t37 * t15 + t34 * t30 - t32 * t35, -t34 * qJD(3) + t37 * qJDD(3), -t35 * qJD(3) + qJDD(3) * t77, 0, -t54 * t15 + t41 * t35 - t40 * t77 + t74, t54 * t14 - t41 * t34 + t40 * t37 + t121, -t1 * t77 + t13 * t15 + t6 * t30 + t7 * t35 + t74, -t10 * t34 - t11 * t35 + t14 * t78 - t19 * t15 + t2 * t77 + t3 * t37 - t8 * t30 + t9 * t32 - t87, -t1 * t37 + t13 * t14 - t6 * t32 + t7 * t34 - t121, t1 * t13 + t10 * t9 + t11 * t8 - t3 * t78 + t2 * t19 + t7 * t6 + (-g(1) * t110 - g(2) * t79) * t67 + (g(1) * t79 - g(2) * t110) * t66; 0, 0, 0, -t102, t103, -t94, -qJ(2) * t94 - t80, 0, 0, 0, 0, 0, t70, -t4, t70, -t29 - t118, t4, t11 * t30 + (-qJD(4) - t10) * t32 + t71 - t123; 0, 0, 0, 0, 0, 0, 0, t113, -t29 + t118, t5, -t85, qJDD(3), -t41 * t32 + t105 + t76, t93 * qJD(3) + t41 * t30 - t120, -t12 * t30 + t105 + 0.2e1 * t106 - t73, pkin(3) * t14 - t15 * qJ(4) + (t11 - t17) * t32 + (t10 - t104) * t30, 0.2e1 * t100 + t12 * t32 - t7 * t30 + (0.2e1 * qJD(4) - t93) * qJD(3) + t120, -t3 * pkin(3) - g(3) * t84 + t2 * qJ(4) - t10 * t17 + t104 * t11 - t7 * t12 + t87 * (pkin(3) * t56 - qJ(4) * t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t113, t5, -qJD(3) ^ 2 - t118, -t11 * qJD(3) - t106 + t73;];
tau_reg = t18;
