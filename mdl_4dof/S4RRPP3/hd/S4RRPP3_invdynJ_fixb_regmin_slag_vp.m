% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPP3
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% tau_reg [4x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:52
% EndTime: 2019-12-31 16:57:54
% DurationCPUTime: 0.53s
% Computational Cost: add. (802->170), mult. (1872->223), div. (0->0), fcn. (1185->8), ass. (0->89)
t68 = sin(pkin(6));
t69 = cos(pkin(6));
t71 = sin(qJ(2));
t73 = cos(qJ(2));
t45 = t68 * t73 + t69 * t71;
t38 = t45 * qJD(1);
t34 = t38 ^ 2;
t106 = qJD(1) * t71;
t109 = t69 * t73;
t97 = qJD(1) * t109;
t35 = t68 * t106 - t97;
t118 = -t35 ^ 2 - t34;
t72 = sin(qJ(1));
t74 = cos(qJ(1));
t117 = g(1) * t72 - g(2) * t74;
t90 = g(1) * t74 + g(2) * t72;
t64 = qJ(2) + pkin(6);
t61 = sin(t64);
t62 = cos(t64);
t111 = t73 * pkin(2);
t60 = pkin(1) + t111;
t49 = -t60 * qJD(1) + qJD(3);
t9 = t35 * pkin(3) - t38 * qJ(4) + t49;
t108 = qJ(3) + pkin(5);
t91 = qJD(2) * t108;
t83 = -t71 * qJD(3) - t73 * t91;
t94 = t108 * t71;
t19 = qJDD(2) * pkin(2) + t83 * qJD(1) - qJDD(1) * t94;
t32 = t73 * qJD(3) - t71 * t91;
t50 = t108 * t73;
t25 = t32 * qJD(1) + qJDD(1) * t50;
t4 = t69 * t19 - t68 * t25;
t96 = -qJDD(4) + t4;
t116 = -g(3) * t62 - t9 * t38 + t90 * t61 + t96;
t100 = qJD(1) * qJD(2);
t93 = t71 * t100;
t81 = t45 * qJDD(1) - t68 * t93;
t112 = g(3) * t73;
t48 = qJD(1) * t50;
t110 = t68 * t48;
t41 = t69 * t48;
t5 = t68 * t19 + t69 * t25;
t47 = qJD(1) * t94;
t43 = qJD(2) * pkin(2) - t47;
t22 = t68 * t43 + t41;
t66 = t71 ^ 2;
t107 = -t73 ^ 2 + t66;
t105 = qJD(2) * t71;
t104 = qJDD(2) * pkin(3);
t27 = -t69 * t47 - t110;
t103 = qJD(4) - t27;
t102 = t71 * qJDD(1);
t101 = t73 * qJDD(1);
t99 = qJDD(2) * qJ(4) + t5;
t98 = pkin(2) * t105;
t92 = t73 * t100;
t88 = -t69 * t101 + t68 * t102;
t87 = t62 * pkin(3) + t61 * qJ(4);
t21 = t69 * t43 - t110;
t86 = -0.2e1 * pkin(1) * t100 - pkin(5) * qJDD(2);
t37 = t45 * qJD(2);
t82 = pkin(2) * t93 - t60 * qJDD(1) + qJDD(3);
t75 = qJD(2) ^ 2;
t80 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t75 + t117;
t76 = qJD(1) ^ 2;
t79 = pkin(1) * t76 - pkin(5) * qJDD(1) + t90;
t11 = t68 * t32 - t69 * t83;
t12 = t69 * t32 + t68 * t83;
t23 = qJD(1) * t37 + t88;
t24 = t69 * t92 + t81;
t28 = t68 * t50 + t69 * t94;
t29 = t69 * t50 - t68 * t94;
t78 = t11 * t38 - t12 * t35 - t29 * t23 + t28 * t24 - t90;
t77 = t23 * pkin(3) - t24 * qJ(4) + t82;
t59 = -t69 * pkin(2) - pkin(3);
t57 = t68 * pkin(2) + qJ(4);
t52 = t74 * t60;
t44 = t68 * t71 - t109;
t40 = qJD(2) * t109 - t68 * t105;
t26 = -t68 * t47 + t41;
t20 = t44 * pkin(3) - t45 * qJ(4) - t60;
t16 = qJD(2) * qJ(4) + t22;
t13 = -qJD(2) * pkin(3) + qJD(4) - t21;
t10 = pkin(2) * t106 + t38 * pkin(3) + t35 * qJ(4);
t8 = t37 * pkin(3) - t40 * qJ(4) - t45 * qJD(4) + t98;
t3 = -t96 - t104;
t2 = qJD(2) * qJD(4) + t99;
t1 = -t38 * qJD(4) + t77;
t6 = [qJDD(1), t117, t90, t66 * qJDD(1) + 0.2e1 * t71 * t92, -0.2e1 * t107 * t100 + 0.2e1 * t71 * t101, qJDD(2) * t71 + t75 * t73, qJDD(2) * t73 - t75 * t71, 0, t86 * t71 + t80 * t73, -t80 * t71 + t86 * t73, -t21 * t40 - t22 * t37 - t4 * t45 - t5 * t44 + t78, t5 * t29 + t22 * t12 - t4 * t28 - t21 * t11 - t82 * t60 + t49 * t98 - g(1) * (t108 * t74 - t72 * t60) - g(2) * (t108 * t72 + t52), -t11 * qJD(2) - t28 * qJDD(2) + t1 * t44 + t117 * t62 + t20 * t23 + t8 * t35 + t9 * t37, t13 * t40 - t16 * t37 - t2 * t44 + t3 * t45 + t78, t12 * qJD(2) + t29 * qJDD(2) - t1 * t45 + t117 * t61 - t20 * t24 - t8 * t38 - t9 * t40, -g(2) * t52 + t1 * t20 + t13 * t11 + t16 * t12 + t2 * t29 + t3 * t28 + t9 * t8 + (-g(1) * t108 - g(2) * t87) * t74 + (-g(1) * (-t60 - t87) - g(2) * t108) * t72; 0, 0, 0, -t71 * t76 * t73, t107 * t76, t102, t101, qJDD(2), t79 * t71 - t112, g(3) * t71 + t79 * t73, (t22 - t26) * t38 + (-t21 + t27) * t35 + (-t23 * t68 - t24 * t69) * pkin(2), t21 * t26 - t22 * t27 + (-t112 + t4 * t69 + t5 * t68 + (-qJD(1) * t49 + t90) * t71) * pkin(2), t26 * qJD(2) - t10 * t35 + (pkin(3) - t59) * qJDD(2) + t116, -t57 * t23 + t59 * t24 + (t16 - t26) * t38 + (t13 - t103) * t35, -g(3) * t61 + t57 * qJDD(2) + t10 * t38 - t9 * t35 - t90 * t62 + (0.2e1 * qJD(4) - t27) * qJD(2) + t99, t2 * t57 + t3 * t59 - t9 * t10 - t13 * t26 - g(3) * (t87 + t111) + t103 * t16 + t90 * (pkin(2) * t71 + pkin(3) * t61 - qJ(4) * t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t21 * t38 + t22 * t35 - t117 + t82, 0.2e1 * t38 * qJD(2) + t88, t118, (t35 - t97) * qJD(2) - t81, t16 * t35 + (-qJD(4) - t13) * t38 + t77 - t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t35 - qJDD(2), (t35 + t97) * qJD(2) + t81, -t34 - t75, -t16 * qJD(2) - t104 - t116;];
tau_reg = t6;
