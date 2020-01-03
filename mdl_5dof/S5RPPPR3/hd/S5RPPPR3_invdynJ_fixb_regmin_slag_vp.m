% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:03
% EndTime: 2019-12-31 17:44:05
% DurationCPUTime: 0.52s
% Computational Cost: add. (478->136), mult. (897->176), div. (0->0), fcn. (633->10), ass. (0->88)
t71 = sin(pkin(7));
t55 = t71 * pkin(1) + qJ(3);
t38 = qJD(1) * qJD(3) + qJDD(1) * t55;
t70 = sin(pkin(8));
t67 = t70 ^ 2;
t72 = cos(pkin(8));
t116 = t72 ^ 2 + t67;
t69 = qJ(1) + pkin(7);
t64 = sin(t69);
t114 = g(1) * t64;
t65 = cos(t69);
t59 = g(2) * t65;
t115 = t59 - t114;
t74 = sin(qJ(5));
t76 = cos(qJ(5));
t41 = t70 * t74 + t72 * t76;
t32 = t41 * qJD(1);
t106 = t72 * t74;
t42 = t70 * t76 - t106;
t95 = qJD(1) * t106;
t83 = qJD(5) * t95 - t41 * qJDD(1);
t113 = t72 * pkin(3);
t112 = -pkin(6) + t55;
t18 = t70 * qJDD(2) + t72 * t38;
t111 = t18 * t70 - g(3);
t47 = t55 * qJD(1);
t109 = (t70 * qJD(2) + t72 * t47) * t72;
t12 = t18 * t72;
t110 = qJD(3) * t109 + t55 * t12;
t108 = t64 * t72;
t105 = qJD(1) * t70;
t104 = qJD(4) * t70;
t103 = qJDD(2) - g(3);
t73 = cos(pkin(7));
t60 = -t73 * pkin(1) - pkin(2);
t86 = t70 * qJ(4) - t60;
t35 = -t86 - t113;
t102 = qJDD(1) * t35;
t100 = qJDD(1) * t60;
t62 = t70 * qJDD(1);
t99 = t72 * qJDD(1);
t97 = qJD(1) * qJD(4);
t77 = cos(qJ(1));
t96 = t77 * pkin(1) + t65 * pkin(2) + t64 * qJ(3);
t94 = t76 * t105;
t93 = t70 * t97;
t26 = t72 * qJD(2) - t70 * t47;
t28 = t70 * t38;
t17 = t72 * qJDD(2) - t28;
t92 = -g(1) * t65 - g(2) * t64;
t75 = sin(qJ(1));
t91 = g(1) * t75 - g(2) * t77;
t90 = t76 * t62 - t74 * t99;
t36 = t112 * t70;
t37 = t112 * t72;
t89 = t76 * t36 - t74 * t37;
t88 = t74 * t36 + t76 * t37;
t15 = qJDD(4) - t17;
t87 = -t75 * pkin(1) - t64 * pkin(2) + t65 * qJ(3);
t81 = qJDD(3) + t102;
t8 = t81 - t93;
t85 = -t102 - t8 - t59;
t45 = qJDD(3) + t100;
t84 = -t100 - t45 - t59;
t30 = t41 * qJD(5);
t82 = t38 * t116 + t12 + t92;
t25 = (pkin(3) + pkin(4)) * t72 + t86;
t79 = qJD(1) ^ 2;
t78 = qJD(5) ^ 2;
t48 = g(1) * t108;
t46 = t116 * t79;
t34 = t94 - t95;
t31 = t42 * qJD(5);
t24 = qJD(4) - t26;
t23 = t41 * t65;
t22 = t42 * t65;
t21 = t41 * t64;
t20 = t42 * t64;
t19 = t35 * qJD(1) + qJD(3);
t10 = t25 * qJD(1) - qJD(3);
t9 = -pkin(6) * t99 + t18;
t7 = -pkin(6) * t62 + t15;
t5 = -t31 * qJD(5) - t41 * qJDD(5);
t4 = -t30 * qJD(5) + t42 * qJDD(5);
t3 = qJD(5) * t94 - t83;
t2 = -qJD(1) * t30 + t90;
t1 = t25 * qJDD(1) - qJDD(3) + t93;
t6 = [qJDD(1), t91, g(1) * t77 + g(2) * t75, (t91 + (t71 ^ 2 + t73 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t84 * t72 + t48, (-t84 - t114) * t70, -t17 * t70 + t82, t45 * t60 - g(1) * t87 - g(2) * t96 + (-t26 * qJD(3) - t17 * t55) * t70 + t110, t48 + (t85 + t93) * t72, t15 * t70 + t82, t67 * t97 + (t85 + t114) * t70, t8 * t35 - g(1) * (-pkin(3) * t108 + t87) - g(2) * (t65 * t113 + t96) + (-t115 * qJ(4) + t24 * qJD(3) - t19 * qJD(4) + t15 * t55) * t70 + t110, t2 * t42 - t34 * t30, -t2 * t41 - t42 * t3 + t30 * t32 - t34 * t31, t4, t5, 0, t32 * t104 + t25 * t3 + t1 * t41 + t10 * t31 + t89 * qJDD(5) + g(1) * t21 - g(2) * t23 + (t42 * qJD(3) - t88 * qJD(5)) * qJD(5), t34 * t104 + t25 * t2 + t1 * t42 - t10 * t30 - t88 * qJDD(5) + g(1) * t20 - g(2) * t22 + (-t41 * qJD(3) - t89 * qJD(5)) * qJD(5); 0, 0, 0, t103, 0, 0, 0, t17 * t72 + t111, 0, 0, 0, -t15 * t72 + t111, 0, 0, 0, 0, 0, t5, -t4; 0, 0, 0, 0, -t99, t62, -t46, (t26 * t70 - t109) * qJD(1) + t45 + t115, -t99, -t46, -t62, (-t109 + (-qJD(4) - t24) * t70) * qJD(1) + t81 + t115, 0, 0, 0, 0, 0, (-t34 - t94) * qJD(5) + t83, 0.2e1 * t32 * qJD(5) - t90; 0, 0, 0, 0, 0, 0, 0, 0, -t70 * t79 * t72, t62, -t67 * t79, qJDD(4) + t28 - t103 * t72 + (qJD(1) * t19 + t92) * t70, 0, 0, 0, 0, 0, qJDD(5) * t76 - t32 * t105 - t78 * t74, -qJDD(5) * t74 - t34 * t105 - t78 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t32, -t32 ^ 2 + t34 ^ 2, t90, (t34 - t94) * qJD(5) + t83, qJDD(5), -g(1) * t22 - g(2) * t20 + g(3) * t41 - t10 * t34 + t76 * t7 - t74 * t9, g(1) * t23 + g(2) * t21 + g(3) * t42 + t10 * t32 - t74 * t7 - t76 * t9;];
tau_reg = t6;
