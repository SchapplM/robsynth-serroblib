% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:17
% EndTime: 2019-12-31 17:45:18
% DurationCPUTime: 0.42s
% Computational Cost: add. (447->113), mult. (742->143), div. (0->0), fcn. (498->12), ass. (0->74)
t63 = sin(pkin(7));
t41 = pkin(1) * t63 + qJ(3);
t99 = qJDD(1) * t41;
t62 = sin(pkin(8));
t64 = cos(pkin(8));
t98 = t62 ^ 2 + t64 ^ 2;
t57 = qJ(1) + pkin(7);
t50 = sin(t57);
t52 = cos(t57);
t87 = -g(1) * t50 + g(2) * t52;
t65 = cos(pkin(7));
t44 = -pkin(1) * t65 - pkin(2);
t105 = qJDD(1) * t44;
t66 = sin(qJ(5));
t68 = cos(qJ(5));
t29 = t62 * t68 + t64 * t66;
t97 = qJD(1) * t62;
t89 = t66 * t97;
t73 = -qJD(5) * t89 + t29 * qJDD(1);
t38 = -qJ(4) + t44;
t102 = -pkin(6) + t38;
t101 = t62 * t66;
t100 = t64 * t68;
t33 = t41 * qJD(1);
t94 = qJD(1) * qJD(4);
t16 = qJDD(1) * t38 + qJDD(3) - t94;
t10 = t64 * qJDD(2) + t62 * t16;
t96 = t62 * qJDD(1);
t95 = t64 * qJDD(1);
t59 = qJD(3) * qJD(1);
t69 = cos(qJ(1));
t93 = t69 * pkin(1) + t52 * pkin(2) + t50 * qJ(3);
t92 = -t59 - t99;
t91 = qJD(1) * t100;
t90 = qJD(5) * t100;
t31 = qJD(4) + t33;
t67 = sin(qJ(1));
t88 = -pkin(1) * t67 + t52 * qJ(3);
t86 = t98 * qJDD(1);
t9 = -qJDD(2) * t62 + t64 * t16;
t27 = qJDD(4) - t92;
t85 = -g(1) * t52 - g(2) * t50;
t83 = g(1) * t67 - g(2) * t69;
t82 = -t66 * t96 + t68 * t95;
t81 = t10 * t62 + t9 * t64;
t80 = t98 * (t38 * qJD(1) + qJD(3));
t20 = t102 * t62;
t21 = t102 * t64;
t79 = t20 * t68 + t21 * t66;
t78 = t20 * t66 - t21 * t68;
t77 = -t100 + t101;
t25 = t29 * qJD(5);
t3 = -qJD(5) * t25 - qJDD(5) * t77;
t26 = -qJD(5) * t101 + t90;
t4 = -qJD(5) * t26 - qJDD(5) * t29;
t22 = t29 * qJD(1);
t76 = qJDD(3) + t105;
t75 = t85 + t99;
t74 = -t81 - t87;
t72 = t27 + t75 + t59;
t70 = qJD(1) ^ 2;
t61 = qJDD(2) - g(3);
t56 = pkin(8) + qJ(5);
t51 = cos(t56);
t49 = sin(t56);
t32 = pkin(4) * t62 + t41;
t24 = -t89 + t91;
t19 = pkin(4) * t97 + t31;
t15 = pkin(4) * t96 + t27;
t6 = -pkin(6) * t96 + t10;
t5 = -pkin(6) * t95 + t9;
t2 = qJD(1) * t90 + t73;
t1 = -qJD(1) * t25 + t82;
t7 = [qJDD(1), t83, g(1) * t69 + g(2) * t67, (t83 + (t63 ^ 2 + t65 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(3) + t87 + 0.2e1 * t105, 0.2e1 * t59 + t75 + t99, -t92 * t41 + t33 * qJD(3) + t76 * t44 - g(1) * (-pkin(2) * t50 + t88) - g(2) * t93, t72 * t62, t72 * t64, -t38 * t86 + t98 * t94 + t74, t27 * t41 + t31 * qJD(3) - g(1) * ((-pkin(2) - qJ(4)) * t50 + t88) - g(2) * (qJ(4) * t52 + t93) + t81 * t38 - t80 * qJD(4), -t1 * t77 - t24 * t25, -t1 * t29 + t2 * t77 + t22 * t25 - t24 * t26, t3, t4, 0, qJD(3) * t22 + t32 * t2 + t15 * t29 + t19 * t26 - t78 * qJDD(5) + t85 * t49 + (t77 * qJD(4) - t79 * qJD(5)) * qJD(5), qJD(3) * t24 + t32 * t1 - t15 * t77 - t19 * t25 - t79 * qJDD(5) + t85 * t51 + (t29 * qJD(4) + t78 * qJD(5)) * qJD(5); 0, 0, 0, t61, 0, 0, t61, 0, 0, 0, t10 * t64 - t62 * t9 - g(3), 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, qJDD(1), -t70, -qJD(1) * t33 + t76 + t87, -t70 * t62, -t70 * t64, -t86, -qJD(1) * t31 - t74, 0, 0, 0, 0, 0, -qJD(1) * t22 + t3, -qJD(1) * t24 + t4; 0, 0, 0, 0, 0, 0, 0, t96, t95, -t98 * t70, t80 * qJD(1) + t27 + t85, 0, 0, 0, 0, 0, (t24 + t91) * qJD(5) + t73, -0.2e1 * qJD(5) * t22 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t22, -t22 ^ 2 + t24 ^ 2, t82, (t24 - t91) * qJD(5) - t73, qJDD(5), g(3) * t49 - t19 * t24 + t68 * t5 + t51 * t87 - t66 * t6, g(3) * t51 + t19 * t22 - t49 * t87 - t66 * t5 - t68 * t6;];
tau_reg = t7;
