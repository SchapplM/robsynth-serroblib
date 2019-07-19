% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPR2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:37
% EndTime: 2019-07-18 18:16:38
% DurationCPUTime: 0.34s
% Computational Cost: add. (765->123), mult. (884->148), div. (0->0), fcn. (430->8), ass. (0->77)
t59 = qJ(1) + qJ(2);
t54 = sin(t59);
t55 = cos(t59);
t102 = -g(1) * t55 - g(2) * t54;
t90 = pkin(1) * qJDD(1);
t66 = -pkin(2) - pkin(3);
t58 = qJD(1) + qJD(2);
t64 = cos(qJ(2));
t91 = pkin(1) * qJD(1);
t82 = t64 * t91;
t72 = qJD(3) - t82;
t20 = t66 * t58 + t72;
t61 = sin(qJ(2));
t83 = t61 * t91;
t27 = t58 * qJ(3) + t83;
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t5 = t63 * t20 - t60 * t27;
t53 = qJD(4) - t58;
t101 = t5 * t53;
t57 = qJDD(1) + qJDD(2);
t100 = t57 * pkin(2);
t6 = t60 * t20 + t63 * t27;
t99 = t6 * t53;
t62 = sin(qJ(1));
t98 = t62 * pkin(1);
t28 = -t60 * qJ(3) + t63 * t66;
t97 = t63 * qJD(3) + t28 * qJD(4) - (t60 * t61 + t63 * t64) * t91;
t29 = t63 * qJ(3) + t60 * t66;
t96 = -t60 * qJD(3) - t29 * qJD(4) - (-t60 * t64 + t61 * t63) * t91;
t80 = qJD(2) * t91;
t95 = t61 * t90 + t64 * t80;
t94 = t55 * pkin(2) + t54 * qJ(3);
t93 = -g(1) * t54 + g(2) * t55;
t92 = -t61 * t80 + t64 * t90;
t89 = qJD(2) * t61;
t88 = qJD(2) * t64;
t87 = qJD(4) * t60;
t86 = qJD(4) * t63;
t65 = cos(qJ(1));
t85 = t65 * pkin(1) + t94;
t84 = pkin(1) * t89;
t81 = t58 * t89;
t47 = -t64 * pkin(1) - pkin(2);
t79 = -qJDD(3) + t92;
t39 = t55 * qJ(3);
t78 = -t54 * pkin(2) + t39;
t11 = t66 * t57 - t79;
t50 = t57 * qJ(3);
t51 = t58 * qJD(3);
t12 = t50 + t51 + t95;
t2 = t63 * t11 - t60 * t12 - t20 * t87 - t27 * t86;
t77 = t53 ^ 2;
t76 = t92 - t93;
t75 = t95 + t102;
t74 = t66 * t54 + t39;
t73 = g(1) * t62 - g(2) * t65;
t37 = -pkin(3) + t47;
t40 = t61 * pkin(1) + qJ(3);
t16 = t63 * t37 - t60 * t40;
t17 = t60 * t37 + t63 * t40;
t71 = -qJDD(3) + t76;
t15 = -t79 - t100;
t1 = t60 * t11 + t63 * t12 + t20 * t86 - t27 * t87;
t23 = -t54 * t60 - t55 * t63;
t24 = -t54 * t63 + t55 * t60;
t70 = g(1) * t24 - g(2) * t23 + t2;
t69 = t58 * t82 - t75;
t68 = -g(1) * t23 - g(2) * t24 - t1;
t52 = -qJDD(4) + t57;
t41 = t55 * pkin(3);
t34 = pkin(1) * t88 + qJD(3);
t30 = t58 * t83;
t26 = -t58 * pkin(2) + t72;
t4 = -t17 * qJD(4) - t60 * t34 + t63 * t84;
t3 = t16 * qJD(4) + t63 * t34 + t60 * t84;
t7 = [0, 0, 0, 0, 0, qJDD(1), t73, g(1) * t65 + g(2) * t62, 0, 0, 0, 0, 0, 0, 0, t57, (t57 * t64 - t81) * pkin(1) + t76, (-t57 * t61 - t58 * t88) * pkin(1) - t75, 0, (t73 + (t61 ^ 2 + t64 ^ 2) * t90) * pkin(1), 0, 0, 0, t57, 0, 0, -pkin(1) * t81 + (pkin(2) - t47) * t57 + t71, 0, t34 * t58 + t40 * t57 + t102 + t12, t12 * t40 + t27 * t34 + t15 * t47 + t26 * t84 - g(1) * (t78 - t98) - g(2) * t85, 0, 0, 0, 0, 0, t52, -t16 * t52 + t4 * t53 - t70, t17 * t52 - t3 * t53 - t68, 0, t1 * t17 + t6 * t3 + t2 * t16 + t5 * t4 - g(1) * (t74 - t98) - g(2) * (t41 + t85); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t30 + t76, t69, 0, 0, 0, 0, 0, t57, 0, 0, t30 + t71 + 0.2e1 * t100, 0, 0.2e1 * t50 + 0.2e1 * t51 - t69, t12 * qJ(3) + t27 * qJD(3) - t15 * pkin(2) - g(1) * t78 - g(2) * t94 + (-t26 * t61 - t27 * t64) * t91, 0, 0, 0, 0, 0, t52, -t28 * t52 + t96 * t53 - t70, t29 * t52 - t97 * t53 - t68, 0, t1 * t29 + t2 * t28 - g(1) * t74 - g(2) * (t41 + t94) + t97 * t6 + t96 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, -t58 ^ 2, -t27 * t58 + t15 + t93, 0, 0, 0, 0, 0, 0, -t63 * t52 - t60 * t77, t60 * t52 - t63 * t77, 0, (t2 + t99) * t63 + (t1 - t101) * t60 + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t70 + t99, t68 + t101, 0, 0;];
tau_reg  = t7;
