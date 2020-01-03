% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRPR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:58
% EndTime: 2019-12-31 16:20:58
% DurationCPUTime: 0.39s
% Computational Cost: add. (534->109), mult. (1160->146), div. (0->0), fcn. (851->8), ass. (0->73)
t88 = qJDD(2) * pkin(2);
t61 = pkin(6) + qJ(2);
t55 = sin(t61);
t57 = cos(t61);
t99 = g(1) * t55 - g(2) * t57;
t71 = -qJDD(3) + t88 + t99;
t83 = qJD(2) * qJD(3);
t98 = qJ(3) * qJDD(2) + t83;
t64 = cos(pkin(7));
t67 = cos(qJ(4));
t63 = sin(pkin(7));
t66 = sin(qJ(4));
t93 = t66 * t63;
t29 = -t67 * t64 + t93;
t30 = t67 * t63 + t66 * t64;
t25 = t30 * qJD(2);
t97 = t25 ^ 2;
t89 = qJD(2) * t64;
t80 = t67 * t89;
t81 = qJD(2) * t93;
t23 = -t80 + t81;
t27 = t29 * qJD(4);
t28 = t30 * qJD(4);
t85 = t64 * qJDD(2);
t86 = t63 * qJDD(2);
t75 = t66 * t86 - t67 * t85;
t9 = qJD(2) * t28 + t75;
t95 = t27 * t23 - t30 * t9;
t94 = t25 * t23;
t91 = pkin(5) + qJ(3);
t87 = qJ(3) * qJD(2);
t32 = t63 * qJD(1) + t64 * t87;
t58 = t63 ^ 2;
t59 = t64 ^ 2;
t90 = t58 + t59;
t82 = qJD(4) * t80 + t66 * t85 + t67 * t86;
t19 = t63 * qJDD(1) + t98 * t64;
t48 = t64 * pkin(3) + pkin(2);
t35 = t91 * t63;
t51 = t64 * qJDD(1);
t16 = t51 + (-t91 * qJDD(2) - t83) * t63;
t17 = pkin(5) * t85 + t19;
t78 = t67 * t16 - t66 * t17;
t77 = g(1) * t57 + g(2) * t55;
t8 = qJD(4) * t81 - t82;
t74 = t25 * t28 - t29 * t8;
t73 = t66 * t16 + t67 * t17;
t53 = t64 * qJD(1);
t20 = -qJD(2) * t35 + t53;
t21 = pkin(5) * t89 + t32;
t6 = t67 * t20 - t66 * t21;
t7 = t66 * t20 + t67 * t21;
t72 = (-t63 * t87 + t53) * t63 - t32 * t64;
t36 = t91 * t64;
t14 = -t67 * t35 - t66 * t36;
t15 = -t66 * t35 + t67 * t36;
t33 = -t48 * qJDD(2) + qJDD(3);
t70 = t71 + t88;
t18 = -t98 * t63 + t51;
t69 = -t18 * t63 + t19 * t64 - t77;
t62 = qJDD(1) - g(3);
t60 = pkin(7) + qJ(4);
t56 = cos(t60);
t54 = sin(t60);
t34 = -t48 * qJD(2) + qJD(3);
t22 = t23 ^ 2;
t11 = -t28 * qJD(4) - t29 * qJDD(4);
t10 = -t27 * qJD(4) + t30 * qJDD(4);
t5 = -t30 * qJD(3) - t15 * qJD(4);
t4 = -t29 * qJD(3) + t14 * qJD(4);
t2 = -t7 * qJD(4) + t78;
t1 = t6 * qJD(4) + t73;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t64 + t19 * t63 - g(3), 0, 0, 0, 0, 0, 0, t11, -t10, t74 + t95, t1 * t30 - t2 * t29 - t7 * t27 - t6 * t28 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t99, t77, 0, 0, t58 * qJDD(2), 0.2e1 * t63 * t85, 0, t59 * qJDD(2), 0, 0, t70 * t64, -t70 * t63, t98 * t90 + t69, t71 * pkin(2) + t69 * qJ(3) - t72 * qJD(3), -t25 * t27 - t8 * t30, -t74 + t95, t10, t23 * t28 + t9 * t29, t11, 0, t5 * qJD(4) + t14 * qJDD(4) + t34 * t28 + t33 * t29 - t48 * t9 + t56 * t99, -t4 * qJD(4) - t15 * qJDD(4) - t34 * t27 + t33 * t30 + t48 * t8 - t54 * t99, -t1 * t29 + t14 * t8 - t15 * t9 - t2 * t30 - t4 * t23 - t5 * t25 + t6 * t27 - t7 * t28 - t77, t1 * t15 + t7 * t4 + t2 * t14 + t6 * t5 - t33 * t48 - g(1) * (-t55 * t48 + t57 * t91) - g(2) * (t57 * t48 + t55 * t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t86, -t90 * qJD(2) ^ 2, t72 * qJD(2) - t71, 0, 0, 0, 0, 0, 0, 0.2e1 * t25 * qJD(4) + t75, (-t23 - t81) * qJD(4) + t82, -t22 - t97, t7 * t23 + t6 * t25 + t33 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t22 + t97, (t23 - t81) * qJD(4) + t82, -t94, -t75, qJDD(4), -g(3) * t56 - t34 * t25 + t77 * t54 + t78, g(3) * t54 + t34 * t23 + t77 * t56 - t73, 0, 0;];
tau_reg = t3;
