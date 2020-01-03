% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPPR3
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
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:59
% EndTime: 2019-12-31 16:38:00
% DurationCPUTime: 0.50s
% Computational Cost: add. (668->129), mult. (1378->170), div. (0->0), fcn. (967->12), ass. (0->85)
t70 = cos(pkin(6));
t105 = t70 * pkin(1);
t53 = -pkin(2) - t105;
t95 = qJDD(1) * t53;
t40 = qJDD(3) + t95;
t66 = qJ(1) + pkin(6);
t59 = sin(t66);
t61 = cos(t66);
t88 = -g(1) * t59 + g(2) * t61;
t108 = -t40 - t88;
t67 = sin(pkin(7));
t72 = sin(qJ(4));
t100 = t72 * t67;
t69 = cos(pkin(7));
t74 = cos(qJ(4));
t38 = -t74 * t69 + t100;
t97 = pkin(1) * qJDD(1);
t39 = t74 * t67 + t72 * t69;
t30 = t39 * qJD(1);
t68 = sin(pkin(6));
t49 = t68 * pkin(1) + qJ(3);
t36 = qJD(1) * qJD(3) + qJDD(1) * t49;
t107 = t30 ^ 2;
t73 = sin(qJ(1));
t104 = t73 * pkin(1);
t103 = pkin(5) + t49;
t33 = t39 * qJD(4);
t93 = t69 * qJDD(1);
t94 = t67 * qJDD(1);
t83 = t72 * t94 - t74 * t93;
t11 = qJD(1) * t33 + t83;
t96 = qJD(1) * t69;
t89 = t74 * t96;
t90 = qJD(1) * t100;
t28 = -t89 + t90;
t32 = t38 * qJD(4);
t102 = -t39 * t11 + t32 * t28;
t101 = t30 * t28;
t21 = t67 * qJDD(2) + t69 * t36;
t42 = t49 * qJD(1);
t23 = t67 * qJD(2) + t69 * t42;
t63 = t67 ^ 2;
t64 = t69 ^ 2;
t98 = t63 + t64;
t91 = qJD(4) * t89 + t72 * t93 + t74 * t94;
t52 = t69 * pkin(3) + pkin(2);
t55 = t69 * qJDD(2);
t16 = t55 + (-pkin(5) * qJDD(1) - t36) * t67;
t17 = pkin(5) * t93 + t21;
t87 = t74 * t16 - t72 * t17;
t86 = g(1) * t61 + g(2) * t59;
t75 = cos(qJ(1));
t84 = g(1) * t73 - g(2) * t75;
t10 = qJD(4) * t90 - t91;
t82 = -t38 * t10 + t30 * t33;
t81 = t72 * t16 + t74 * t17;
t57 = t69 * qJD(2);
t18 = t57 + (-pkin(5) * qJD(1) - t42) * t67;
t19 = pkin(5) * t96 + t23;
t3 = t74 * t18 - t72 * t19;
t4 = t72 * t18 + t74 * t19;
t20 = -t67 * t36 + t55;
t80 = -t20 * t67 + t21 * t69;
t79 = (-t67 * t42 + t57) * t67 - t23 * t69;
t34 = t103 * t67;
t35 = t103 * t69;
t8 = -t74 * t34 - t72 * t35;
t9 = -t72 * t34 + t74 * t35;
t41 = -t52 - t105;
t78 = -t95 + t108;
t24 = t41 * qJDD(1) + qJDD(3);
t71 = -pkin(5) - qJ(3);
t65 = pkin(7) + qJ(4);
t62 = t75 * pkin(1);
t60 = cos(t65);
t58 = sin(t65);
t27 = t28 ^ 2;
t26 = t41 * qJD(1) + qJD(3);
t13 = -t33 * qJD(4) - t38 * qJDD(4);
t12 = -t32 * qJD(4) + t39 * qJDD(4);
t6 = -t39 * qJD(3) - t9 * qJD(4);
t5 = -t38 * qJD(3) + t8 * qJD(4);
t2 = -t4 * qJD(4) + t87;
t1 = t3 * qJD(4) + t81;
t7 = [0, 0, 0, 0, 0, qJDD(1), t84, g(1) * t75 + g(2) * t73, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t70 * t97 - t88, -0.2e1 * t68 * t97 + t86, 0, (t84 + (t68 ^ 2 + t70 ^ 2) * t97) * pkin(1), t63 * qJDD(1), 0.2e1 * t67 * t93, 0, t64 * qJDD(1), 0, 0, t78 * t69, -t78 * t67, t36 * t98 + t80 - t86, t40 * t53 - g(1) * (-t59 * pkin(2) + t61 * qJ(3) - t104) - g(2) * (t61 * pkin(2) + t59 * qJ(3) + t62) + t80 * t49 - t79 * qJD(3), -t10 * t39 - t30 * t32, -t82 + t102, t12, t11 * t38 + t28 * t33, t13, 0, t6 * qJD(4) + t8 * qJDD(4) + t41 * t11 + t24 * t38 + t26 * t33 - t60 * t88, -t5 * qJD(4) - t9 * qJDD(4) - t41 * t10 + t24 * t39 - t26 * t32 + t58 * t88, -t1 * t38 + t8 * t10 - t9 * t11 - t2 * t39 - t5 * t28 + t3 * t32 - t6 * t30 - t4 * t33 - t86, t1 * t9 + t4 * t5 + t2 * t8 + t3 * t6 + t24 * t41 - g(1) * (-t59 * t52 - t61 * t71 - t104) - g(2) * (t61 * t52 - t59 * t71 + t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t69 + t21 * t67 - g(3), 0, 0, 0, 0, 0, 0, t13, -t12, t82 + t102, t1 * t39 - t2 * t38 - t3 * t33 - t4 * t32 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t94, -t98 * qJD(1) ^ 2, t79 * qJD(1) - t108, 0, 0, 0, 0, 0, 0, 0.2e1 * t30 * qJD(4) + t83, (-t28 - t90) * qJD(4) + t91, -t27 - t107, t4 * t28 + t3 * t30 + t24 + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t27 + t107, (t28 - t90) * qJD(4) + t91, -t101, -t83, qJDD(4), -g(3) * t60 - t26 * t30 + t86 * t58 + t87, g(3) * t58 + t26 * t28 + t86 * t60 - t81, 0, 0;];
tau_reg = t7;
