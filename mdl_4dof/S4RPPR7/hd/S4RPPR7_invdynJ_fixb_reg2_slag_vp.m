% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPPR7
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:44
% EndTime: 2019-12-31 16:41:44
% DurationCPUTime: 0.48s
% Computational Cost: add. (660->127), mult. (1229->160), div. (0->0), fcn. (785->8), ass. (0->83)
t57 = -pkin(1) - qJ(3);
t34 = t57 * qJD(1) + qJD(2);
t54 = sin(pkin(6));
t48 = t54 ^ 2;
t55 = cos(pkin(6));
t49 = t55 ^ 2;
t88 = t48 + t49;
t103 = t34 * t88;
t59 = sin(qJ(1));
t61 = cos(qJ(1));
t102 = g(1) * t59 - g(2) * t61;
t51 = (qJDD(1) * qJ(2));
t52 = (qJD(1) * qJD(2));
t101 = t51 + t52;
t32 = qJDD(3) + t101;
t75 = g(1) * t61 + g(2) * t59;
t66 = -t75 + t32;
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t92 = t60 * t55;
t100 = -t58 * t54 + t92;
t26 = t60 * t54 + t58 * t55;
t86 = qJD(1) * t54;
t82 = t58 * t86;
t65 = -qJD(4) * t82 + t26 * qJDD(1);
t99 = -qJD(1) * qJD(3) + qJDD(1) * t57;
t81 = qJD(1) * t92;
t20 = t81 - t82;
t98 = t20 ^ 2;
t97 = 2 * t52;
t96 = t54 * pkin(3);
t95 = -pkin(5) + t57;
t18 = t26 * qJD(1);
t94 = t20 * t18;
t22 = t26 * qJD(4);
t91 = -t22 * qJD(4) + qJDD(4) * t100;
t90 = t61 * pkin(1) + t59 * qJ(2);
t87 = pkin(1) * qJDD(1);
t40 = qJD(1) * qJ(2) + qJD(3);
t85 = t54 * qJDD(1);
t84 = t55 * qJDD(1);
t28 = qJDD(2) + t99;
t80 = t88 * t28;
t77 = -pkin(5) * qJDD(1) + t28;
t12 = t77 * t54;
t13 = t77 * t55;
t79 = -t58 * t12 + t60 * t13;
t78 = -pkin(5) * qJD(1) + t34;
t76 = qJDD(2) - t87;
t73 = -t58 * t85 + t60 * t84;
t23 = t100 * qJD(4);
t8 = qJD(4) * t81 + t65;
t72 = t23 * t18 + t26 * t8;
t7 = qJD(1) * t22 - t73;
t71 = -t100 * t7 - t20 * t22;
t70 = t60 * t12 + t58 * t13;
t14 = t78 * t54;
t15 = t78 * t55;
t6 = t60 * t14 + t58 * t15;
t5 = -t58 * t14 + t60 * t15;
t29 = t95 * t54;
t30 = t95 * t55;
t10 = t60 * t29 + t58 * t30;
t9 = -t58 * t29 + t60 * t30;
t68 = -t23 * qJD(4) - t26 * qJDD(4);
t1 = t5 * qJD(4) + t70;
t2 = -t6 * qJD(4) + t79;
t64 = t1 * t26 + t100 * t2 - t5 * t22 + t6 * t23 - t102;
t63 = t66 + t101;
t62 = qJD(1) ^ 2;
t56 = -pkin(5) - qJ(3);
t50 = pkin(6) + qJ(4);
t44 = t61 * qJ(2);
t42 = cos(t50);
t41 = sin(t50);
t39 = pkin(3) * t85;
t38 = qJ(2) + t96;
t31 = pkin(3) * t86 + t40;
t25 = t39 + t32;
t17 = t18 ^ 2;
t4 = -qJD(3) * t100 - t10 * qJD(4);
t3 = -t26 * qJD(3) + t9 * qJD(4);
t11 = [0, 0, 0, 0, 0, qJDD(1), t102, t75, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - 0.2e1 * t87 - t102, (2 * t51) + t97 - t75, -t76 * pkin(1) - g(1) * (-t59 * pkin(1) + t44) - g(2) * t90 + (t51 + t97) * qJ(2), t49 * qJDD(1), -0.2e1 * t54 * t84, 0, t48 * qJDD(1), 0, 0, t63 * t54, t63 * t55, t102 + t88 * (-t28 - t99), t32 * qJ(2) + t40 * qJD(2) - g(1) * (t57 * t59 + t44) - g(2) * (t61 * qJ(3) + t90) + t57 * t80 - qJD(3) * t103, t71, -t100 * t8 + t22 * t18 - t20 * t23 + t7 * t26, t91, t72, t68, 0, qJD(2) * t18 + t4 * qJD(4) + t9 * qJDD(4) + t31 * t23 + t25 * t26 + t38 * t8 - t75 * t41, qJD(2) * t20 - t3 * qJD(4) - t10 * qJDD(4) + t100 * t25 - t31 * t22 - t38 * t7 - t75 * t42, -t10 * t8 - t3 * t18 - t4 * t20 + t9 * t7 - t64, t1 * t10 + t6 * t3 + t2 * t9 + t5 * t4 + t25 * t38 + t31 * qJD(2) - g(1) * (t61 * t96 + t44 + (-pkin(1) + t56) * t59) - g(2) * (-t61 * t56 + t59 * t96 + t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t62, -t62 * qJ(2) - t102 + t76, 0, 0, 0, 0, 0, 0, -t62 * t54, -t62 * t55, -t88 * qJDD(1), -t40 * qJD(1) - t102 + t80, 0, 0, 0, 0, 0, 0, -qJD(1) * t18 + t91, -qJD(1) * t20 + t68, -t71 - t72, -t31 * qJD(1) + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t84, -t88 * t62, qJD(1) * t103 + t66, 0, 0, 0, 0, 0, 0, (t20 + t81) * qJD(4) + t65, -0.2e1 * t18 * qJD(4) + t73, -t17 - t98, t6 * t18 + t5 * t20 + t39 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t17 + t98, t73, -t94, (t20 - t81) * qJD(4) - t65, qJDD(4), g(3) * t41 - t102 * t42 - t31 * t20 + t79, g(3) * t42 + t102 * t41 + t31 * t18 - t70, 0, 0;];
tau_reg = t11;
