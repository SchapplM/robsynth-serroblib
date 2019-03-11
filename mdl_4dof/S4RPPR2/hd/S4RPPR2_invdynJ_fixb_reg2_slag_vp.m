% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPPR2
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
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:34
% EndTime: 2019-03-08 18:28:34
% DurationCPUTime: 0.32s
% Computational Cost: add. (486->102), mult. (735->124), div. (0->0), fcn. (424->8), ass. (0->64)
t54 = sin(qJ(1));
t56 = cos(qJ(1));
t70 = pkin(6) + qJ(4);
t63 = sin(t70);
t64 = cos(t70);
t19 = -t54 * t63 - t56 * t64;
t57 = -pkin(1) - pkin(2);
t35 = t57 * qJDD(1) + qJDD(2);
t52 = cos(pkin(6));
t31 = t52 * t35;
t51 = sin(pkin(6));
t75 = t51 * qJ(2);
t65 = -pkin(3) - t75;
t71 = qJD(1) * qJD(2);
t67 = t51 * t71;
t12 = t65 * qJDD(1) + t31 - t67;
t66 = t52 * t71;
t72 = qJ(2) * qJDD(1);
t78 = t51 * t35 + t52 * t72;
t15 = t66 + t78;
t53 = sin(qJ(4));
t55 = cos(qJ(4));
t39 = t57 * qJD(1) + qJD(2);
t34 = t52 * t39;
t16 = t65 * qJD(1) + t34;
t73 = qJ(2) * qJD(1);
t18 = t51 * t39 + t52 * t73;
t6 = t53 * t16 + t55 * t18;
t2 = -t6 * qJD(4) + t55 * t12 - t53 * t15;
t20 = -t54 * t64 + t56 * t63;
t89 = g(1) * t20 - g(2) * t19 + t2;
t5 = t55 * t16 - t53 * t18;
t88 = qJD(1) - qJD(4);
t83 = t54 * t51;
t81 = t56 * t51;
t27 = t55 * t51 + t53 * t52;
t80 = t88 * t27;
t25 = -t53 * t51 + t55 * t52;
t79 = t88 * t25;
t77 = t56 * pkin(1) + t54 * qJ(2);
t76 = g(1) * t54 - g(2) * t56;
t74 = pkin(1) * qJDD(1);
t69 = 0.2e1 * t71;
t32 = t52 * t57 - t75;
t62 = qJDD(2) - t74;
t61 = g(1) * t56 + g(2) * t54;
t60 = (-t51 * t73 + t34) * t51 - t18 * t52;
t29 = -pkin(3) + t32;
t33 = t52 * qJ(2) + t51 * t57;
t9 = t55 * t29 - t53 * t33;
t10 = t53 * t29 + t55 * t33;
t1 = t5 * qJD(4) + t53 * t12 + t55 * t15;
t59 = -g(1) * t19 - g(2) * t20 - t1;
t58 = qJD(1) ^ 2;
t50 = qJDD(3) + g(3);
t48 = qJDD(1) - qJDD(4);
t44 = t56 * qJ(2);
t41 = t52 * pkin(3) + pkin(2);
t28 = t56 * t52 + t83;
t26 = -t54 * t52 + t81;
t14 = t31 + (-t71 - t72) * t51;
t4 = -t27 * qJD(2) - t10 * qJD(4);
t3 = t25 * qJD(2) + t9 * qJD(4);
t7 = [0, 0, 0, 0, 0, qJDD(1), t76, t61, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t74 + t76, 0, -t61 + t69 + 0.2e1 * t72, -t62 * pkin(1) - g(1) * (-t54 * pkin(1) + t44) - g(2) * t77 + (t69 + t72) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t67 - g(1) * t26 - g(2) * t28 - t31 + (-t32 + t75) * qJDD(1), -g(1) * t28 + g(2) * t26 + t33 * qJDD(1) + 0.2e1 * t66 + t78, 0, t15 * t33 + t14 * t32 - g(1) * (t57 * t54 + t44) - g(2) * (t56 * pkin(2) + t77) - t60 * qJD(2), 0, 0, 0, 0, 0, t48, -t4 * t88 - t9 * t48 - t89, t10 * t48 + t3 * t88 - t59, 0, t1 * t10 + t6 * t3 + t2 * t9 + t5 * t4 - g(1) * (pkin(3) * t81 + t44 + (-pkin(1) - t41) * t54) - g(2) * (pkin(3) * t83 + t56 * t41 + t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t58, -t58 * qJ(2) + t62 - t76, 0, 0, 0, 0, 0, 0, -t52 * qJDD(1) - t51 * t58, t51 * qJDD(1) - t52 * t58, 0, t60 * qJD(1) + t14 * t52 + t15 * t51 - t76, 0, 0, 0, 0, 0, 0, -t25 * t48 - t80 * t88, t27 * t48 - t79 * t88, 0, t1 * t27 + t2 * t25 + t80 * t5 - t79 * t6 - t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t6 * t88 + t89, -t5 * t88 + t59, 0, 0;];
tau_reg  = t7;
