% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:57
% EndTime: 2019-12-31 16:25:58
% DurationCPUTime: 0.32s
% Computational Cost: add. (238->93), mult. (480->119), div. (0->0), fcn. (287->6), ass. (0->67)
t76 = pkin(2) + pkin(5);
t82 = t76 * qJDD(2);
t25 = cos(qJ(2));
t59 = t25 * qJD(1);
t42 = qJD(3) - t59;
t22 = sin(qJ(4));
t18 = t22 ^ 2;
t24 = cos(qJ(4));
t19 = t24 ^ 2;
t65 = t18 + t19;
t81 = (-t76 * qJD(2) + t42) * t65;
t23 = sin(qJ(2));
t50 = qJD(1) * qJD(2);
t13 = t23 * t50;
t52 = t25 * qJDD(1);
t37 = qJDD(3) + t13 - t52;
t1 = t37 - t82;
t43 = t65 * t1;
t28 = qJD(2) ^ 2;
t80 = t25 * qJDD(2) - t28 * t23;
t20 = sin(pkin(6));
t21 = cos(pkin(6));
t40 = g(1) * t21 + g(2) * t20;
t17 = g(3) * t25;
t57 = qJD(2) * qJ(3);
t60 = t23 * qJD(1);
t8 = t57 + t60;
t62 = t8 * qJD(2);
t79 = t40 * t23 - t1 - t17 + t62;
t78 = (t57 + t8 - t60) * qJD(4) - qJDD(4) * t76;
t58 = qJDD(1) - g(3);
t77 = -t58 * t23 + t40 * t25;
t73 = t8 * t25;
t47 = qJDD(2) * qJ(3);
t2 = t47 + t23 * qJDD(1) + (qJD(3) + t59) * qJD(2);
t72 = t2 * qJ(3);
t71 = t20 * t23;
t70 = t21 * t23;
t69 = t22 * t24;
t67 = t25 * pkin(2) + t23 * qJ(3);
t66 = t18 - t19;
t27 = qJD(4) ^ 2;
t64 = t27 + t28;
t63 = qJ(3) * t25;
t61 = qJDD(2) * pkin(2);
t56 = qJDD(2) * t23;
t55 = qJDD(4) * t22;
t53 = t24 * qJDD(2);
t49 = qJD(2) * qJD(3);
t48 = qJD(2) * qJD(4);
t46 = t28 * t69;
t45 = -g(1) * t70 - g(2) * t71 + t17;
t44 = t2 * t23 - g(3);
t41 = t48 * t69;
t39 = g(1) * t20 - g(2) * t21;
t38 = (-qJD(2) * pkin(2) + t42) * t23 + t73;
t36 = t45 - t62;
t34 = -t45 + t52;
t32 = t64 * t25 + t56;
t31 = -qJDD(4) * t25 + 0.2e1 * t23 * t48;
t3 = t37 - t61;
t29 = -g(3) * t23 + t47 + t49 + t76 * t27 + t2 + (-t40 - t50) * t25;
t14 = qJDD(4) * t24;
t10 = t21 * t63;
t9 = t20 * t63;
t5 = t28 * t25 + t56;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, t80, -t5, 0, -g(3) + (t23 ^ 2 + t25 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, 0, -t80, t5, t38 * qJD(2) - t3 * t25 + t44, 0, 0, 0, 0, 0, 0, t32 * t22 + t31 * t24, -t31 * t22 + t32 * t24, t80 * t65, -t25 * t43 + (t23 * t81 + t73) * qJD(2) + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t34, t77, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, qJDD(3) - t34 - 0.2e1 * t61, 0.2e1 * t47 + 0.2e1 * t49 - t77, t72 + t8 * qJD(3) - t3 * pkin(2) - g(1) * (-pkin(2) * t70 + t10) - g(2) * (-pkin(2) * t71 + t9) - g(3) * t67 - t38 * qJD(1), t19 * qJDD(2) - 0.2e1 * t41, -0.2e1 * t22 * t53 + 0.2e1 * t66 * t48, -t27 * t22 + t14, t18 * qJDD(2) + 0.2e1 * t41, -t27 * t24 - t55, 0, t29 * t22 + t78 * t24, -t78 * t22 + t29 * t24, -t45 + t65 * (-t1 + t13 + t82), t72 - g(1) * t10 - g(2) * t9 - g(3) * (t25 * pkin(5) + t67) + t42 * t8 - t76 * t43 + (-qJD(1) * t81 + t40 * t76) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t28, t3 + t36, 0, 0, 0, 0, 0, 0, -t64 * t22 + t14, -t64 * t24 - t55, -t65 * qJDD(2), t43 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t66 * t28, t53, -t46, -t22 * qJDD(2), qJDD(4), t39 * t22 - t79 * t24, t79 * t22 + t39 * t24, 0, 0;];
tau_reg = t4;
