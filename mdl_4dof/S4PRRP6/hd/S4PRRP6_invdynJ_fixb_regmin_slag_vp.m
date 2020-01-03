% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRP6
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:44
% EndTime: 2019-12-31 16:30:45
% DurationCPUTime: 0.50s
% Computational Cost: add. (360->125), mult. (806->172), div. (0->0), fcn. (502->6), ass. (0->78)
t34 = sin(pkin(6));
t35 = cos(pkin(6));
t52 = g(1) * t35 + g(2) * t34;
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t48 = t38 * pkin(3) + t36 * qJ(4) + pkin(2);
t92 = t48 * qJD(2);
t51 = pkin(3) * t36 - qJ(4) * t38;
t9 = t51 * qJD(3) - t36 * qJD(4);
t91 = qJD(2) * t9 - t48 * qJDD(2);
t37 = sin(qJ(2));
t64 = qJD(1) * qJD(2);
t90 = (t52 + t64) * t37;
t71 = qJDD(3) * pkin(3);
t25 = qJD(2) * pkin(5) + t37 * qJD(1);
t82 = t38 * t25;
t89 = qJD(3) * t82 - t71;
t86 = g(3) * t37;
t85 = t36 * t25;
t39 = cos(qJ(2));
t84 = t36 * t39;
t83 = t37 * t38;
t81 = t38 * t39;
t41 = qJD(2) ^ 2;
t80 = t41 * t39;
t32 = t36 ^ 2;
t33 = t38 ^ 2;
t79 = t32 - t33;
t78 = t32 + t33;
t40 = qJD(3) ^ 2;
t77 = t40 + t41;
t76 = qJD(2) * pkin(2);
t74 = pkin(5) * qJDD(3);
t73 = qJD(2) * t36;
t70 = t39 * qJD(1);
t69 = qJDD(1) - g(3);
t68 = qJD(3) * qJ(4);
t67 = qJDD(3) * t36;
t66 = t38 * qJDD(2);
t65 = t39 * qJDD(2);
t63 = qJD(2) * qJD(3);
t62 = qJDD(3) * qJ(4);
t61 = t36 * t41 * t38;
t60 = t36 * t63;
t30 = t37 * t64;
t26 = -t70 - t76;
t59 = t26 - t76;
t6 = -t70 - t92;
t58 = t6 - t92;
t57 = t36 * qJD(3) * t70 + t38 * t30 + t52 * t83;
t56 = t78 * qJDD(2);
t55 = -qJD(3) * pkin(3) + qJD(4);
t54 = -t39 * qJDD(1) + t30;
t53 = pkin(5) * t40 + g(3) * t39;
t10 = t55 + t85;
t15 = t68 + t82;
t50 = t10 * t36 + t15 * t38;
t12 = t34 * t81 - t35 * t36;
t14 = t34 * t36 + t35 * t81;
t17 = qJDD(2) * pkin(5) + t37 * qJDD(1) + t39 * t64;
t8 = t38 * t17;
t49 = g(1) * t14 + g(2) * t12 - t8;
t11 = t34 * t84 + t35 * t38;
t13 = -t34 * t38 + t35 * t84;
t7 = t36 * t17;
t47 = g(1) * t13 + g(2) * t11 + t36 * t86 - t7;
t45 = -qJDD(4) + t47;
t44 = 0.2e1 * qJDD(2) * pkin(2) - t53 - t54;
t1 = t54 + t91;
t43 = -t1 - t53 - t91;
t4 = t62 + t8 + (qJD(4) - t85) * qJD(3);
t5 = qJDD(4) + t7 + t89;
t42 = t5 * t36 + t4 * t38 + (t10 * t38 - t15 * t36) * qJD(3);
t31 = t36 * qJDD(2);
t19 = t51 * qJD(2);
t3 = (-0.2e1 * t60 + t66) * t39 + (-t77 * t38 - t67) * t37;
t2 = (-qJDD(3) * t37 - 0.2e1 * t39 * t63) * t38 + (t77 * t37 - t65) * t36;
t16 = [t69, 0, -t41 * t37 + t65, -qJDD(2) * t37 - t80, 0, 0, 0, 0, 0, t3, t2, t3, t37 * t56 + t78 * t80, -t2, -g(3) + (t50 * qJD(2) - t1) * t39 + (qJD(2) * t6 + t42) * t37; 0, qJDD(2), t52 * t37 + t69 * t39, -t69 * t37 + t52 * t39, t32 * qJDD(2) + 0.2e1 * t38 * t60, 0.2e1 * t36 * t66 - 0.2e1 * t79 * t63, t40 * t38 + t67, qJDD(3) * t38 - t40 * t36, 0, (t59 * qJD(3) - t74) * t36 + t44 * t38 + t57, (-t74 + (t59 + t70) * qJD(3)) * t38 + (-t44 - t90) * t36, (t58 * qJD(3) - t74) * t36 + t43 * t38 + t57, -t86 + pkin(5) * t56 + (-t78 * t64 - t52) * t39 + t42, (t74 + (-t58 - t70) * qJD(3)) * t38 + (t43 + t90) * t36, -t1 * t48 + t6 * t9 + t42 * pkin(5) + (-t52 * pkin(5) - g(3) * t48 - t50 * qJD(1)) * t39 + (-g(3) * pkin(5) - t6 * qJD(1) + t52 * t48) * t37; 0, 0, 0, 0, -t61, t79 * t41, t31, t66, qJDD(3), -t26 * t73 + t47, (-qJD(2) * t26 + t86) * t38 + t49, 0.2e1 * t71 + (t19 * t38 - t36 * t6) * qJD(2) + t45, -t51 * qJDD(2) + ((t15 - t68) * t36 + (-t10 + t55) * t38) * qJD(2), -g(3) * t83 + 0.2e1 * t62 + 0.2e1 * qJD(3) * qJD(4) + (t19 * t36 + t38 * t6) * qJD(2) - t49, t4 * qJ(4) - t5 * pkin(3) - t6 * t19 - t10 * t82 - g(1) * (-t13 * pkin(3) + t14 * qJ(4)) - g(2) * (-t11 * pkin(3) + t12 * qJ(4)) + t51 * t86 + (qJD(4) + t85) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t61, t31, -t32 * t41 - t40, -t15 * qJD(3) + t6 * t73 - t45 + t89;];
tau_reg = t16;
