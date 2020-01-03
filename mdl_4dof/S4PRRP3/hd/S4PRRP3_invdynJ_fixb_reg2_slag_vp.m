% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRP3
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:58
% EndTime: 2019-12-31 16:26:59
% DurationCPUTime: 0.43s
% Computational Cost: add. (363->131), mult. (783->165), div. (0->0), fcn. (415->4), ass. (0->83)
t48 = cos(qJ(3));
t38 = t48 * qJDD(2);
t47 = sin(qJ(3));
t68 = qJD(2) * qJD(3);
t66 = t47 * t68;
t93 = -t66 + t38;
t78 = qJD(3) * pkin(3);
t40 = t48 * qJD(1);
t81 = qJ(4) + pkin(5);
t64 = qJD(2) * t81;
t8 = -t47 * t64 + t40;
t6 = t8 + t78;
t92 = t6 - t8;
t43 = t47 ^ 2;
t91 = pkin(3) * t43;
t42 = pkin(6) + qJ(2);
t34 = sin(t42);
t90 = g(1) * t34;
t35 = cos(t42);
t89 = g(1) * t35;
t88 = g(2) * t34;
t87 = g(2) * t35;
t86 = g(3) * t48;
t85 = t48 * pkin(3);
t84 = t34 * t48;
t83 = t35 * t47;
t50 = qJD(2) ^ 2;
t82 = t48 * t50;
t44 = t48 ^ 2;
t80 = t43 - t44;
t79 = t43 + t44;
t77 = pkin(5) * qJDD(2);
t76 = qJD(2) * t48;
t75 = qJD(3) * t47;
t74 = qJDD(2) * pkin(2);
t73 = qJDD(3) * pkin(3);
t72 = t47 * qJD(1);
t32 = pkin(2) + t85;
t16 = -qJD(2) * t32 + qJD(4);
t71 = -qJD(4) - t16;
t70 = qJDD(2) * t32;
t36 = t47 * qJDD(2);
t69 = qJD(1) * qJD(3);
t67 = qJD(2) * qJD(4);
t21 = t81 * t48;
t65 = t48 * t68;
t63 = qJD(3) * t81;
t3 = t93 * pkin(5) + t47 * qJDD(1) + t48 * t69;
t7 = pkin(3) * t66 + qJDD(4) - t70;
t62 = t7 - t70;
t61 = t47 * t65;
t60 = t48 * t63;
t59 = t88 + t89;
t58 = -t87 + t90;
t10 = qJD(2) * t21 + t72;
t57 = t10 * t48 - t47 * t6;
t49 = qJD(3) ^ 2;
t56 = pkin(5) * t49 - 0.2e1 * t74;
t39 = t48 * qJDD(1);
t55 = g(1) * t83 + t47 * t88 + t39 - t86;
t15 = pkin(5) * t76 + t72;
t54 = g(2) * t84 + g(3) * t47 + t48 * t89 - t3;
t53 = -0.2e1 * pkin(2) * t68 - pkin(5) * qJDD(3);
t52 = -t81 * qJDD(2) - t69;
t14 = -t47 * qJD(2) * pkin(5) + t40;
t4 = -t47 * t69 + t39 + (-t65 - t36) * pkin(5);
t51 = t3 * t48 - t4 * t47 + (-t14 * t48 - t15 * t47) * qJD(3) - t59;
t45 = qJDD(1) - g(3);
t28 = t47 * t82;
t25 = g(1) * t84;
t24 = g(2) * t83;
t20 = t81 * t47;
t19 = t80 * t50;
t18 = qJDD(3) * t48 - t49 * t47;
t17 = qJDD(3) * t47 + t49 * t48;
t13 = t44 * qJDD(2) - 0.2e1 * t61;
t12 = t43 * qJDD(2) + 0.2e1 * t61;
t11 = -t47 * qJD(4) - t60;
t9 = t48 * qJD(4) - t47 * t63;
t5 = 0.2e1 * t47 * t38 - 0.2e1 * t80 * t68;
t2 = t93 * qJ(4) + t48 * t67 + t3;
t1 = t73 + t39 - qJD(2) * t60 + (t52 - t67) * t47;
t22 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, t18, -t17, 0, t3 * t47 + t4 * t48 - g(3) + (-t14 * t47 + t15 * t48) * qJD(3), 0, 0, 0, 0, 0, 0, t18, -t17, 0, t57 * qJD(3) + t1 * t48 + t2 * t47 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t58, t59, 0, 0, t12, t5, t17, t13, t18, 0, t25 + t53 * t47 + (-t56 - t87) * t48, t24 + t53 * t48 + (t56 - t90) * t47, t79 * t77 + t51, (t58 + t74) * pkin(2) + t51 * pkin(5), t12, t5, t17, t13, t18, 0, -t20 * qJDD(3) + t25 + (-t62 - t87) * t48 + (t11 + (t16 + (-t32 - t85) * qJD(2)) * t47) * qJD(3), -t21 * qJDD(3) + t24 + (t62 - t90) * t47 + (t16 * t48 - t9 + (-t32 * t48 + t91) * qJD(2)) * qJD(3), (-qJD(3) * t6 + qJDD(2) * t21 + t2 + (qJD(3) * t20 + t9) * qJD(2)) * t48 + (-t10 * qJD(3) + qJDD(2) * t20 - t1 + (-qJD(3) * t21 - t11) * qJD(2)) * t47 - t59, t2 * t21 + t10 * t9 - t1 * t20 + t6 * t11 - t7 * t32 + t16 * pkin(3) * t75 - g(1) * (-t34 * t32 + t35 * t81) - g(2) * (t35 * t32 + t34 * t81); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t19, t36, t28, t38, qJDD(3), (pkin(2) * t50 - t77) * t47 + t55, pkin(2) * t82 + t14 * qJD(3) + t54, 0, 0, -t28, t19, t36, t28, t38, qJDD(3), 0.2e1 * t73 + (-t48 * t64 + t10) * qJD(3) + (pkin(3) * t82 + t71 * qJD(2) + t52) * t47 + t55, -t50 * t91 - qJ(4) * t38 + t8 * qJD(3) + (qJ(4) * t75 + t71 * t48) * qJD(2) + t54, -pkin(3) * t36 + (-t78 + t92) * t76, t92 * t10 + (-t86 + t1 + (-qJD(2) * t16 + t59) * t47) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38 + 0.2e1 * t66, t36 + 0.2e1 * t65, -t79 * t50, -t57 * qJD(2) - t58 + t7;];
tau_reg = t22;
