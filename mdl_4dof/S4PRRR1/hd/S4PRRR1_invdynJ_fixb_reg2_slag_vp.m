% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRR1
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:22
% EndTime: 2019-03-08 18:25:22
% DurationCPUTime: 0.29s
% Computational Cost: add. (326->82), mult. (520->109), div. (0->0), fcn. (278->10), ass. (0->61)
t35 = qJD(2) + qJD(3);
t40 = cos(qJ(3));
t62 = pkin(2) * qJD(2);
t13 = t35 * pkin(3) + t40 * t62;
t37 = sin(qJ(4));
t38 = sin(qJ(3));
t56 = t38 * t62;
t51 = qJD(4) * t56;
t14 = t37 * t51;
t39 = cos(qJ(4));
t57 = qJDD(2) * t38;
t61 = qJD(3) * t40;
t45 = (qJD(2) * t61 + t57) * pkin(2);
t66 = t40 * pkin(2);
t27 = qJDD(2) * t66;
t33 = qJDD(2) + qJDD(3);
t8 = t33 * pkin(3) - qJD(3) * t56 + t27;
t1 = t37 * t8 + (qJD(4) * t13 + t45) * t39 - t14;
t34 = pkin(7) + qJ(2);
t32 = qJ(3) + t34;
t25 = qJ(4) + t32;
t21 = sin(t25);
t22 = cos(t25);
t72 = g(1) * t22 + g(2) * t21;
t71 = g(1) * t21 - g(2) * t22;
t23 = sin(t32);
t24 = cos(t32);
t70 = g(1) * t23 - g(2) * t24;
t28 = qJDD(4) + t33;
t69 = pkin(3) * t28;
t65 = t37 * t38;
t64 = t38 * t39;
t63 = g(1) * t24 + g(2) * t23;
t59 = qJD(4) * t37;
t58 = qJD(4) * t39;
t55 = t38 * t58;
t53 = qJD(2) * (-qJD(3) + t35);
t52 = qJD(3) * (-qJD(2) - t35);
t50 = t27 + t70;
t29 = sin(t34);
t30 = cos(t34);
t49 = g(1) * t29 - g(2) * t30;
t48 = -t37 * t40 - t64;
t47 = t39 * t40 - t65;
t31 = qJD(4) + t35;
t44 = (-pkin(3) * t31 - t13) * qJD(4) - t45;
t43 = -t1 + t72;
t7 = t39 * t8;
t2 = -t13 * t59 + t7 + (-t37 * t57 + (-t37 * t61 - t55) * qJD(2)) * pkin(2);
t42 = t2 + t71;
t36 = qJDD(1) - g(3);
t26 = pkin(3) + t66;
t12 = pkin(2) * t64 + t37 * t26;
t11 = -pkin(2) * t65 + t39 * t26;
t10 = t47 * t62;
t9 = t48 * t62;
t6 = t37 * t13 + t39 * t56;
t5 = t39 * t13 - t37 * t56;
t4 = -t26 * t59 + (t48 * qJD(3) - t55) * pkin(2);
t3 = t26 * t58 + (t47 * qJD(3) - t38 * t59) * pkin(2);
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t49, g(1) * t30 + g(2) * t29, 0, 0, 0, 0, 0, 0, 0, t33 (t33 * t40 + t38 * t52) * pkin(2) + t50 ((-qJDD(2) - t33) * t38 + t40 * t52) * pkin(2) + t63, 0 (t49 + (t38 ^ 2 + t40 ^ 2) * qJDD(2) * pkin(2)) * pkin(2), 0, 0, 0, 0, 0, t28, t11 * t28 + t4 * t31 + t42, -t12 * t28 - t3 * t31 + t43, 0, t1 * t12 + t6 * t3 + t2 * t11 + t5 * t4 - g(1) * (-pkin(2) * t29 - pkin(3) * t23) - g(2) * (pkin(2) * t30 + pkin(3) * t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t38 * pkin(2) * t53 + t50 (t40 * t53 - t57) * pkin(2) + t63, 0, 0, 0, 0, 0, 0, 0, t28, -t9 * t31 + t7 + (-t51 + t69) * t39 + t44 * t37 + t71, t10 * t31 + t14 + (-t8 - t69) * t37 + t44 * t39 + t72, 0, -t6 * t10 - t5 * t9 + (t1 * t37 + t2 * t39 + (-t37 * t5 + t39 * t6) * qJD(4) + t70) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t6 * t31 + t42, t5 * t31 + t43, 0, 0;];
tau_reg  = t15;
