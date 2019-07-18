% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:27
% EndTime: 2019-07-18 13:27:27
% DurationCPUTime: 0.29s
% Computational Cost: add. (294->81), mult. (520->109), div. (0->0), fcn. (278->10), ass. (0->61)
t31 = qJD(2) + qJD(3);
t38 = cos(qJ(3));
t62 = pkin(1) * qJD(2);
t13 = pkin(2) * t31 + t38 * t62;
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t35 = sin(qJ(3));
t57 = qJDD(2) * t35;
t61 = qJD(3) * t38;
t44 = (qJD(2) * t61 + t57) * pkin(1);
t65 = t38 * pkin(1);
t24 = qJDD(2) * t65;
t30 = qJDD(2) + qJDD(3);
t55 = t35 * t62;
t8 = pkin(2) * t30 - qJD(3) * t55 + t24;
t71 = t34 * t8 + (qJD(4) * t13 + t44) * t37;
t33 = qJ(2) + qJ(3);
t29 = qJ(4) + t33;
t21 = sin(t29);
t22 = cos(t29);
t69 = g(1) * t22 + g(3) * t21;
t27 = sin(t33);
t28 = cos(t33);
t68 = g(1) * t28 + g(3) * t27;
t25 = qJDD(4) + t30;
t67 = pkin(2) * t25;
t64 = t34 * t35;
t63 = t35 * t37;
t59 = qJD(4) * t34;
t58 = qJD(4) * t37;
t56 = t24 + t68;
t54 = t35 * t58;
t53 = -g(1) * t27 + g(3) * t28;
t52 = qJD(2) * (-qJD(3) + t31);
t51 = qJD(3) * (-qJD(2) - t31);
t50 = qJD(4) * t55;
t14 = t34 * t50;
t49 = -g(1) * t21 + g(3) * t22 + t14;
t36 = sin(qJ(2));
t39 = cos(qJ(2));
t48 = g(1) * t39 + g(3) * t36;
t47 = -t34 * t38 - t63;
t46 = t37 * t38 - t64;
t26 = qJD(4) + t31;
t43 = (-pkin(2) * t26 - t13) * qJD(4) - t44;
t7 = t37 * t8;
t2 = -t13 * t59 + t7 + (-t34 * t57 + (-t34 * t61 - t54) * qJD(2)) * pkin(1);
t42 = t49 - t71;
t41 = t2 + t69;
t32 = qJDD(1) + g(2);
t23 = pkin(2) + t65;
t12 = pkin(1) * t63 + t23 * t34;
t11 = -pkin(1) * t64 + t23 * t37;
t10 = t46 * t62;
t9 = t47 * t62;
t6 = t13 * t34 + t37 * t55;
t5 = t13 * t37 - t34 * t55;
t4 = -t23 * t59 + (qJD(3) * t47 - t54) * pkin(1);
t3 = t23 * t58 + (qJD(3) * t46 - t35 * t59) * pkin(1);
t1 = -t14 + t71;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t48, -g(1) * t36 + g(3) * t39, 0, 0, 0, 0, 0, 0, 0, t30, (t30 * t38 + t35 * t51) * pkin(1) + t56, ((-qJDD(2) - t30) * t35 + t38 * t51) * pkin(1) + t53, 0, (t48 + (t35 ^ 2 + t38 ^ 2) * qJDD(2) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t25, t11 * t25 + t4 * t26 + t41, -t12 * t25 - t3 * t26 + t42, 0, t1 * t12 + t6 * t3 + t2 * t11 + t5 * t4 - g(1) * (-pkin(1) * t39 - pkin(2) * t28) - g(3) * (-pkin(1) * t36 - pkin(2) * t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, pkin(1) * t35 * t52 + t56, (t38 * t52 - t57) * pkin(1) + t53, 0, 0, 0, 0, 0, 0, 0, t25, -t9 * t26 + t7 + (-t50 + t67) * t37 + t43 * t34 + t69, t10 * t26 + (-t8 - t67) * t34 + t43 * t37 + t49, 0, -t6 * t10 - t5 * t9 + (t1 * t34 + t2 * t37 + (-t34 * t5 + t37 * t6) * qJD(4) + t68) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t6 * t26 + t41, t5 * t26 + t42, 0, 0;];
tau_reg  = t15;
