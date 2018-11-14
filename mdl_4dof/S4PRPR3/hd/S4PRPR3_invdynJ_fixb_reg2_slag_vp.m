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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:11:27
% EndTime: 2018-11-14 14:11:27
% DurationCPUTime: 0.26s
% Computational Cost: add. (418->90), mult. (838->125), div. (0->0), fcn. (660->10), ass. (0->58)
t52 = cos(qJ(2));
t41 = t52 * qJDD(1);
t50 = sin(qJ(2));
t60 = qJD(1) * qJD(2);
t64 = qJDD(2) * pkin(2);
t28 = -t50 * t60 + t41 + t64;
t48 = cos(pkin(6));
t21 = t48 * t28;
t47 = sin(pkin(6));
t61 = t50 * qJDD(1);
t56 = t52 * t60 + t61;
t11 = -t56 * t47 + t21;
t10 = qJDD(2) * pkin(3) + t11;
t12 = t47 * t28 + t56 * t48;
t63 = t52 * qJD(1);
t32 = qJD(2) * pkin(2) + t63;
t65 = qJD(1) * t50;
t17 = t48 * t32 - t47 * t65;
t16 = qJD(2) * pkin(3) + t17;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t18 = t47 * t32 + t48 * t65;
t66 = t49 * t18;
t1 = (qJD(4) * t16 + t12) * t51 - qJD(4) * t66 + t49 * t10;
t69 = pkin(2) * t47;
t38 = t48 * pkin(2) + pkin(3);
t22 = t51 * t38 - t49 * t69;
t30 = t47 * t52 + t48 * t50;
t25 = t30 * qJD(1);
t29 = -t47 * t50 + t48 * t52;
t27 = t29 * qJD(1);
t68 = t22 * qJD(4) + t49 * t25 - t51 * t27;
t23 = t49 * t38 + t51 * t69;
t67 = -t23 * qJD(4) + t51 * t25 + t49 * t27;
t62 = qJDD(1) - g(1);
t45 = qJ(2) + pkin(6);
t57 = -g(1) * t52 + g(2) * t50;
t6 = t49 * t16 + t51 * t18;
t13 = t51 * t29 - t49 * t30;
t14 = t49 * t29 + t51 * t30;
t2 = -t6 * qJD(4) + t51 * t10 - t49 * t12;
t42 = qJ(4) + t45;
t36 = sin(t42);
t37 = cos(t42);
t55 = g(1) * t36 + g(2) * t37 - t1;
t54 = -g(1) * t37 + g(2) * t36 + t2;
t53 = qJD(2) ^ 2;
t46 = qJDD(3) + g(3);
t44 = qJD(2) + qJD(4);
t43 = qJDD(2) + qJDD(4);
t40 = cos(t45);
t39 = sin(t45);
t26 = t29 * qJD(2);
t24 = t30 * qJD(2);
t5 = t51 * t16 - t66;
t4 = -t14 * qJD(4) - t51 * t24 - t49 * t26;
t3 = t13 * qJD(4) - t49 * t24 + t51 * t26;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, t52 * qJDD(2) - t53 * t50, -qJDD(2) * t50 - t53 * t52, 0, -g(1) + (t50 ^ 2 + t52 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t24 * qJD(2) + t29 * qJDD(2), -t26 * qJD(2) - t30 * qJDD(2), 0, t11 * t29 + t12 * t30 - t17 * t24 + t18 * t26 - g(1), 0, 0, 0, 0, 0, 0, t13 * t43 + t4 * t44, -t14 * t43 - t3 * t44, 0, t1 * t14 + t2 * t13 + t6 * t3 + t5 * t4 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t41 + t57, g(2) * t52 - t62 * t50, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t48 * t64 - t47 * t61 - g(1) * t40 + g(2) * t39 + t21 + (-t47 * t63 + t25) * qJD(2), -t48 * t61 + g(1) * t39 + g(2) * t40 + (-t28 - t64) * t47 + (-t48 * t63 + t27) * qJD(2), 0, t17 * t25 - t18 * t27 + (t11 * t48 + t12 * t47 + t57) * pkin(2), 0, 0, 0, 0, 0, t43, t22 * t43 + t67 * t44 + t54, -t23 * t43 - t68 * t44 + t55, 0, t1 * t23 + t2 * t22 - g(1) * (t52 * pkin(2) + pkin(3) * t40) - g(2) * (-t50 * pkin(2) - pkin(3) * t39) + t68 * t6 + t67 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t6 * t44 + t54, t5 * t44 + t55, 0, 0;];
tau_reg  = t7;
