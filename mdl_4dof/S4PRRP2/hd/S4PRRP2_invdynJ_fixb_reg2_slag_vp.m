% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRP2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:24:11
% EndTime: 2019-03-08 18:24:11
% DurationCPUTime: 0.20s
% Computational Cost: add. (350->74), mult. (650->91), div. (0->0), fcn. (463->6), ass. (0->52)
t36 = qJ(2) + qJ(3);
t32 = sin(t36);
t33 = cos(t36);
t61 = g(1) * t33 + g(2) * t32;
t60 = g(1) * t32 - g(2) * t33;
t35 = qJD(2) + qJD(3);
t40 = cos(qJ(2));
t31 = t40 * qJDD(1);
t38 = sin(qJ(2));
t53 = qJD(1) * qJD(2);
t20 = qJDD(2) * pkin(2) - t38 * t53 + t31;
t37 = sin(qJ(3));
t56 = qJD(1) * t38;
t51 = qJD(3) * t56;
t24 = t37 * t51;
t39 = cos(qJ(3));
t25 = qJD(2) * pkin(2) + t40 * qJD(1);
t46 = -t38 * qJDD(1) - t40 * t53;
t45 = qJD(3) * t25 - t46;
t6 = t37 * t20 + t45 * t39 - t24;
t34 = qJDD(2) + qJDD(3);
t59 = pkin(2) * t34;
t57 = t39 * pkin(2);
t55 = qJD(3) * t37;
t54 = qJDD(1) - g(2);
t21 = -t37 * t38 + t39 * t40;
t11 = t35 * t21;
t16 = t37 * t25 + t39 * t56;
t22 = t37 * t40 + t39 * t38;
t52 = t16 * t11 + t6 * t22 - g(2);
t19 = t21 * qJD(1);
t49 = t6 * t37 * pkin(2) + (qJD(3) * t57 - t19) * t16;
t48 = t39 * t51;
t47 = g(1) * t38 - g(2) * t40;
t15 = t39 * t25 - t37 * t56;
t44 = (-pkin(2) * t35 - t25) * qJD(3) + t46;
t17 = t39 * t20;
t7 = -t45 * t37 + t17 - t48;
t18 = t22 * qJD(1);
t43 = t18 * t35 + t44 * t37 + t17 + t60;
t42 = t16 * t35 + t60 + t7;
t41 = qJD(2) ^ 2;
t30 = t34 * pkin(3);
t29 = pkin(3) + t57;
t14 = t35 * pkin(3) + t15;
t12 = t35 * t22;
t9 = -t12 * t35 + t21 * t34;
t8 = -t11 * t35 - t22 * t34;
t4 = t30 + t7;
t2 = t15 * t35 - t6 + t61;
t1 = t19 * t35 + t24 + (-t20 - t59) * t37 + t44 * t39 + t61;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, t40 * qJDD(2) - t41 * t38, -qJDD(2) * t38 - t41 * t40, 0, -g(2) + (t38 ^ 2 + t40 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t9, t8, 0, -t15 * t12 + t7 * t21 + t52, 0, 0, 0, 0, 0, 0, t9, t8, 0, -t14 * t12 + t4 * t21 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t31 + t47, g(1) * t40 - t54 * t38, 0, 0, 0, 0, 0, 0, 0, t34 (-t51 + t59) * t39 + t43, t1, 0, t15 * t18 + (-t15 * t55 + t39 * t7 + t47) * pkin(2) + t49, 0, 0, 0, 0, 0, t34, t29 * t34 + t30 + t43 - t48, t1, 0, t4 * t29 - g(1) * (-t38 * pkin(2) - pkin(3) * t32) - g(2) * (t40 * pkin(2) + pkin(3) * t33) + (-pkin(2) * t55 + t18) * t14 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t42, t2, 0, 0, 0, 0, 0, 0, 0, t34, 0.2e1 * t30 + t42, t2, 0 (t14 - t15) * t16 + (t4 + t60) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) - g(3);];
tau_reg  = t3;
