% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPPR1
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:33
% EndTime: 2019-03-08 18:27:33
% DurationCPUTime: 0.20s
% Computational Cost: add. (336->74), mult. (463->88), div. (0->0), fcn. (239->8), ass. (0->46)
t40 = cos(pkin(6));
t29 = -t40 * pkin(1) - pkin(2);
t22 = -pkin(3) + t29;
t16 = t22 * qJD(1) + qJD(3);
t39 = sin(pkin(6));
t25 = t39 * pkin(1) + qJ(3);
t18 = t25 * qJD(1);
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t3 = t43 * t16 - t41 * t18;
t54 = qJD(1) - qJD(4);
t66 = t54 * t3;
t15 = t22 * qJDD(1) + qJDD(3);
t55 = qJD(3) * qJD(1);
t56 = pkin(1) * qJDD(1);
t30 = t39 * t56;
t57 = qJDD(1) * qJ(3) + t30;
t17 = t55 + t57;
t4 = t41 * t16 + t43 * t18;
t2 = -qJD(4) * t4 + t43 * t15 - t41 * t17;
t65 = -t4 * t54 + t2;
t36 = qJ(1) + pkin(6);
t31 = sin(t36);
t32 = cos(t36);
t12 = -t31 * t41 - t32 * t43;
t13 = -t31 * t43 + t32 * t41;
t64 = g(1) * t13 - g(2) * t12;
t63 = t29 * qJDD(1);
t58 = g(1) * t31 - g(2) * t32;
t44 = cos(qJ(1));
t53 = t44 * pkin(1) + t32 * pkin(2) + t31 * qJ(3);
t42 = sin(qJ(1));
t52 = -t42 * pkin(1) + t32 * qJ(3);
t50 = t54 ^ 2;
t49 = g(1) * t32 + g(2) * t31;
t48 = g(1) * t42 - g(2) * t44;
t7 = t43 * t22 - t41 * t25;
t8 = t41 * t22 + t43 * t25;
t1 = qJD(4) * t3 + t41 * t15 + t43 * t17;
t47 = qJDD(3) + t63;
t46 = -g(1) * t12 - g(2) * t13 - t1;
t38 = qJDD(2) - g(3);
t34 = qJDD(1) - qJDD(4);
t6 = -t41 * qJD(3) - t8 * qJD(4);
t5 = t43 * qJD(3) + t7 * qJD(4);
t9 = [0, 0, 0, 0, 0, qJDD(1), t48, g(1) * t44 + g(2) * t42, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t40 * t56 + t58, -0.2e1 * t30 + t49, 0 (t48 + (t39 ^ 2 + t40 ^ 2) * t56) * pkin(1), 0, 0, 0, qJDD(1), 0, 0, -qJDD(3) + t58 - 0.2e1 * t63, 0, t25 * qJDD(1) - t49 + 0.2e1 * t55 + t57, t17 * t25 + t18 * qJD(3) + t47 * t29 - g(1) * (-t31 * pkin(2) + t52) - g(2) * t53, 0, 0, 0, 0, 0, t34, -t7 * t34 - t54 * t6 - t2 - t64, t8 * t34 + t5 * t54 - t46, 0, t1 * t8 + t4 * t5 + t2 * t7 + t3 * t6 - g(1) * ((-pkin(2) - pkin(3)) * t31 + t52) - g(2) * (t32 * pkin(3) + t53); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -qJD(1) ^ 2, -t18 * qJD(1) + t47 - t58, 0, 0, 0, 0, 0, 0, -t43 * t34 - t41 * t50, t41 * t34 - t43 * t50, 0, t65 * t43 + (t1 + t66) * t41 - t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t64 + t65, t46 - t66, 0, 0;];
tau_reg  = t9;
