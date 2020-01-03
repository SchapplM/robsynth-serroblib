% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRP6
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
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:15
% EndTime: 2019-12-31 16:46:16
% DurationCPUTime: 0.30s
% Computational Cost: add. (262->95), mult. (475->121), div. (0->0), fcn. (217->4), ass. (0->59)
t29 = sin(qJ(3));
t65 = t29 * pkin(3);
t43 = qJ(2) + t65;
t10 = t43 * qJD(1) + qJD(4);
t30 = sin(qJ(1));
t32 = cos(qJ(1));
t63 = g(1) * t30 - g(2) * t32;
t71 = -t10 * qJD(1) - t63;
t35 = qJD(1) ^ 2;
t37 = -t35 * qJ(2) - t63;
t33 = -pkin(1) - pkin(5);
t14 = t33 * qJD(1) + qJD(2);
t70 = -qJ(4) * qJD(1) + t14;
t24 = qJDD(1) * qJ(2);
t41 = g(1) * t32 + g(2) * t30;
t25 = qJD(1) * qJD(2);
t45 = 0.2e1 * t25;
t69 = 0.2e1 * t24 + t45 - t41;
t59 = qJD(3) * pkin(3);
t31 = cos(qJ(3));
t7 = t70 * t31;
t3 = t7 + t59;
t67 = t3 - t7;
t66 = g(3) * t29;
t64 = t32 * pkin(1) + t30 * qJ(2);
t26 = t29 ^ 2;
t27 = t31 ^ 2;
t62 = -t26 - t27;
t61 = t26 - t27;
t34 = qJD(3) ^ 2;
t60 = -t34 - t35;
t57 = qJ(4) - t33;
t56 = pkin(1) * qJDD(1);
t12 = t57 * t31;
t55 = qJD(3) * t12;
t54 = qJD(3) * t29;
t53 = qJD(3) * t31;
t50 = qJDD(3) * t29;
t49 = t29 * qJDD(1);
t48 = t31 * qJDD(1);
t47 = qJD(1) * qJD(3);
t46 = qJD(1) * qJD(4);
t44 = t31 * t47;
t42 = qJDD(2) - t56;
t39 = qJDD(4) + t24 + t25 + (t44 + t49) * pkin(3);
t38 = 0.2e1 * qJ(2) * t47 + qJDD(3) * t33;
t36 = -t33 * t34 + t69;
t28 = -qJ(4) - pkin(5);
t20 = t32 * qJ(2);
t18 = qJDD(3) * t31;
t13 = t33 * qJDD(1) + qJDD(2);
t11 = t57 * t29;
t8 = t31 * t13;
t6 = t70 * t29;
t5 = -t29 * qJD(4) - t55;
t4 = -t31 * qJD(4) + t57 * t54;
t2 = t70 * t53 + (-qJ(4) * qJDD(1) + t13 - t46) * t29;
t1 = -t31 * t46 - t14 * t54 + qJDD(3) * pkin(3) + t8 + (t29 * t47 - t48) * qJ(4);
t9 = [qJDD(1), t63, t41, qJDD(2) - 0.2e1 * t56 - t63, t69, -t42 * pkin(1) - g(1) * (-t30 * pkin(1) + t20) - g(2) * t64 + (t45 + t24) * qJ(2), t27 * qJDD(1) - 0.2e1 * t29 * t44, -0.2e1 * t29 * t48 + 0.2e1 * t61 * t47, -t34 * t29 + t18, -t34 * t31 - t50, 0, t36 * t29 + t38 * t31, -t38 * t29 + t36 * t31, (-qJD(3) * t6 + qJDD(1) * t12 - t1 + (qJD(3) * t11 - t4) * qJD(1)) * t31 + (qJD(3) * t3 + qJDD(1) * t11 - t2 + (-t5 - t55) * qJD(1)) * t29 + t63, -t2 * t11 + t6 * t5 - t1 * t12 + t3 * t4 + t39 * t43 + t10 * (pkin(3) * t53 + qJD(2)) - g(1) * (t32 * t65 + t20 + (-pkin(1) + t28) * t30) - g(2) * (-t32 * t28 + t30 * t65 + t64); 0, 0, 0, qJDD(1), -t35, t42 + t37, 0, 0, 0, 0, 0, t60 * t29 + t18, t60 * t31 - t50, t62 * qJDD(1), t1 * t31 + t2 * t29 + (-t29 * t3 + t31 * t6) * qJD(3) + t71; 0, 0, 0, 0, 0, 0, t31 * t35 * t29, -t61 * t35, t48, -t49, qJDD(3), t37 * t31 + t66 + t8, g(3) * t31 + (-t13 - t37) * t29, -pkin(3) * t48 + (t59 - t67) * t29 * qJD(1), t67 * t6 + (t71 * t31 + t1 + t66) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t35, (t29 * t6 + t3 * t31) * qJD(1) + t39 - t41;];
tau_reg = t9;
