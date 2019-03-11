% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRP2
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
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:53
% EndTime: 2019-03-08 18:30:54
% DurationCPUTime: 0.19s
% Computational Cost: add. (231->74), mult. (337->80), div. (0->0), fcn. (169->4), ass. (0->44)
t46 = qJD(1) - qJD(3);
t35 = -pkin(1) - pkin(2);
t16 = t35 * qJD(1) + qJD(2);
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t49 = qJ(2) * qJD(1);
t7 = t16 * t31 + t33 * t49;
t63 = t7 * t46;
t62 = -qJD(3) * t49 + t35 * qJDD(1) + qJDD(2);
t47 = qJD(1) * qJD(2);
t40 = qJD(3) * t16 + t47;
t48 = qJ(2) * qJDD(1);
t61 = -t40 - t48;
t23 = t33 * t35;
t60 = -qJ(2) * t31 + t23;
t32 = sin(qJ(1));
t34 = cos(qJ(1));
t56 = t31 * t34;
t10 = -t32 * t33 + t56;
t55 = t32 * t31;
t9 = -t34 * t33 - t55;
t59 = -g(1) * t10 + g(2) * t9;
t54 = t34 * pkin(1) + t32 * qJ(2);
t53 = g(1) * t32 - g(2) * t34;
t51 = pkin(1) * qJDD(1);
t45 = 0.2e1 * t47;
t43 = t46 ^ 2;
t42 = qJDD(2) - t51;
t41 = g(1) * t34 + g(2) * t32;
t6 = t33 * t16 - t31 * t49;
t14 = qJ(2) * t33 + t31 * t35;
t39 = t62 * t33;
t2 = t62 * t31 - t61 * t33;
t38 = -g(1) * t9 - g(2) * t10 - t2;
t37 = t61 * t31 + t39;
t36 = qJD(1) ^ 2;
t29 = qJDD(1) - qJDD(3);
t25 = t34 * qJ(2);
t22 = pkin(3) * t33 + pkin(2);
t5 = -t31 * qJD(2) - t14 * qJD(3);
t4 = t33 * qJD(2) + t60 * qJD(3);
t3 = -pkin(3) * t46 + t6;
t1 = -t29 * pkin(3) + t37;
t8 = [qJDD(1), t53, t41, -qJDD(2) + 0.2e1 * t51 + t53, -t41 + t45 + 0.2e1 * t48, -t42 * pkin(1) - g(1) * (-pkin(1) * t32 + t25) - g(2) * t54 + (t45 + t48) * qJ(2), t29, -t23 * t29 - t5 * t46 + ((qJDD(1) + t29) * qJ(2) + t40) * t31 - t39 + t59, t14 * t29 + t4 * t46 - t38, t2 * t14 + t7 * t4 + t1 * (-pkin(3) + t60) + t3 * t5 - g(1) * (pkin(3) * t56 + t25 + (-pkin(1) - t22) * t32) - g(2) * (pkin(3) * t55 + t22 * t34 + t54); 0, 0, 0, -qJDD(1), -t36, -qJ(2) * t36 + t42 - t53, 0, -t33 * t29 - t31 * t43, t31 * t29 - t33 * t43 (t1 - t63) * t33 + (t46 * t3 + t2) * t31 - t53; 0, 0, 0, 0, 0, 0, -t29, t37 - t59 - t63, -t46 * t6 + t38 (-t6 + t3) * t7 + (t1 - t59) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) + g(3);];
tau_reg  = t8;
