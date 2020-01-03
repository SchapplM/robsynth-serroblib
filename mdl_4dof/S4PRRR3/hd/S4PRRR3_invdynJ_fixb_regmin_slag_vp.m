% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRR3
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
% tau_reg [4x14]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:40
% EndTime: 2019-12-31 16:31:40
% DurationCPUTime: 0.19s
% Computational Cost: add. (216->60), mult. (317->86), div. (0->0), fcn. (168->8), ass. (0->50)
t27 = pkin(7) + qJ(2);
t24 = qJ(3) + t27;
t18 = cos(t24);
t34 = cos(qJ(3));
t57 = t34 * pkin(2);
t63 = -g(2) * t18 + qJDD(2) * t57;
t32 = sin(qJ(3));
t52 = pkin(2) * qJD(2);
t46 = t32 * t52;
t26 = qJDD(2) + qJDD(3);
t60 = t26 * pkin(3);
t62 = qJD(3) * t46 - t60 - t63;
t17 = sin(t24);
t54 = g(1) * t18 + g(2) * t17;
t15 = g(1) * t17;
t28 = qJD(2) + qJD(3);
t59 = t28 * pkin(3);
t58 = t32 * pkin(2);
t31 = sin(qJ(4));
t33 = cos(qJ(4));
t45 = t34 * t52;
t9 = -t45 - t59;
t56 = t9 * qJD(4) * t31 + t33 * t15;
t55 = t33 * t26;
t29 = t31 ^ 2;
t53 = -t33 ^ 2 + t29;
t51 = qJD(3) * t34;
t50 = t33 * qJD(4);
t49 = qJDD(1) - g(3);
t48 = qJDD(2) * t32;
t47 = t62 * t31 + t9 * t50;
t43 = qJD(2) * (-qJD(3) + t28);
t42 = qJD(3) * (-qJD(2) - t28);
t41 = t15 + t63;
t40 = -t9 * t28 - t26 * pkin(6) - (qJD(2) * t51 + t48) * pkin(2) + t54;
t35 = qJD(4) ^ 2;
t39 = pkin(6) * t35 - t28 * t46 - t60;
t19 = pkin(6) + t58;
t20 = -pkin(3) - t57;
t38 = qJD(3) * t28 * t58 + t19 * t35 + t20 * t26;
t37 = -pkin(6) * qJDD(4) + (t45 - t59) * qJD(4);
t36 = -qJDD(4) * t19 + (-pkin(2) * t51 + t20 * t28) * qJD(4);
t25 = t28 ^ 2;
t23 = cos(t27);
t22 = sin(t27);
t11 = qJDD(4) * t33 - t35 * t31;
t10 = qJDD(4) * t31 + t35 * t33;
t5 = 0.2e1 * t31 * t28 * t50 + t29 * t26;
t1 = -0.2e1 * t53 * t28 * qJD(4) + 0.2e1 * t31 * t55;
t2 = [t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10; 0, qJDD(2), g(1) * t22 - g(2) * t23, g(1) * t23 + g(2) * t22, t26, (t26 * t34 + t32 * t42) * pkin(2) + t41, ((-qJDD(2) - t26) * t32 + t34 * t42) * pkin(2) + t54, t5, t1, t10, t11, 0, t36 * t31 + (-t38 - t62) * t33 + t56, t36 * t33 + (t38 - t15) * t31 + t47; 0, 0, 0, 0, t26, t43 * t58 + t41, (t34 * t43 - t48) * pkin(2) + t54, t5, t1, t10, t11, 0, t37 * t31 + (-t39 - t62) * t33 + t56, t37 * t33 + (t39 - t15) * t31 + t47; 0, 0, 0, 0, 0, 0, 0, -t31 * t25 * t33, t53 * t25, t31 * t26, t55, qJDD(4), t40 * t31 + t49 * t33, -t49 * t31 + t40 * t33;];
tau_reg = t2;
