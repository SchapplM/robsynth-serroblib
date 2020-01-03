% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRPR5
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [4x12]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:19
% EndTime: 2019-12-31 16:23:20
% DurationCPUTime: 0.24s
% Computational Cost: add. (192->62), mult. (426->100), div. (0->0), fcn. (328->10), ass. (0->51)
t37 = sin(qJ(2));
t39 = cos(qJ(2));
t53 = qJD(1) * qJD(2);
t64 = t37 * qJDD(1) + t39 * t53;
t34 = cos(pkin(7));
t59 = qJD(1) * t37;
t23 = t34 * t59;
t32 = sin(pkin(7));
t58 = qJD(1) * t39;
t11 = t32 * t58 + t23;
t24 = pkin(2) * t32 + pkin(5);
t25 = -pkin(2) * t34 - pkin(3);
t29 = qJ(2) + pkin(7);
t26 = sin(t29);
t27 = cos(t29);
t28 = t39 * qJDD(1);
t14 = qJDD(2) * pkin(2) - t37 * t53 + t28;
t3 = t34 * t14 - t64 * t32;
t40 = qJD(4) ^ 2;
t33 = sin(pkin(6));
t35 = cos(pkin(6));
t50 = g(1) * t35 + g(2) * t33;
t63 = t50 * t26 - g(3) * t27 + qJD(2) * t11 - t24 * t40 + t3 + (pkin(3) - t25) * qJDD(2);
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t61 = t36 * t38;
t30 = t36 ^ 2;
t60 = -t38 ^ 2 + t30;
t15 = t32 * t37 - t34 * t39;
t57 = qJD(2) * t15;
t56 = qJDD(1) - g(3);
t54 = t38 * qJDD(2);
t52 = qJD(2) * qJD(4);
t4 = t32 * t14 + t64 * t34;
t16 = t32 * t39 + t34 * t37;
t49 = g(1) * t33 - g(2) * t35 - qJDD(3);
t21 = qJD(2) * pkin(2) + t58;
t7 = t21 * t34 - t32 * t59;
t10 = t16 * qJD(2);
t47 = qJD(2) * t10 + qJDD(2) * t15 + t16 * t40;
t46 = 0.2e1 * t57 * qJD(4) - qJDD(4) * t16;
t45 = -g(3) * t39 + t50 * t37;
t13 = t15 * qJD(1);
t5 = -qJD(2) * pkin(3) - t7;
t44 = -qJDD(4) * t24 + (qJD(2) * t25 - t13 + t5) * qJD(4);
t42 = -qJDD(2) * pkin(5) + g(3) * t26 - t5 * qJD(2) + t50 * t27 - t4;
t41 = qJD(2) ^ 2;
t19 = qJDD(4) * t38 - t36 * t40;
t18 = qJDD(4) * t36 + t38 * t40;
t8 = t32 * t21 + t23;
t1 = [t56, 0, qJDD(2) * t39 - t37 * t41, -qJDD(2) * t37 - t39 * t41, -t10 * t7 - t15 * t3 + t16 * t4 - t57 * t8 - g(3), 0, 0, 0, 0, 0, t46 * t36 - t47 * t38, t47 * t36 + t46 * t38; 0, qJDD(2), t28 + t45, -t56 * t37 + t50 * t39, t7 * t11 + t8 * t13 + (t3 * t34 + t32 * t4 + t45) * pkin(2), qJDD(2) * t30 + 0.2e1 * t52 * t61, 0.2e1 * t36 * t54 - 0.2e1 * t60 * t52, t18, t19, 0, t44 * t36 + t63 * t38, -t63 * t36 + t44 * t38; 0, 0, 0, 0, -t49, 0, 0, 0, 0, 0, t19, -t18; 0, 0, 0, 0, 0, -t41 * t61, t60 * t41, t36 * qJDD(2), t54, qJDD(4), t42 * t36 - t49 * t38, t49 * t36 + t42 * t38;];
tau_reg = t1;
