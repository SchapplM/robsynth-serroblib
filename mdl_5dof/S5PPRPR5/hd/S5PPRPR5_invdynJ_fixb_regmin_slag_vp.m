% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% tau_reg [5x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:31
% EndTime: 2019-12-31 17:33:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (176->70), mult. (335->88), div. (0->0), fcn. (235->6), ass. (0->52)
t51 = qJD(3) * qJ(4);
t23 = sin(qJ(3));
t53 = t23 * qJD(2);
t12 = t51 + t53;
t26 = -pkin(3) - pkin(6);
t61 = (t12 + t51 - t53) * qJD(5) + qJDD(5) * t26;
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t60 = t22 * t24;
t20 = t24 ^ 2;
t59 = t22 ^ 2 - t20;
t27 = qJD(5) ^ 2;
t28 = qJD(3) ^ 2;
t58 = t27 + t28;
t57 = cos(pkin(7));
t56 = sin(pkin(7));
t55 = qJDD(3) * pkin(3);
t54 = t12 * qJD(3);
t25 = cos(qJ(3));
t52 = t25 * qJD(2);
t21 = qJDD(1) - g(3);
t50 = qJDD(3) * t23;
t49 = qJDD(5) * t22;
t47 = t23 * qJDD(2);
t46 = t24 * qJDD(3);
t45 = t25 * qJDD(2);
t44 = qJD(3) * qJD(5);
t43 = qJDD(3) * qJ(4);
t42 = qJD(4) - t52;
t4 = -t56 * t23 - t57 * t25;
t5 = t57 * t23 - t56 * t25;
t41 = g(1) * t5 - g(2) * t4;
t40 = g(1) * t4 + g(2) * t5;
t39 = -g(1) * t56 + g(2) * t57;
t38 = (-qJD(3) * pkin(3) + t42) * t23 + t12 * t25;
t16 = qJD(3) * t53;
t37 = qJDD(4) + t16 - t45;
t35 = t58 * t25 + t50;
t34 = -qJDD(5) * t25 + 0.2e1 * t23 * t44;
t33 = -t40 - t47;
t32 = t41 + t45;
t31 = -t26 * qJDD(3) - t37 + t41 + t54;
t30 = qJDD(4) - t32;
t2 = t43 + t47 + (qJD(4) + t52) * qJD(3);
t29 = t42 * qJD(3) - t26 * t27 + t2 + t40 + t43;
t17 = qJDD(5) * t24;
t10 = -t25 * qJDD(3) + t28 * t23;
t9 = t28 * t25 + t50;
t8 = -t27 * t22 + t17;
t7 = t27 * t24 + t49;
t3 = t37 - t55;
t1 = [t21, t21, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, t7, t8; 0, qJDD(2) + t39, 0, -t10, -t9, t10, t9, t38 * qJD(3) + t2 * t23 - t3 * t25 + t39, 0, 0, 0, 0, 0, t35 * t22 + t34 * t24, -t34 * t22 + t35 * t24; 0, 0, qJDD(3), t32, t33, t30 - 0.2e1 * t55, 0.2e1 * qJD(3) * qJD(4) - t33 + 0.2e1 * t43, t2 * qJ(4) + t12 * qJD(4) - t3 * pkin(3) - g(1) * (-t5 * pkin(3) - t4 * qJ(4)) - g(2) * (t4 * pkin(3) - t5 * qJ(4)) - t38 * qJD(2), t20 * qJDD(3) - 0.2e1 * t44 * t60, -0.2e1 * t22 * t46 + 0.2e1 * t59 * t44, t8, -t7, 0, t29 * t22 + t61 * t24, -t61 * t22 + t29 * t24; 0, 0, 0, 0, 0, qJDD(3), -t28, t16 + t30 - t54 - t55, 0, 0, 0, 0, 0, -t58 * t22 + t17, -t58 * t24 - t49; 0, 0, 0, 0, 0, 0, 0, 0, t28 * t60, -t59 * t28, t46, -t22 * qJDD(3), qJDD(5), t21 * t22 - t31 * t24, t21 * t24 + t31 * t22;];
tau_reg = t1;
