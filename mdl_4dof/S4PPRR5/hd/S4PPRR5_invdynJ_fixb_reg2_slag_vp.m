% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPRR5
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:54
% EndTime: 2019-12-31 16:19:55
% DurationCPUTime: 0.34s
% Computational Cost: add. (243->83), mult. (547->118), div. (0->0), fcn. (363->6), ass. (0->57)
t24 = cos(qJ(3));
t13 = t24 * qJDD(2);
t22 = sin(qJ(3));
t10 = t24 * qJD(1) + t22 * qJD(2);
t54 = t10 * qJD(3);
t4 = -t22 * qJDD(1) + t13 - t54;
t56 = qJDD(3) * pkin(3);
t2 = -t56 - t4;
t19 = sin(pkin(6));
t20 = cos(pkin(6));
t42 = -g(1) * t19 + g(2) * t20;
t33 = t42 * t24;
t66 = g(3) * t22 - t2 + t33;
t25 = qJD(4) ^ 2;
t65 = pkin(5) * t25 - t54 - t56 - t66;
t21 = sin(qJ(4));
t16 = t21 ^ 2;
t23 = cos(qJ(4));
t17 = t23 ^ 2;
t60 = t16 + t17;
t62 = g(3) * t24;
t30 = t22 * t42 - t62;
t51 = t22 * qJDD(3);
t26 = qJD(3) ^ 2;
t59 = t25 + t26;
t64 = t59 * t24 + t51;
t61 = t16 - t17;
t58 = qJD(3) * pkin(3);
t53 = t24 * qJD(2);
t9 = -t22 * qJD(1) + t53;
t57 = t9 * qJD(3);
t55 = qJDD(3) * pkin(5);
t18 = qJDD(1) - g(3);
t52 = qJDD(4) * t21;
t50 = t23 * qJDD(3);
t49 = t24 * qJDD(3);
t48 = qJD(1) * qJD(3);
t47 = qJD(3) * qJD(4);
t46 = t21 * t26 * t23;
t45 = -qJD(3) * t53 - t24 * qJDD(1) - t22 * qJDD(2);
t3 = -t22 * t48 - t45;
t1 = t3 + t55;
t43 = t60 * t1;
t41 = t22 * t60;
t40 = t60 * t24;
t39 = t21 * t47;
t38 = qJDD(3) * t60;
t36 = t23 * t39;
t35 = g(1) * t20 + g(2) * t19;
t31 = -qJDD(4) * t24 + 0.2e1 * t22 * t47;
t5 = -t9 - t58;
t29 = -pkin(5) * qJDD(4) + (t5 + t9 - t58) * qJD(4);
t27 = -t5 * qJD(3) - t1 - t30;
t8 = -t24 * t26 - t51;
t7 = t22 * t26 - t49;
t6 = qJD(3) * pkin(5) + t10;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, t8, t7, 0, -t4 * t22 + t3 * t24 - g(3) + (-t10 * t22 - t24 * t9) * qJD(3), 0, 0, 0, 0, 0, 0, t31 * t21 - t64 * t23, t64 * t21 + t31 * t23, t24 * t38 - t26 * t41, t2 * t22 - g(3) + t1 * t40 + (t24 * t5 - t6 * t41) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t42, 0, 0, 0, 0, 0, 0, -t7, t8, 0, t3 * t22 + t4 * t24 + (t10 * t24 - t22 * t9) * qJD(3) + t42, 0, 0, 0, 0, 0, 0, (-0.2e1 * t39 + t50) * t24 + (-t59 * t23 - t52) * t22, (-qJDD(4) * t22 - 0.2e1 * t24 * t47) * t23 + (t59 * t22 - t49) * t21, t22 * t38 + t26 * t40, -t2 * t24 + t22 * t43 + (t22 * t5 + t6 * t40) * qJD(3) + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t18 * t22 + t13 + t33, t62 + t57 + (-t42 + t48) * t22 + t45, 0, 0, t16 * qJDD(3) + 0.2e1 * t36, 0.2e1 * t21 * t50 - 0.2e1 * t61 * t47, t25 * t23 + t52, t17 * qJDD(3) - 0.2e1 * t36, qJDD(4) * t23 - t25 * t21, 0, t29 * t21 - t65 * t23, t65 * t21 + t29 * t23, t30 + t60 * (t1 + t55 - t57), -t5 * t10 - t60 * t9 * t6 + t66 * pkin(3) + (t43 + t30) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t61 * t26, t21 * qJDD(3), t46, t50, qJDD(4), t27 * t21 - t35 * t23, t35 * t21 + t27 * t23, 0, 0;];
tau_reg = t11;
