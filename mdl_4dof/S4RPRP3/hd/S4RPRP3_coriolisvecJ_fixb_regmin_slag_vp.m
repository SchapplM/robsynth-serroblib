% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tauc_reg [4x13]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:48
% EndTime: 2019-12-31 16:42:48
% DurationCPUTime: 0.19s
% Computational Cost: add. (165->44), mult. (404->74), div. (0->0), fcn. (198->4), ass. (0->37)
t13 = sin(pkin(6)) * pkin(1) + pkin(5);
t34 = qJ(4) + t13;
t42 = 2 * qJD(3);
t20 = sin(qJ(3));
t21 = cos(qJ(3));
t27 = t34 * qJD(1);
t4 = t21 * qJD(2) - t27 * t20;
t5 = t20 * qJD(2) + t27 * t21;
t35 = qJD(3) * pkin(3);
t3 = t4 + t35;
t41 = t3 - t4;
t40 = t21 * t5;
t39 = t20 * t21;
t22 = qJD(3) ^ 2;
t38 = t22 * t20;
t37 = t22 * t21;
t16 = t20 ^ 2;
t17 = t21 ^ 2;
t36 = t16 - t17;
t14 = -cos(pkin(6)) * pkin(1) - pkin(2);
t12 = qJD(1) * t14;
t33 = qJD(1) * t21;
t32 = t20 * qJD(1);
t30 = qJD(1) * qJD(4);
t29 = qJD(1) * t42;
t28 = qJD(3) * t34;
t25 = t12 * t42;
t24 = (-t21 * pkin(3) + t14) * qJD(1);
t23 = qJD(1) ^ 2;
t10 = t34 * t21;
t9 = t34 * t20;
t8 = qJD(4) + t24;
t7 = -t20 * qJD(4) - t21 * t28;
t6 = t21 * qJD(4) - t20 * t28;
t2 = -qJD(3) * t5 - t20 * t30;
t1 = t4 * qJD(3) + t21 * t30;
t11 = [0, 0, 0, 0, t29 * t39, -t36 * t29, t37, -t38, 0, -t13 * t37 + t20 * t25, t13 * t38 + t21 * t25, t1 * t21 - t2 * t20 + (-t20 * t5 - t21 * t3) * qJD(3) + (-t20 * t7 + t21 * t6 + (-t10 * t20 + t21 * t9) * qJD(3)) * qJD(1), t1 * t10 - t2 * t9 + t3 * t7 + t5 * t6 + (t8 + t24) * t20 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t37, 0, t1 * t20 + t2 * t21 + (-t20 * t3 + t40) * qJD(3); 0, 0, 0, 0, -t23 * t39, t36 * t23, 0, 0, 0, -t12 * t32, -t12 * t33, (-t35 + t41) * t33, t41 * t5 + (-t8 * t32 + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t16 - t17) * t23, (-t40 + (t3 + t35) * t20) * qJD(1);];
tauc_reg = t11;
