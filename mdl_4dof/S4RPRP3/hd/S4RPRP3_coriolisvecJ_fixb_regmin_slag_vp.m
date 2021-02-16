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
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
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
% StartTime: 2021-01-15 10:20:36
% EndTime: 2021-01-15 10:20:37
% DurationCPUTime: 0.24s
% Computational Cost: add. (215->62), mult. (524->97), div. (0->0), fcn. (252->4), ass. (0->47)
t50 = 2 * qJD(3);
t15 = sin(pkin(6)) * pkin(1) + pkin(5);
t13 = t15 * qJD(1);
t22 = sin(qJ(3));
t35 = t22 * qJD(3);
t11 = t13 * t35;
t23 = cos(qJ(3));
t29 = qJ(4) * qJD(1) + t13;
t5 = t22 * qJD(2) + t29 * t23;
t17 = t23 * qJD(2);
t4 = -t29 * t22 + t17;
t41 = qJD(3) * pkin(3);
t3 = t4 + t41;
t49 = t3 - t4;
t18 = t22 ^ 2;
t48 = pkin(3) * t18;
t47 = t23 * pkin(3);
t46 = t23 * t5;
t25 = qJD(1) ^ 2;
t45 = t23 * t25;
t24 = qJD(3) ^ 2;
t44 = t24 * t22;
t43 = t24 * t23;
t19 = t23 ^ 2;
t42 = t18 - t19;
t40 = qJ(4) + t15;
t16 = -cos(pkin(6)) * pkin(1) - pkin(2);
t12 = t16 - t47;
t38 = qJD(1) * t12;
t8 = qJD(4) + t38;
t39 = -qJD(4) - t8;
t14 = qJD(1) * t16;
t37 = t14 * qJD(1);
t34 = t22 * qJD(4);
t33 = t23 * qJD(4);
t32 = qJD(1) * t50;
t31 = qJ(4) * t35;
t30 = qJD(3) * t40;
t28 = t23 * t32;
t27 = t14 * t50;
t10 = t40 * t23;
t9 = t40 * t22;
t7 = -t23 * t30 - t34;
t6 = -t22 * t30 + t33;
t2 = -qJD(1) * t34 - t5 * qJD(3);
t1 = qJD(3) * t17 - t11 + (-t31 + t33) * qJD(1);
t20 = [0, 0, 0, 0, t22 * t28, -t42 * t32, t43, -t44, 0, -t15 * t43 + t22 * t27, t15 * t44 + t23 * t27, (t7 + (t8 + (t12 - 0.2e1 * t47) * qJD(1)) * t22) * qJD(3), (t23 * t8 - t6 + (t12 * t23 + 0.2e1 * t48) * qJD(1)) * qJD(3), t1 * t23 - t2 * t22 + (-t22 * t5 - t23 * t3) * qJD(3) + (-t22 * t7 + t23 * t6 + (-t10 * t22 + t23 * t9) * qJD(3)) * qJD(1), t1 * t10 - t2 * t9 + t3 * t7 + t5 * t6 + (t8 + t38) * pkin(3) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t43, -t44, -t43, 0, t1 * t22 + t2 * t23 + (-t22 * t3 + t46) * qJD(3); 0, 0, 0, 0, -t22 * t45, t42 * t25, 0, 0, 0, -t22 * t37, -t23 * t37, (pkin(3) * t45 + t39 * qJD(1)) * t22, -t25 * t48 + t11 + (t4 - t17) * qJD(3) + (t39 * t23 + t31) * qJD(1), (-t41 + t49) * t23 * qJD(1), t49 * t5 + (-t8 * t22 * qJD(1) + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t32, t28, (-t18 - t19) * t25, (-t46 + (t3 + t41) * t22) * qJD(1);];
tauc_reg = t20;
