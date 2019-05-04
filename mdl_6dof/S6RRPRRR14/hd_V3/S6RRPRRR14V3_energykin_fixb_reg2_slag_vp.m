% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR14V3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energykin_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:10:04
% EndTime: 2019-04-12 15:10:05
% DurationCPUTime: 0.18s
% Computational Cost: add. (208->41), mult. (569->100), div. (0->0), fcn. (400->8), ass. (0->39)
t34 = sin(qJ(4));
t35 = sin(qJ(2));
t46 = qJD(1) * t35;
t50 = cos(qJ(4));
t52 = t50 * qJD(2) - t34 * t46;
t37 = qJD(2) ^ 2;
t30 = t37 / 0.2e1;
t38 = qJD(1) ^ 2;
t51 = t38 / 0.2e1;
t49 = cos(qJ(5));
t48 = cos(qJ(6));
t47 = t35 ^ 2 * t38;
t16 = t52 * qJ(3);
t33 = sin(qJ(5));
t10 = -t49 * qJD(3) + t33 * t16;
t45 = t10 ^ 2 / 0.2e1;
t20 = t34 * qJD(2) + t50 * t46;
t14 = t20 * qJ(3);
t44 = t14 ^ 2 / 0.2e1;
t43 = qJD(1) * qJD(2);
t36 = cos(qJ(2));
t42 = t35 * t38 * t36;
t24 = t47 / 0.2e1;
t40 = t36 * t43;
t22 = t36 * qJD(1) - qJD(4);
t7 = t33 * t20 + t49 * t22;
t32 = sin(qJ(6));
t29 = qJD(3) ^ 2 / 0.2e1;
t26 = t35 * t43;
t25 = t36 ^ 2 * t51;
t17 = qJD(5) - t52;
t12 = t33 * qJD(3) + t49 * t16;
t9 = t49 * t20 - t33 * t22;
t6 = qJD(6) + t7;
t5 = t48 * t12 + t32 * t14;
t4 = -t32 * t12 + t48 * t14;
t3 = t32 * t17 + t48 * t9;
t1 = -t48 * t17 + t32 * t9;
t2 = [0, 0, 0, 0, 0, t51, 0, 0, 0, 0, t24, t42, t26, t25, t40, t30, 0, 0, 0, 0, t24, t26, -t42, t30, -t40, t25, qJ(3) * t42 - qJD(3) * qJD(2) (qJ(3) * qJD(2) * t36 + qJD(3) * t35) * qJD(1) (t37 + t47) * qJ(3), t29 + (t30 + t24) * qJ(3) ^ 2, t20 ^ 2 / 0.2e1, t20 * t52, -t20 * t22, t52 ^ 2 / 0.2e1, -t52 * t22, t22 ^ 2 / 0.2e1, -qJD(3) * t52 + t14 * t22, qJD(3) * t20 + t16 * t22, t14 * t20 + t16 * t52, t16 ^ 2 / 0.2e1 + t44 + t29, t9 ^ 2 / 0.2e1, -t9 * t7, t9 * t17, t7 ^ 2 / 0.2e1, -t7 * t17, t17 ^ 2 / 0.2e1, -t10 * t17 + t14 * t7, -t12 * t17 + t14 * t9, t10 * t9 - t12 * t7, t12 ^ 2 / 0.2e1 + t45 + t44, t3 ^ 2 / 0.2e1, -t3 * t1, t3 * t6, t1 ^ 2 / 0.2e1, -t1 * t6, t6 ^ 2 / 0.2e1, t10 * t1 + t4 * t6, t10 * t3 - t5 * t6, -t5 * t1 - t4 * t3, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t45;];
T_reg  = t2;
