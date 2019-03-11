% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_energykin_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:30
% EndTime: 2019-03-09 03:59:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (607->58), mult. (1255->131), div. (0->0), fcn. (801->8), ass. (0->44)
t45 = sin(pkin(10));
t46 = cos(pkin(10));
t49 = sin(qJ(3));
t50 = cos(qJ(3));
t30 = (t45 * t50 + t46 * t49) * qJD(1);
t51 = qJD(1) ^ 2;
t41 = t51 / 0.2e1;
t58 = cos(qJ(5));
t57 = cos(qJ(6));
t35 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t53 = -qJ(4) * qJD(1) + t35;
t26 = qJD(3) * pkin(3) + t53 * t50;
t28 = t53 * t49;
t17 = t45 * t26 + t46 * t28;
t13 = qJD(3) * pkin(8) + t17;
t32 = (-t45 * t49 + t46 * t50) * qJD(1);
t33 = qJD(4) + (pkin(3) * t49 + qJ(2)) * qJD(1);
t18 = t30 * pkin(4) - t32 * pkin(8) + t33;
t48 = sin(qJ(5));
t6 = t58 * t13 + t48 * t18;
t56 = t51 * qJ(2);
t55 = qJD(3) * t35;
t54 = qJD(1) * qJD(3);
t5 = -t48 * t13 + t58 * t18;
t16 = t46 * t26 - t45 * t28;
t12 = -qJD(3) * pkin(4) - t16;
t29 = qJD(5) + t30;
t47 = sin(qJ(6));
t44 = t50 ^ 2;
t43 = t49 ^ 2;
t40 = qJD(3) ^ 2 / 0.2e1;
t38 = qJ(2) ^ 2 * t41;
t37 = -pkin(1) * qJD(1) + qJD(2);
t27 = qJD(6) + t29;
t22 = t48 * qJD(3) + t58 * t32;
t20 = -t58 * qJD(3) + t48 * t32;
t10 = -t47 * t20 + t57 * t22;
t8 = t57 * t20 + t47 * t22;
t7 = t20 * pkin(5) + t12;
t4 = -t20 * pkin(9) + t6;
t3 = t29 * pkin(5) - t22 * pkin(9) + t5;
t2 = t47 * t3 + t57 * t4;
t1 = t57 * t3 - t47 * t4;
t9 = [0, 0, 0, 0, 0, t41, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, t37 * qJD(1), t56, t38 + t37 ^ 2 / 0.2e1, t44 * t41, -t50 * t51 * t49, t50 * t54, t43 * t41, -t49 * t54, t40, t49 * t56 + t50 * t55, -t49 * t55 + t50 * t56 (-t43 - t44) * t35 * qJD(1), t38 + (t43 / 0.2e1 + t44 / 0.2e1) * t35 ^ 2, t32 ^ 2 / 0.2e1, -t32 * t30, t32 * qJD(3), t30 ^ 2 / 0.2e1, -t30 * qJD(3), t40, t16 * qJD(3) + t33 * t30, -t17 * qJD(3) + t33 * t32, -t16 * t32 - t17 * t30, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t20, t22 * t29, t20 ^ 2 / 0.2e1, -t20 * t29, t29 ^ 2 / 0.2e1, t12 * t20 + t5 * t29, t12 * t22 - t6 * t29, -t6 * t20 - t5 * t22, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t27, t8 ^ 2 / 0.2e1, -t8 * t27, t27 ^ 2 / 0.2e1, t1 * t27 + t7 * t8, t7 * t10 - t2 * t27, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t9;
