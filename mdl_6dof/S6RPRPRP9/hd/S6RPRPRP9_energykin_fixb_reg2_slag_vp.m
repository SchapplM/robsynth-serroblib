% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP9_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:26
% EndTime: 2019-03-09 03:29:26
% DurationCPUTime: 0.14s
% Computational Cost: add. (456->53), mult. (941->113), div. (0->0), fcn. (558->6), ass. (0->44)
t43 = qJD(1) ^ 2;
t36 = t43 / 0.2e1;
t40 = sin(qJ(5));
t53 = cos(qJ(5));
t41 = sin(qJ(3));
t42 = cos(qJ(3));
t23 = (pkin(3) * t41 - qJ(4) * t42 + qJ(2)) * qJD(1);
t30 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t24 = qJD(3) * qJ(4) + t41 * t30;
t39 = sin(pkin(9));
t49 = cos(pkin(9));
t11 = t49 * t23 - t39 * t24;
t48 = qJD(1) * t42;
t27 = t39 * qJD(3) + t49 * t48;
t46 = t41 * qJD(1);
t7 = pkin(4) * t46 - t27 * pkin(8) + t11;
t12 = t39 * t23 + t49 * t24;
t25 = -t49 * qJD(3) + t39 * t48;
t9 = -t25 * pkin(8) + t12;
t4 = t40 * t7 + t53 * t9;
t14 = t53 * t25 + t40 * t27;
t16 = -t40 * t25 + t53 * t27;
t52 = t16 * t14;
t31 = qJD(5) + t46;
t51 = t31 * t14;
t50 = t43 * qJ(2);
t47 = qJD(3) * t30;
t45 = t14 ^ 2 / 0.2e1;
t44 = qJD(1) * qJD(3);
t3 = -t40 * t9 + t53 * t7;
t22 = -qJD(3) * pkin(3) - t42 * t30 + qJD(4);
t17 = t25 * pkin(4) + t22;
t38 = t42 ^ 2;
t37 = t41 ^ 2;
t35 = qJ(2) ^ 2 * t36;
t33 = -pkin(1) * qJD(1) + qJD(2);
t32 = t37 * t36;
t28 = t31 ^ 2 / 0.2e1;
t13 = t16 ^ 2 / 0.2e1;
t10 = t16 * t31;
t5 = t14 * pkin(5) - t16 * qJ(6) + t17;
t2 = t31 * qJ(6) + t4;
t1 = -t31 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t36, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, t33 * qJD(1), t50, t35 + t33 ^ 2 / 0.2e1, t38 * t36, -t42 * t43 * t41, t42 * t44, t32, -t41 * t44, qJD(3) ^ 2 / 0.2e1, t41 * t50 + t42 * t47, -t41 * t47 + t42 * t50 (-t37 - t38) * t30 * qJD(1), t35 + (t37 / 0.2e1 + t38 / 0.2e1) * t30 ^ 2, t27 ^ 2 / 0.2e1, -t27 * t25, t27 * t46, t25 ^ 2 / 0.2e1, -t25 * t46, t32, t11 * t46 + t22 * t25, -t12 * t46 + t22 * t27, -t11 * t27 - t12 * t25, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t13, -t52, t10, t45, -t51, t28, t17 * t14 + t3 * t31, t17 * t16 - t4 * t31, -t4 * t14 - t3 * t16, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t13, t10, t52, t28, t51, t45, -t1 * t31 + t5 * t14, t1 * t16 - t2 * t14, -t5 * t16 + t2 * t31, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
