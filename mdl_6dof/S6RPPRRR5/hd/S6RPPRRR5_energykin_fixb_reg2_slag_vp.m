% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:59
% EndTime: 2019-03-09 02:29:00
% DurationCPUTime: 0.16s
% Computational Cost: add. (296->45), mult. (562->109), div. (0->0), fcn. (293->6), ass. (0->36)
t34 = sin(qJ(5));
t35 = sin(qJ(4));
t36 = cos(qJ(5));
t37 = cos(qJ(4));
t16 = (t34 * t37 + t35 * t36) * qJD(1);
t22 = (pkin(1) + qJ(3)) * qJD(1) - qJD(2);
t38 = qJD(1) ^ 2;
t29 = t38 / 0.2e1;
t46 = cos(qJ(6));
t26 = qJD(1) * qJ(2) + qJD(3);
t21 = -pkin(7) * qJD(1) + t26;
t40 = -pkin(8) * qJD(1) + t21;
t13 = qJD(4) * pkin(4) + t40 * t37;
t14 = t40 * t35;
t7 = t34 * t13 + t36 * t14;
t44 = qJD(4) * t21;
t43 = t22 * qJD(1);
t42 = t22 ^ 2 / 0.2e1;
t41 = qJD(1) * qJD(4);
t6 = t36 * t13 - t34 * t14;
t19 = t35 * qJD(1) * pkin(4) + t22;
t33 = sin(qJ(6));
t32 = t37 ^ 2;
t31 = t35 ^ 2;
t28 = qJD(4) + qJD(5);
t27 = -qJD(1) * pkin(1) + qJD(2);
t18 = (-t34 * t35 + t36 * t37) * qJD(1);
t15 = qJD(6) + t16;
t10 = t46 * t18 + t33 * t28;
t8 = t33 * t18 - t46 * t28;
t5 = t16 * pkin(5) - t18 * pkin(9) + t19;
t4 = t28 * pkin(9) + t7;
t3 = -t28 * pkin(5) - t6;
t2 = t33 * t5 + t46 * t4;
t1 = -t33 * t4 + t46 * t5;
t9 = [0, 0, 0, 0, 0, t29, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, t27 * qJD(1), t38 * qJ(2), qJ(2) ^ 2 * t29 + t27 ^ 2 / 0.2e1, t29, 0, 0, 0, 0, 0, 0, t26 * qJD(1), t43, t42 + t26 ^ 2 / 0.2e1, t32 * t29, -t37 * t38 * t35, t37 * t41, t31 * t29, -t35 * t41, qJD(4) ^ 2 / 0.2e1, t35 * t43 + t37 * t44, -t35 * t44 + t37 * t43 (-t31 - t32) * t21 * qJD(1), t42 + (t31 / 0.2e1 + t32 / 0.2e1) * t21 ^ 2, t18 ^ 2 / 0.2e1, -t18 * t16, t18 * t28, t16 ^ 2 / 0.2e1, -t16 * t28, t28 ^ 2 / 0.2e1, t19 * t16 + t6 * t28, t19 * t18 - t7 * t28, -t7 * t16 - t6 * t18, t7 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, t10 * t15, t8 ^ 2 / 0.2e1, -t8 * t15, t15 ^ 2 / 0.2e1, t1 * t15 + t3 * t8, t3 * t10 - t2 * t15, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t9;
