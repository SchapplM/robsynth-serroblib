% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:11
% EndTime: 2019-03-09 04:52:11
% DurationCPUTime: 0.15s
% Computational Cost: add. (298->50), mult. (588->98), div. (0->0), fcn. (296->4), ass. (0->39)
t37 = qJD(1) ^ 2;
t31 = t37 / 0.2e1;
t46 = -pkin(4) - pkin(5);
t45 = cos(qJ(4));
t34 = sin(qJ(4));
t36 = cos(qJ(3));
t42 = qJD(1) * t36;
t19 = -t45 * qJD(3) + t34 * t42;
t21 = t34 * qJD(3) + t45 * t42;
t9 = t21 * t19;
t35 = sin(qJ(3));
t26 = t35 * qJD(1) + qJD(4);
t10 = t21 * t26;
t44 = t26 * t19;
t14 = (pkin(3) * t35 - pkin(8) * t36 + qJ(2)) * qJD(1);
t25 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t15 = qJD(3) * pkin(8) + t35 * t25;
t8 = t34 * t14 + t45 * t15;
t43 = t37 * qJ(2);
t41 = qJD(3) * t25;
t17 = t19 ^ 2 / 0.2e1;
t22 = t26 ^ 2 / 0.2e1;
t40 = qJD(1) * qJD(3);
t5 = t26 * qJ(5) + t8;
t7 = t45 * t14 - t34 * t15;
t16 = -qJD(3) * pkin(3) - t36 * t25;
t39 = qJD(5) - t7;
t38 = t21 * qJ(5) - t16;
t33 = t36 ^ 2;
t32 = t35 ^ 2;
t29 = qJ(2) ^ 2 * t31;
t28 = -pkin(1) * qJD(1) + qJD(2);
t18 = t21 ^ 2 / 0.2e1;
t6 = t19 * pkin(4) - t38;
t4 = -t26 * pkin(4) + t39;
t3 = t46 * t19 + qJD(6) + t38;
t2 = t19 * qJ(6) + t5;
t1 = -t21 * qJ(6) + t46 * t26 + t39;
t11 = [0, 0, 0, 0, 0, t31, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, t28 * qJD(1), t43, t29 + t28 ^ 2 / 0.2e1, t33 * t31, -t36 * t37 * t35, t36 * t40, t32 * t31, -t35 * t40, qJD(3) ^ 2 / 0.2e1, t35 * t43 + t36 * t41, -t35 * t41 + t36 * t43 (-t32 - t33) * t25 * qJD(1), t29 + (t32 / 0.2e1 + t33 / 0.2e1) * t25 ^ 2, t18, -t9, t10, t17, -t44, t22, t16 * t19 + t7 * t26, t16 * t21 - t8 * t26, -t8 * t19 - t7 * t21, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1, t18, t10, t9, t22, t44, t17, t6 * t19 - t4 * t26, -t5 * t19 + t4 * t21, -t6 * t21 + t5 * t26, t5 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t18, t9, -t10, t17, -t44, t22, -t1 * t26 - t3 * t19, t2 * t26 + t3 * t21, -t1 * t21 + t2 * t19, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t11;
