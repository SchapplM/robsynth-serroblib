% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR6_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:32
% EndTime: 2019-03-09 01:51:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (180->46), mult. (339->94), div. (0->0), fcn. (119->4), ass. (0->41)
t31 = qJD(1) ^ 2;
t24 = t31 / 0.2e1;
t46 = cos(qJ(6));
t45 = -pkin(1) - qJ(3);
t29 = sin(qJ(4));
t44 = qJD(1) * t29;
t21 = qJ(2) * qJD(1) + qJD(3);
t13 = -pkin(7) * qJD(1) + t21;
t43 = qJD(4) * t13;
t14 = -t45 * qJD(1) - qJD(2);
t42 = t14 * qJD(1);
t30 = cos(qJ(4));
t41 = t30 * qJD(1);
t40 = t14 ^ 2 / 0.2e1;
t39 = pkin(4) * t44 - qJD(2);
t38 = qJD(4) * qJ(5);
t37 = qJD(1) * qJD(4);
t36 = t30 * t31 * t29;
t35 = t29 * t37;
t34 = t30 * t37;
t33 = pkin(5) * qJD(1) - t13;
t32 = -qJ(5) * t30 - t45;
t28 = sin(qJ(6));
t27 = t30 ^ 2;
t26 = t29 ^ 2;
t23 = qJD(4) ^ 2 / 0.2e1;
t22 = -qJD(1) * pkin(1) + qJD(2);
t19 = t27 * t24;
t18 = t26 * t24;
t17 = qJD(6) + t41;
t11 = t46 * qJD(4) + t28 * t44;
t9 = t28 * qJD(4) - t46 * t44;
t8 = -t29 * t13 - t38;
t7 = -qJD(4) * pkin(4) - t30 * t13 + qJD(5);
t6 = t32 * qJD(1) + t39;
t5 = -t33 * t29 + t38;
t4 = qJD(5) + t33 * t30 + (-pkin(4) - pkin(8)) * qJD(4);
t3 = (pkin(8) * t29 + t32) * qJD(1) + t39;
t2 = t28 * t4 + t46 * t3;
t1 = -t28 * t3 + t46 * t4;
t10 = [0, 0, 0, 0, 0, t24, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, t22 * qJD(1), t31 * qJ(2), qJ(2) ^ 2 * t24 + t22 ^ 2 / 0.2e1, t24, 0, 0, 0, 0, 0, 0, t21 * qJD(1), t42, t40 + t21 ^ 2 / 0.2e1, t19, -t36, t34, t18, -t35, t23, t29 * t42 + t30 * t43, t14 * t41 - t29 * t43 (-t26 - t27) * t13 * qJD(1), t40 + (t26 / 0.2e1 + t27 / 0.2e1) * t13 ^ 2, t23, -t34, t35, t19, -t36, t18 (t29 * t8 + t30 * t7) * qJD(1), t7 * qJD(4) - t6 * t44, -t8 * qJD(4) - t6 * t41, t6 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t11 ^ 2 / 0.2e1, -t11 * t9, t11 * t17, t9 ^ 2 / 0.2e1, -t9 * t17, t17 ^ 2 / 0.2e1, t1 * t17 + t5 * t9, t5 * t11 - t2 * t17, -t1 * t11 - t2 * t9, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg  = t10;
