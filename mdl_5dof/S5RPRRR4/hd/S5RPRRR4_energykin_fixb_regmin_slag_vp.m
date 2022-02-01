% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:45
% EndTime: 2022-01-23 09:34:45
% DurationCPUTime: 0.06s
% Computational Cost: add. (87->18), mult. (167->51), div. (0->0), fcn. (81->8), ass. (0->24)
t17 = qJD(1) + qJD(3);
t16 = qJD(4) + t17;
t15 = t16 ^ 2;
t34 = t15 / 0.2e1;
t19 = cos(pkin(9));
t14 = (pkin(1) * t19 + pkin(2)) * qJD(1);
t22 = sin(qJ(3));
t25 = cos(qJ(3));
t18 = sin(pkin(9));
t30 = pkin(1) * qJD(1) * t18;
t28 = t25 * t14 - t22 * t30;
t10 = t17 * pkin(3) + t28;
t12 = t22 * t14 + t25 * t30;
t21 = sin(qJ(4));
t24 = cos(qJ(4));
t29 = t24 * t10 - t21 * t12;
t33 = (-t16 * pkin(4) - t29) * t16;
t32 = t21 * t10 + t24 * t12;
t31 = qJD(5) * t16;
t26 = qJD(1) ^ 2;
t23 = cos(qJ(5));
t20 = sin(qJ(5));
t8 = t16 * pkin(8) + t32;
t1 = [t26 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t18 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t26, t17 ^ 2 / 0.2e1, t28 * t17, -t12 * t17, t34, t29 * t16, -t32 * t16, t20 ^ 2 * t34, t20 * t15 * t23, t20 * t31, t23 * t31, qJD(5) ^ 2 / 0.2e1, -t23 * t33 + (t23 * qJD(2) - t20 * t8) * qJD(5), t20 * t33 - (t20 * qJD(2) + t23 * t8) * qJD(5);];
T_reg = t1;
