% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:52
% EndTime: 2022-01-20 11:30:52
% DurationCPUTime: 0.05s
% Computational Cost: add. (107->18), mult. (163->49), div. (0->0), fcn. (81->8), ass. (0->23)
t19 = qJD(1) + qJD(2);
t18 = qJD(3) + t19;
t17 = t18 ^ 2;
t34 = t17 / 0.2e1;
t32 = pkin(1) * qJD(1);
t29 = cos(qJ(2)) * t32;
t16 = t19 * pkin(2) + t29;
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t30 = sin(qJ(2)) * t32;
t28 = t26 * t16 - t23 * t30;
t12 = t18 * pkin(3) + t28;
t14 = t23 * t16 + t26 * t30;
t20 = sin(pkin(9));
t21 = cos(pkin(9));
t9 = t21 * t12 - t20 * t14;
t33 = (-t18 * pkin(4) - t9) * t18;
t10 = t20 * t12 + t21 * t14;
t31 = qJD(5) * t18;
t25 = cos(qJ(5));
t22 = sin(qJ(5));
t8 = t18 * pkin(8) + t10;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t19 ^ 2 / 0.2e1, t19 * t29, -t19 * t30, t34, t28 * t18, -t14 * t18, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t22 ^ 2 * t34, t22 * t17 * t25, t22 * t31, t25 * t31, qJD(5) ^ 2 / 0.2e1, -t25 * t33 + (t25 * qJD(4) - t22 * t8) * qJD(5), t22 * t33 - (t22 * qJD(4) + t25 * t8) * qJD(5);];
T_reg = t1;
