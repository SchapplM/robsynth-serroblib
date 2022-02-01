% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:19
% EndTime: 2022-01-20 10:34:19
% DurationCPUTime: 0.05s
% Computational Cost: add. (101->18), mult. (163->49), div. (0->0), fcn. (81->8), ass. (0->23)
t17 = qJD(1) + qJD(2);
t16 = qJD(4) + t17;
t15 = t16 ^ 2;
t33 = t15 / 0.2e1;
t30 = pkin(1) * qJD(1);
t27 = cos(qJ(2)) * t30;
t14 = t17 * pkin(2) + t27;
t18 = sin(pkin(9));
t19 = cos(pkin(9));
t28 = sin(qJ(2)) * t30;
t12 = t18 * t14 + t19 * t28;
t21 = sin(qJ(4));
t24 = cos(qJ(4));
t11 = t19 * t14 - t18 * t28;
t9 = t17 * pkin(3) + t11;
t26 = -t21 * t12 + t24 * t9;
t32 = (-t16 * pkin(4) - t26) * t16;
t31 = t24 * t12 + t21 * t9;
t29 = qJD(5) * t16;
t23 = cos(qJ(5));
t20 = sin(qJ(5));
t7 = t16 * pkin(8) + t31;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t17 ^ 2 / 0.2e1, t17 * t27, -t17 * t28, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t33, t26 * t16, -t31 * t16, t20 ^ 2 * t33, t20 * t15 * t23, t20 * t29, t23 * t29, qJD(5) ^ 2 / 0.2e1, -t23 * t32 + (t23 * qJD(3) - t20 * t7) * qJD(5), t20 * t32 - (t20 * qJD(3) + t23 * t7) * qJD(5);];
T_reg = t1;
