% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:41
% EndTime: 2021-01-15 16:33:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (117->28), mult. (276->69), div. (0->0), fcn. (169->6), ass. (0->29)
t30 = qJD(2) ^ 2;
t41 = t30 / 0.2e1;
t40 = cos(qJ(4));
t26 = sin(qJ(3));
t27 = sin(qJ(2));
t21 = qJD(2) * pkin(6) + t27 * qJD(1);
t31 = pkin(7) * qJD(2) + t21;
t16 = qJD(3) * pkin(3) - t31 * t26;
t28 = cos(qJ(3));
t17 = t31 * t28;
t25 = sin(qJ(4));
t39 = t25 * t16 + t40 * t17;
t38 = qJD(2) * t26;
t37 = qJD(2) * t28;
t36 = qJD(3) * t21;
t29 = cos(qJ(2));
t35 = t29 * qJD(1);
t34 = qJD(1) * qJD(2);
t33 = qJD(2) * qJD(3);
t32 = t40 * t16 - t25 * t17;
t20 = -t35 + (-pkin(3) * t28 - pkin(2)) * qJD(2);
t24 = qJD(3) + qJD(4);
t22 = -qJD(2) * pkin(2) - t35;
t19 = (t25 * t28 + t40 * t26) * qJD(2);
t18 = t25 * t38 - t40 * t37;
t12 = t18 * pkin(4) + qJD(5) + t20;
t11 = -t18 * qJ(5) + t39;
t10 = t24 * pkin(4) - t19 * qJ(5) + t32;
t1 = [qJD(1) ^ 2 / 0.2e1, t41, t29 * t34, -t27 * t34, t26 ^ 2 * t41, t26 * t30 * t28, t26 * t33, t28 * t33, qJD(3) ^ 2 / 0.2e1, -t22 * t37 - t26 * t36, t22 * t38 - t28 * t36, t19 ^ 2 / 0.2e1, -t19 * t18, t19 * t24, -t18 * t24, t24 ^ 2 / 0.2e1, t20 * t18 + t32 * t24, t20 * t19 - t39 * t24, t10 * t24 + t12 * t18, -t11 * t24 + t12 * t19, -t10 * t19 - t11 * t18, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1;];
T_reg = t1;
