% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:25:52
% EndTime: 2021-01-15 19:25:52
% DurationCPUTime: 0.06s
% Computational Cost: add. (126->31), mult. (253->68), div. (0->0), fcn. (124->4), ass. (0->24)
t27 = qJD(1) ^ 2;
t35 = t27 / 0.2e1;
t34 = cos(qJ(4));
t25 = sin(qJ(3));
t26 = cos(qJ(3));
t15 = (pkin(3) * t25 - pkin(7) * t26 + qJ(2)) * qJD(1);
t20 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t16 = qJD(3) * pkin(7) + t25 * t20;
t24 = sin(qJ(4));
t33 = t24 * t15 + t34 * t16;
t32 = t27 * qJ(2);
t31 = qJD(1) * t26;
t30 = qJD(3) * t20;
t29 = qJD(1) * qJD(3);
t28 = t34 * t15 - t24 * t16;
t17 = -qJD(3) * pkin(3) - t26 * t20;
t22 = -qJD(1) * pkin(1) + qJD(2);
t21 = t25 * qJD(1) + qJD(4);
t19 = t24 * qJD(3) + t34 * t31;
t18 = -t34 * qJD(3) + t24 * t31;
t11 = t18 * pkin(4) + qJD(5) + t17;
t10 = -t18 * qJ(5) + t33;
t9 = t21 * pkin(4) - t19 * qJ(5) + t28;
t1 = [t35, 0, 0, t22 * qJD(1), t32, qJ(2) ^ 2 * t35 + t22 ^ 2 / 0.2e1, t26 ^ 2 * t35, -t26 * t27 * t25, t26 * t29, -t25 * t29, qJD(3) ^ 2 / 0.2e1, t25 * t32 + t26 * t30, -t25 * t30 + t26 * t32, t19 ^ 2 / 0.2e1, -t19 * t18, t19 * t21, -t18 * t21, t21 ^ 2 / 0.2e1, t17 * t18 + t28 * t21, t17 * t19 - t33 * t21, t11 * t18 + t9 * t21, -t10 * t21 + t11 * t19, -t10 * t18 - t9 * t19, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1;];
T_reg = t1;
