% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:26
% EndTime: 2021-01-15 18:08:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (133->31), mult. (306->74), div. (0->0), fcn. (166->6), ass. (0->28)
t40 = qJD(1) ^ 2;
t51 = t40 / 0.2e1;
t50 = cos(qJ(4));
t35 = sin(pkin(8));
t29 = (pkin(1) * t35 + pkin(6)) * qJD(1);
t38 = sin(qJ(3));
t39 = cos(qJ(3));
t48 = t38 * qJD(2) + t39 * t29;
t23 = qJD(3) * pkin(7) + t48;
t36 = cos(pkin(8));
t44 = -pkin(1) * t36 - pkin(2);
t24 = (-pkin(3) * t39 - pkin(7) * t38 + t44) * qJD(1);
t37 = sin(qJ(4));
t49 = t50 * t23 + t37 * t24;
t47 = qJD(1) * t38;
t46 = t39 * qJD(1);
t45 = qJD(1) * qJD(3);
t43 = -t37 * t23 + t50 * t24;
t42 = t39 * qJD(2) - t38 * t29;
t22 = -qJD(3) * pkin(3) - t42;
t31 = -qJD(4) + t46;
t30 = t44 * qJD(1);
t28 = t37 * qJD(3) + t50 * t47;
t27 = -t50 * qJD(3) + t37 * t47;
t18 = t27 * pkin(4) + qJD(5) + t22;
t17 = -t27 * qJ(5) + t49;
t16 = -t31 * pkin(4) - t28 * qJ(5) + t43;
t1 = [t51, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t35 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t40, t38 ^ 2 * t51, t38 * t40 * t39, t38 * t45, t39 * t45, qJD(3) ^ 2 / 0.2e1, t42 * qJD(3) - t30 * t46, -t48 * qJD(3) + t30 * t47, t28 ^ 2 / 0.2e1, -t28 * t27, -t28 * t31, t27 * t31, t31 ^ 2 / 0.2e1, t22 * t27 - t43 * t31, t22 * t28 + t49 * t31, -t16 * t31 + t18 * t27, t17 * t31 + t18 * t28, -t16 * t28 - t17 * t27, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1;];
T_reg = t1;
