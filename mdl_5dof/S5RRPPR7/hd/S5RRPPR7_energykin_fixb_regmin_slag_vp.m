% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:00
% EndTime: 2021-01-15 19:59:00
% DurationCPUTime: 0.08s
% Computational Cost: add. (177->39), mult. (452->83), div. (0->0), fcn. (275->6), ass. (0->35)
t45 = qJD(1) ^ 2;
t56 = t45 / 0.2e1;
t55 = pkin(3) + pkin(7);
t44 = cos(qJ(2));
t54 = t44 * t45;
t53 = pkin(6) + qJ(3);
t42 = sin(qJ(2));
t52 = qJD(1) * t42;
t35 = qJD(2) * pkin(2) - t53 * t52;
t51 = qJD(1) * t44;
t36 = t53 * t51;
t39 = sin(pkin(8));
t40 = cos(pkin(8));
t25 = t39 * t35 + t40 * t36;
t50 = qJD(1) * qJD(2);
t49 = t42 * t50;
t48 = t44 * t50;
t24 = t40 * t35 - t39 * t36;
t47 = qJD(4) - t24;
t23 = -qJD(2) * qJ(4) - t25;
t32 = (t39 * t44 + t40 * t42) * qJD(1);
t37 = qJD(3) + (-pkin(2) * t44 - pkin(1)) * qJD(1);
t46 = -t32 * qJ(4) + t37;
t43 = cos(qJ(5));
t41 = sin(qJ(5));
t31 = t39 * t52 - t40 * t51;
t30 = qJD(5) + t32;
t27 = t43 * qJD(2) + t41 * t31;
t26 = t41 * qJD(2) - t43 * t31;
t22 = -qJD(2) * pkin(3) + t47;
t21 = t31 * pkin(3) + t46;
t20 = -t31 * pkin(4) - t23;
t19 = t32 * pkin(4) - t55 * qJD(2) + t47;
t18 = t55 * t31 + t46;
t1 = [t56, 0, 0, t42 ^ 2 * t56, t42 * t54, t49, t48, qJD(2) ^ 2 / 0.2e1, pkin(1) * t54 - pkin(6) * t49, -t45 * pkin(1) * t42 - pkin(6) * t48, t24 * qJD(2) + t37 * t31, -t25 * qJD(2) + t37 * t32, -t24 * t32 - t25 * t31, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t22 * t32 + t23 * t31, t22 * qJD(2) - t21 * t31, -t23 * qJD(2) - t21 * t32, t21 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1, -t27 * t26, t27 * t30, -t26 * t30, t30 ^ 2 / 0.2e1, (-t41 * t18 + t43 * t19) * t30 + t20 * t26, -(t43 * t18 + t41 * t19) * t30 + t20 * t27;];
T_reg = t1;
