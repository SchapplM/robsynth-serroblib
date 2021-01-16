% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:10:23
% EndTime: 2021-01-16 00:10:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (254->37), mult. (565->82), div. (0->0), fcn. (377->6), ass. (0->35)
t46 = qJD(1) ^ 2;
t59 = t46 / 0.2e1;
t58 = pkin(7) + pkin(6);
t57 = cos(qJ(4));
t45 = cos(qJ(2));
t56 = t45 * t46;
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t52 = qJD(1) * t45;
t43 = sin(qJ(2));
t53 = qJD(1) * t43;
t31 = t42 * t53 - t44 * t52;
t32 = (t42 * t45 + t43 * t44) * qJD(1);
t37 = (-pkin(2) * t45 - pkin(1)) * qJD(1);
t23 = t31 * pkin(3) - t32 * pkin(8) + t37;
t40 = qJD(2) + qJD(3);
t35 = qJD(2) * pkin(2) - t58 * t53;
t36 = t58 * t52;
t54 = t42 * t35 + t44 * t36;
t26 = t40 * pkin(8) + t54;
t41 = sin(qJ(4));
t55 = t41 * t23 + t57 * t26;
t51 = qJD(1) * qJD(2);
t50 = t43 * t51;
t49 = t45 * t51;
t48 = t57 * t23 - t41 * t26;
t47 = t44 * t35 - t42 * t36;
t25 = -t40 * pkin(3) - t47;
t30 = qJD(4) + t31;
t28 = t57 * t32 + t41 * t40;
t27 = t41 * t32 - t57 * t40;
t20 = t27 * pkin(4) + qJD(5) + t25;
t19 = -t27 * qJ(5) + t55;
t18 = t30 * pkin(4) - t28 * qJ(5) + t48;
t1 = [t59, 0, 0, t43 ^ 2 * t59, t43 * t56, t50, t49, qJD(2) ^ 2 / 0.2e1, pkin(1) * t56 - pkin(6) * t50, -t46 * pkin(1) * t43 - pkin(6) * t49, t32 ^ 2 / 0.2e1, -t32 * t31, t32 * t40, -t31 * t40, t40 ^ 2 / 0.2e1, t37 * t31 + t47 * t40, t37 * t32 - t54 * t40, t28 ^ 2 / 0.2e1, -t28 * t27, t28 * t30, -t27 * t30, t30 ^ 2 / 0.2e1, t25 * t27 + t48 * t30, t25 * t28 - t55 * t30, t18 * t30 + t20 * t27, -t19 * t30 + t20 * t28, -t18 * t28 - t19 * t27, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1;];
T_reg = t1;
