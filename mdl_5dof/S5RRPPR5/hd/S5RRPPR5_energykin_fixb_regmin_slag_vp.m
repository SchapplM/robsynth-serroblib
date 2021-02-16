% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR5
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
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:29
% EndTime: 2021-01-15 19:36:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (179->39), mult. (460->83), div. (0->0), fcn. (283->6), ass. (0->35)
t39 = sin(pkin(8));
t40 = cos(pkin(8));
t42 = sin(qJ(2));
t44 = cos(qJ(2));
t29 = (t39 * t44 + t40 * t42) * qJD(1);
t34 = -(pkin(2) * t44 + pkin(1)) * qJD(1) + qJD(3);
t58 = -t29 * qJ(4) + t34;
t45 = qJD(1) ^ 2;
t57 = t45 / 0.2e1;
t56 = -pkin(3) - pkin(4);
t55 = t44 * t45;
t54 = pkin(6) + qJ(3);
t52 = qJD(1) * t42;
t32 = qJD(2) * pkin(2) - t54 * t52;
t51 = qJD(1) * t44;
t33 = t54 * t51;
t25 = t39 * t32 + t40 * t33;
t50 = qJD(1) * qJD(2);
t23 = qJD(2) * qJ(4) + t25;
t48 = t42 * t50;
t47 = t44 * t50;
t24 = t40 * t32 - t39 * t33;
t46 = qJD(4) - t24;
t43 = cos(qJ(5));
t41 = sin(qJ(5));
t36 = qJD(2) - qJD(5);
t28 = t39 * t52 - t40 * t51;
t22 = t41 * t28 + t43 * t29;
t21 = -t43 * t28 + t41 * t29;
t20 = -qJD(2) * pkin(3) + t46;
t19 = t28 * pkin(3) + t58;
t18 = t28 * pkin(7) + t23;
t17 = -t29 * pkin(7) + t56 * qJD(2) + t46;
t16 = t56 * t28 - t58;
t1 = [t57, 0, 0, t42 ^ 2 * t57, t42 * t55, t48, t47, qJD(2) ^ 2 / 0.2e1, pkin(1) * t55 - pkin(6) * t48, -t45 * pkin(1) * t42 - pkin(6) * t47, t24 * qJD(2) + t34 * t28, -t25 * qJD(2) + t34 * t29, -t24 * t29 - t25 * t28, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, -t20 * qJD(2) + t19 * t28, t20 * t29 - t23 * t28, t23 * qJD(2) - t19 * t29, t23 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t22 ^ 2 / 0.2e1, -t22 * t21, -t22 * t36, t21 * t36, t36 ^ 2 / 0.2e1, t16 * t21 - (t43 * t17 - t41 * t18) * t36, t16 * t22 + (t41 * t17 + t43 * t18) * t36;];
T_reg = t1;
