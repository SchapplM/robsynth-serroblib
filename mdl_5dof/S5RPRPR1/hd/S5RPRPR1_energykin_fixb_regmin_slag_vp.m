% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:51
% EndTime: 2021-01-15 11:33:51
% DurationCPUTime: 0.08s
% Computational Cost: add. (154->34), mult. (330->75), div. (0->0), fcn. (188->6), ass. (0->28)
t41 = qJD(1) ^ 2;
t46 = t41 / 0.2e1;
t40 = cos(qJ(3));
t30 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t42 = -qJ(4) * qJD(1) + t30;
t24 = qJD(3) * pkin(3) + t42 * t40;
t38 = sin(qJ(3));
t26 = t42 * t38;
t35 = sin(pkin(8));
t36 = cos(pkin(8));
t17 = t35 * t24 + t36 * t26;
t45 = t41 * qJ(2);
t44 = qJD(3) * t30;
t43 = qJD(1) * qJD(3);
t29 = qJD(4) + (pkin(3) * t38 + qJ(2)) * qJD(1);
t16 = t36 * t24 - t35 * t26;
t39 = cos(qJ(5));
t37 = sin(qJ(5));
t33 = qJD(3) + qJD(5);
t32 = -qJD(1) * pkin(1) + qJD(2);
t28 = (-t35 * t38 + t36 * t40) * qJD(1);
t27 = (t35 * t40 + t36 * t38) * qJD(1);
t20 = t27 * pkin(4) + t29;
t19 = -t37 * t27 + t39 * t28;
t18 = t39 * t27 + t37 * t28;
t15 = -t27 * pkin(7) + t17;
t14 = qJD(3) * pkin(4) - t28 * pkin(7) + t16;
t1 = [t46, 0, 0, t32 * qJD(1), t45, qJ(2) ^ 2 * t46 + t32 ^ 2 / 0.2e1, t40 ^ 2 * t46, -t40 * t41 * t38, t40 * t43, -t38 * t43, qJD(3) ^ 2 / 0.2e1, t38 * t45 + t40 * t44, -t38 * t44 + t40 * t45, t16 * qJD(3) + t29 * t27, -t17 * qJD(3) + t29 * t28, -t16 * t28 - t17 * t27, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t18, t19 * t33, -t18 * t33, t33 ^ 2 / 0.2e1, t20 * t18 + (t39 * t14 - t37 * t15) * t33, t20 * t19 - (t37 * t14 + t39 * t15) * t33;];
T_reg = t1;
