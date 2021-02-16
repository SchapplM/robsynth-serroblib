% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR14
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
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:02
% EndTime: 2021-01-15 12:17:02
% DurationCPUTime: 0.08s
% Computational Cost: add. (150->34), mult. (314->75), div. (0->0), fcn. (172->6), ass. (0->28)
t49 = qJD(1) ^ 2;
t54 = t49 / 0.2e1;
t48 = cos(qJ(3));
t38 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t50 = -qJ(4) * qJD(1) + t38;
t32 = qJD(3) * pkin(3) + t50 * t48;
t46 = sin(qJ(3));
t33 = t50 * t46;
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t26 = t43 * t32 + t44 * t33;
t53 = t49 * qJ(2);
t52 = qJD(3) * t38;
t51 = qJD(1) * qJD(3);
t37 = qJD(4) + (pkin(3) * t46 + qJ(2)) * qJD(1);
t25 = t44 * t32 - t43 * t33;
t35 = (t43 * t48 + t44 * t46) * qJD(1);
t47 = cos(qJ(5));
t45 = sin(qJ(5));
t40 = -qJD(1) * pkin(1) + qJD(2);
t36 = (-t43 * t46 + t44 * t48) * qJD(1);
t34 = qJD(5) + t35;
t29 = t45 * qJD(3) + t47 * t36;
t28 = -t47 * qJD(3) + t45 * t36;
t27 = t35 * pkin(4) - t36 * pkin(7) + t37;
t24 = qJD(3) * pkin(7) + t26;
t23 = -qJD(3) * pkin(4) - t25;
t1 = [t54, 0, 0, t40 * qJD(1), t53, qJ(2) ^ 2 * t54 + t40 ^ 2 / 0.2e1, t48 ^ 2 * t54, -t48 * t49 * t46, t48 * t51, -t46 * t51, qJD(3) ^ 2 / 0.2e1, t46 * t53 + t48 * t52, -t46 * t52 + t48 * t53, t25 * qJD(3) + t37 * t35, -t26 * qJD(3) + t37 * t36, -t25 * t36 - t26 * t35, t26 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1 + t37 ^ 2 / 0.2e1, t29 ^ 2 / 0.2e1, -t29 * t28, t29 * t34, -t28 * t34, t34 ^ 2 / 0.2e1, (-t45 * t24 + t47 * t27) * t34 + t23 * t28, -(t47 * t24 + t45 * t27) * t34 + t23 * t29;];
T_reg = t1;
