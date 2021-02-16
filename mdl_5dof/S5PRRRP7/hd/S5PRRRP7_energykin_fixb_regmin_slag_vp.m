% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:22
% EndTime: 2021-01-15 16:45:22
% DurationCPUTime: 0.07s
% Computational Cost: add. (132->31), mult. (315->72), div. (0->0), fcn. (208->8), ass. (0->31)
t37 = qJD(2) ^ 2;
t50 = t37 / 0.2e1;
t49 = cos(qJ(4));
t34 = sin(qJ(2));
t46 = qJD(1) * sin(pkin(5));
t24 = qJD(2) * pkin(7) + t34 * t46;
t33 = sin(qJ(3));
t35 = cos(qJ(3));
t45 = qJD(1) * cos(pkin(5));
t47 = t35 * t24 + t33 * t45;
t16 = qJD(3) * pkin(8) + t47;
t36 = cos(qJ(2));
t41 = t36 * t46;
t19 = -t41 + (-pkin(3) * t35 - pkin(8) * t33 - pkin(2)) * qJD(2);
t32 = sin(qJ(4));
t48 = t49 * t16 + t32 * t19;
t44 = qJD(2) * t33;
t43 = t35 * qJD(2);
t42 = qJD(2) * qJD(3);
t40 = qJD(2) * t46;
t39 = -t32 * t16 + t49 * t19;
t38 = -t33 * t24 + t35 * t45;
t15 = -qJD(3) * pkin(3) - t38;
t27 = -qJD(4) + t43;
t25 = -qJD(2) * pkin(2) - t41;
t23 = t32 * qJD(3) + t49 * t44;
t22 = -t49 * qJD(3) + t32 * t44;
t13 = t22 * pkin(4) + qJD(5) + t15;
t12 = -t22 * qJ(5) + t48;
t11 = -t27 * pkin(4) - t23 * qJ(5) + t39;
t1 = [qJD(1) ^ 2 / 0.2e1, t50, t36 * t40, -t34 * t40, t33 ^ 2 * t50, t33 * t37 * t35, t33 * t42, t35 * t42, qJD(3) ^ 2 / 0.2e1, t38 * qJD(3) - t25 * t43, -t47 * qJD(3) + t25 * t44, t23 ^ 2 / 0.2e1, -t23 * t22, -t23 * t27, t22 * t27, t27 ^ 2 / 0.2e1, t15 * t22 - t39 * t27, t15 * t23 + t48 * t27, -t11 * t27 + t13 * t22, t12 * t27 + t13 * t23, -t11 * t23 - t12 * t22, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg = t1;
