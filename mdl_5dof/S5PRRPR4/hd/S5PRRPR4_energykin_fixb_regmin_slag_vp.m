% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:40
% EndTime: 2021-01-15 15:52:40
% DurationCPUTime: 0.08s
% Computational Cost: add. (135->32), mult. (332->77), div. (0->0), fcn. (219->8), ass. (0->33)
t40 = qJD(2) ^ 2;
t49 = t40 / 0.2e1;
t35 = sin(qJ(3));
t36 = sin(qJ(2));
t29 = qJD(2) * pkin(6) + t36 * qJD(1);
t41 = qJ(4) * qJD(2) + t29;
t23 = qJD(3) * pkin(3) - t41 * t35;
t38 = cos(qJ(3));
t25 = t41 * t38;
t33 = sin(pkin(9));
t48 = cos(pkin(9));
t16 = t33 * t23 + t48 * t25;
t47 = qJD(2) * t35;
t46 = qJD(2) * t38;
t45 = qJD(3) * t29;
t39 = cos(qJ(2));
t44 = t39 * qJD(1);
t43 = qJD(1) * qJD(2);
t42 = qJD(2) * qJD(3);
t15 = t48 * t23 - t33 * t25;
t28 = -t44 + qJD(4) + (-pkin(3) * t38 - pkin(2)) * qJD(2);
t37 = cos(qJ(5));
t34 = sin(qJ(5));
t32 = qJD(3) + qJD(5);
t30 = -qJD(2) * pkin(2) - t44;
t27 = (t33 * t38 + t48 * t35) * qJD(2);
t26 = t33 * t47 - t48 * t46;
t19 = t26 * pkin(4) + t28;
t18 = -t34 * t26 + t37 * t27;
t17 = t37 * t26 + t34 * t27;
t14 = -t26 * pkin(7) + t16;
t13 = qJD(3) * pkin(4) - t27 * pkin(7) + t15;
t1 = [qJD(1) ^ 2 / 0.2e1, t49, t39 * t43, -t36 * t43, t35 ^ 2 * t49, t35 * t40 * t38, t35 * t42, t38 * t42, qJD(3) ^ 2 / 0.2e1, -t30 * t46 - t35 * t45, t30 * t47 - t38 * t45, t15 * qJD(3) + t28 * t26, -t16 * qJD(3) + t28 * t27, -t15 * t27 - t16 * t26, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1, t18 ^ 2 / 0.2e1, -t18 * t17, t18 * t32, -t17 * t32, t32 ^ 2 / 0.2e1, (t37 * t13 - t34 * t14) * t32 + t19 * t17, -(t34 * t13 + t37 * t14) * t32 + t19 * t18;];
T_reg = t1;
