% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:43:14
% EndTime: 2022-01-20 11:43:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (246->33), mult. (363->77), div. (0->0), fcn. (219->8), ass. (0->33)
t35 = qJD(1) + qJD(2);
t33 = t35 ^ 2;
t52 = t33 / 0.2e1;
t38 = sin(qJ(3));
t51 = t35 * t38;
t41 = cos(qJ(3));
t50 = t35 * t41;
t49 = pkin(1) * qJD(1);
t45 = sin(qJ(2)) * t49;
t30 = t35 * pkin(7) + t45;
t43 = qJ(4) * t35 + t30;
t24 = qJD(3) * pkin(3) - t43 * t38;
t26 = t43 * t41;
t36 = sin(pkin(9));
t48 = cos(pkin(9));
t17 = t36 * t24 + t48 * t26;
t47 = qJD(3) * t38;
t46 = qJD(3) * t41;
t44 = cos(qJ(2)) * t49;
t16 = t48 * t24 - t36 * t26;
t27 = -t44 + qJD(4) + (-pkin(3) * t41 - pkin(2)) * t35;
t40 = cos(qJ(5));
t37 = sin(qJ(5));
t34 = qJD(3) + qJD(5);
t31 = -t35 * pkin(2) - t44;
t29 = (t36 * t41 + t48 * t38) * t35;
t28 = t36 * t51 - t48 * t50;
t20 = t28 * pkin(4) + t27;
t19 = -t37 * t28 + t40 * t29;
t18 = t40 * t28 + t37 * t29;
t15 = -t28 * pkin(8) + t17;
t14 = qJD(3) * pkin(4) - t29 * pkin(8) + t16;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t52, t35 * t44, -t35 * t45, t38 ^ 2 * t52, t38 * t33 * t41, t35 * t47, t35 * t46, qJD(3) ^ 2 / 0.2e1, -t30 * t47 - t31 * t50, -t30 * t46 + t31 * t51, t16 * qJD(3) + t27 * t28, -t17 * qJD(3) + t27 * t29, -t16 * t29 - t17 * t28, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1, t19 ^ 2 / 0.2e1, -t19 * t18, t19 * t34, -t18 * t34, t34 ^ 2 / 0.2e1, t20 * t18 + (t40 * t14 - t37 * t15) * t34, t20 * t19 - (t37 * t14 + t40 * t15) * t34;];
T_reg = t1;
