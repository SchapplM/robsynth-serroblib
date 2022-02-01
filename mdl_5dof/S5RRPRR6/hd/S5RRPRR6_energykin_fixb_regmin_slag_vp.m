% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:24
% EndTime: 2022-01-20 11:17:24
% DurationCPUTime: 0.08s
% Computational Cost: add. (171->32), mult. (256->76), div. (0->0), fcn. (144->8), ass. (0->36)
t33 = qJD(1) + qJD(2);
t49 = pkin(1) * qJD(1);
t45 = sin(qJ(2)) * t49;
t25 = t33 * qJ(3) + t45;
t56 = t25 * t33;
t35 = cos(pkin(9));
t55 = t25 * t35;
t30 = t33 ^ 2;
t34 = sin(pkin(9));
t31 = t34 ^ 2;
t54 = t30 * t31;
t53 = t33 * t34;
t37 = sin(qJ(4));
t52 = t33 * t37;
t51 = t35 * t33;
t44 = cos(qJ(2)) * t49;
t42 = qJD(3) - t44;
t18 = (-pkin(3) * t35 - pkin(7) * t34 - pkin(2)) * t33 + t42;
t40 = cos(qJ(4));
t50 = t37 * t18 + t40 * t55;
t48 = t31 * t56;
t47 = t34 * t52;
t46 = t40 * t53;
t28 = -qJD(4) + t51;
t43 = t40 * t18 - t37 * t55;
t39 = cos(qJ(5));
t36 = sin(qJ(5));
t32 = t35 ^ 2;
t26 = -qJD(5) + t28;
t23 = -t33 * pkin(2) + t42;
t21 = (-t36 * t37 + t39 * t40) * t53;
t20 = (t36 * t40 + t37 * t39) * t53;
t19 = (pkin(4) * t52 + t25) * t34;
t15 = -pkin(8) * t47 + t50;
t14 = -t28 * pkin(4) - pkin(8) * t46 + t43;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t30 / 0.2e1, t33 * t44, -t33 * t45, -t23 * t51, (t31 + t32) * t56, t23 ^ 2 / 0.2e1 + (t32 / 0.2e1 + t31 / 0.2e1) * t25 ^ 2, t40 ^ 2 * t54 / 0.2e1, -t40 * t37 * t54, -t28 * t46, t28 * t47, t28 ^ 2 / 0.2e1, -t43 * t28 + t37 * t48, t50 * t28 + t40 * t48, t21 ^ 2 / 0.2e1, -t21 * t20, -t21 * t26, t20 * t26, t26 ^ 2 / 0.2e1, -(t39 * t14 - t36 * t15) * t26 + t19 * t20, (t36 * t14 + t39 * t15) * t26 + t19 * t21;];
T_reg = t1;
