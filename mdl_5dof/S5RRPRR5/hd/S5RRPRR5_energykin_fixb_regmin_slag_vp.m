% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR5
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:54
% EndTime: 2022-01-20 11:02:54
% DurationCPUTime: 0.16s
% Computational Cost: add. (206->33), mult. (306->72), div. (0->0), fcn. (195->8), ass. (0->31)
t47 = cos(qJ(4));
t31 = qJD(1) + qJD(2);
t33 = cos(pkin(9));
t46 = t31 * t33;
t32 = sin(pkin(9));
t44 = pkin(1) * qJD(1);
t43 = sin(qJ(2)) * t44;
t26 = qJ(3) * t31 + t43;
t41 = pkin(7) * t31 + t26;
t18 = t41 * t32;
t19 = t41 * t33;
t35 = sin(qJ(4));
t45 = -t35 * t18 + t47 * t19;
t42 = cos(qJ(2)) * t44;
t40 = -t47 * t18 - t19 * t35;
t39 = qJD(3) - t42;
t21 = (-pkin(3) * t33 - pkin(2)) * t31 + t39;
t37 = cos(qJ(5));
t34 = sin(qJ(5));
t30 = qJD(4) + qJD(5);
t29 = t33 ^ 2;
t28 = t32 ^ 2;
t24 = -t31 * pkin(2) + t39;
t23 = (t47 * t32 + t33 * t35) * t31;
t22 = t31 * t32 * t35 - t47 * t46;
t14 = t22 * pkin(4) + t21;
t13 = -t22 * t34 + t23 * t37;
t12 = t37 * t22 + t23 * t34;
t11 = -pkin(8) * t22 + t45;
t10 = qJD(4) * pkin(4) - pkin(8) * t23 + t40;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t31 ^ 2 / 0.2e1, t31 * t42, -t31 * t43, -t24 * t46, (t28 + t29) * t31 * t26, t24 ^ 2 / 0.2e1 + (t29 / 0.2e1 + t28 / 0.2e1) * t26 ^ 2, t23 ^ 2 / 0.2e1, -t23 * t22, t23 * qJD(4), -t22 * qJD(4), qJD(4) ^ 2 / 0.2e1, t40 * qJD(4) + t21 * t22, -t45 * qJD(4) + t21 * t23, t13 ^ 2 / 0.2e1, -t13 * t12, t13 * t30, -t12 * t30, t30 ^ 2 / 0.2e1, t14 * t12 + (t10 * t37 - t11 * t34) * t30, t14 * t13 - (t10 * t34 + t11 * t37) * t30;];
T_reg = t1;
