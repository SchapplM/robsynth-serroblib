% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:35
% EndTime: 2022-01-23 09:14:35
% DurationCPUTime: 0.07s
% Computational Cost: add. (125->37), mult. (326->79), div. (0->0), fcn. (213->8), ass. (0->30)
t49 = cos(qJ(4));
t37 = sin(pkin(8));
t31 = (pkin(1) * t37 + qJ(3)) * qJD(1);
t38 = cos(pkin(9));
t34 = t38 * qJD(2);
t36 = sin(pkin(9));
t21 = t34 + (-pkin(6) * qJD(1) - t31) * t36;
t24 = t36 * qJD(2) + t38 * t31;
t47 = qJD(1) * t38;
t22 = pkin(6) * t47 + t24;
t41 = sin(qJ(4));
t48 = t41 * t21 + t49 * t22;
t39 = cos(pkin(8));
t46 = -pkin(1) * t39 - pkin(2);
t45 = t49 * t21 - t41 * t22;
t26 = qJD(3) + (-pkin(3) * t38 + t46) * qJD(1);
t43 = qJD(1) ^ 2;
t42 = cos(qJ(5));
t40 = sin(qJ(5));
t35 = qJD(4) + qJD(5);
t30 = t46 * qJD(1) + qJD(3);
t28 = (t49 * t36 + t38 * t41) * qJD(1);
t27 = t41 * t36 * qJD(1) - t49 * t47;
t23 = -t36 * t31 + t34;
t17 = t27 * pkin(4) + t26;
t16 = -t40 * t27 + t42 * t28;
t15 = t42 * t27 + t40 * t28;
t14 = -t27 * pkin(7) + t48;
t13 = qJD(4) * pkin(4) - t28 * pkin(7) + t45;
t1 = [t43 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t37 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t43, -t30 * t47, (-t23 * t36 + t24 * t38) * qJD(1), t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t30 ^ 2 / 0.2e1, t28 ^ 2 / 0.2e1, -t28 * t27, t28 * qJD(4), -t27 * qJD(4), qJD(4) ^ 2 / 0.2e1, t45 * qJD(4) + t26 * t27, -t48 * qJD(4) + t26 * t28, t16 ^ 2 / 0.2e1, -t16 * t15, t16 * t35, -t15 * t35, t35 ^ 2 / 0.2e1, t17 * t15 + (t42 * t13 - t40 * t14) * t35, t17 * t16 - (t40 * t13 + t42 * t14) * t35;];
T_reg = t1;
