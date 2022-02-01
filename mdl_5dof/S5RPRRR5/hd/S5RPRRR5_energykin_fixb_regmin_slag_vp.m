% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:49:01
% EndTime: 2022-01-20 09:49:01
% DurationCPUTime: 0.16s
% Computational Cost: add. (108->28), mult. (202->69), div. (0->0), fcn. (107->8), ass. (0->30)
t31 = qJD(1) + qJD(3);
t29 = t31 ^ 2;
t49 = t29 / 0.2e1;
t35 = sin(qJ(4));
t48 = t31 * t35;
t38 = cos(qJ(4));
t47 = t31 * t38;
t33 = cos(pkin(9));
t23 = (pkin(1) * t33 + pkin(2)) * qJD(1);
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t32 = sin(pkin(9));
t43 = pkin(1) * qJD(1) * t32;
t45 = t36 * t23 + t39 * t43;
t19 = t31 * pkin(7) + t45;
t46 = t35 * qJD(2) + t38 * t19;
t44 = qJD(4) * t31;
t42 = t39 * t23 - t36 * t43;
t40 = qJD(1) ^ 2;
t37 = cos(qJ(5));
t34 = sin(qJ(5));
t30 = qJD(4) + qJD(5);
t28 = t38 * qJD(2);
t21 = (t34 * t38 + t35 * t37) * t31;
t20 = t34 * t48 - t37 * t47;
t18 = -t31 * pkin(3) - t42;
t16 = (-pkin(4) * t38 - pkin(3)) * t31 - t42;
t15 = pkin(8) * t47 + t46;
t14 = qJD(4) * pkin(4) + t28 + (-pkin(8) * t31 - t19) * t35;
t1 = [t40 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t32 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t40, t49, t42 * t31, -t45 * t31, t35 ^ 2 * t49, t35 * t29 * t38, t35 * t44, t38 * t44, qJD(4) ^ 2 / 0.2e1, -t18 * t47 + (-t35 * t19 + t28) * qJD(4), -t46 * qJD(4) + t18 * t48, t21 ^ 2 / 0.2e1, -t21 * t20, t21 * t30, -t20 * t30, t30 ^ 2 / 0.2e1, t16 * t20 + (t37 * t14 - t34 * t15) * t30, t16 * t21 - (t34 * t14 + t37 * t15) * t30;];
T_reg = t1;
