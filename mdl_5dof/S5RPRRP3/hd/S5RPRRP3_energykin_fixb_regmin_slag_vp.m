% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:30:18
% EndTime: 2022-01-23 09:30:18
% DurationCPUTime: 0.06s
% Computational Cost: add. (136->32), mult. (327->76), div. (0->0), fcn. (187->6), ass. (0->28)
t39 = qJD(1) ^ 2;
t49 = t39 / 0.2e1;
t48 = cos(qJ(4));
t34 = sin(pkin(8));
t28 = (pkin(1) * t34 + pkin(6)) * qJD(1);
t38 = cos(qJ(3));
t32 = t38 * qJD(2);
t37 = sin(qJ(3));
t22 = qJD(3) * pkin(3) + t32 + (-pkin(7) * qJD(1) - t28) * t37;
t44 = qJD(1) * t38;
t46 = t37 * qJD(2) + t38 * t28;
t23 = pkin(7) * t44 + t46;
t36 = sin(qJ(4));
t47 = t36 * t22 + t48 * t23;
t45 = qJD(1) * t37;
t43 = qJD(1) * qJD(3);
t35 = cos(pkin(8));
t42 = -pkin(1) * t35 - pkin(2);
t41 = t48 * t22 - t36 * t23;
t26 = (-pkin(3) * t38 + t42) * qJD(1);
t33 = qJD(3) + qJD(4);
t29 = t42 * qJD(1);
t25 = (t36 * t38 + t48 * t37) * qJD(1);
t24 = t36 * t45 - t48 * t44;
t18 = t24 * pkin(4) + qJD(5) + t26;
t17 = -t24 * qJ(5) + t47;
t16 = t33 * pkin(4) - t25 * qJ(5) + t41;
t1 = [t49, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t34 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t39, t37 ^ 2 * t49, t37 * t39 * t38, t37 * t43, t38 * t43, qJD(3) ^ 2 / 0.2e1, -t29 * t44 + (-t37 * t28 + t32) * qJD(3), -t46 * qJD(3) + t29 * t45, t25 ^ 2 / 0.2e1, -t25 * t24, t25 * t33, -t24 * t33, t33 ^ 2 / 0.2e1, t26 * t24 + t41 * t33, t26 * t25 - t47 * t33, t16 * t33 + t18 * t24, -t17 * t33 + t18 * t25, -t16 * t25 - t17 * t24, t17 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1;];
T_reg = t1;
