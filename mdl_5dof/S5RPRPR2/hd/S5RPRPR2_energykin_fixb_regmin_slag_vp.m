% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:06
% EndTime: 2022-01-23 09:19:06
% DurationCPUTime: 0.14s
% Computational Cost: add. (114->29), mult. (213->65), div. (0->0), fcn. (115->8), ass. (0->27)
t25 = qJD(1) + qJD(3);
t28 = cos(pkin(9));
t40 = t25 * t28;
t29 = cos(pkin(8));
t20 = (pkin(1) * t29 + pkin(2)) * qJD(1);
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t27 = sin(pkin(8));
t38 = pkin(1) * qJD(1) * t27;
t39 = t31 * t20 + t33 * t38;
t15 = t25 * qJ(4) + t39;
t26 = sin(pkin(9));
t11 = t26 * qJD(2) + t28 * t15;
t37 = t33 * t20 - t31 * t38;
t36 = qJD(4) - t37;
t34 = qJD(1) ^ 2;
t32 = cos(qJ(5));
t30 = sin(qJ(5));
t24 = t28 * qJD(2);
t17 = (t26 * t32 + t28 * t30) * t25;
t16 = t30 * t26 * t25 - t32 * t40;
t14 = -t25 * pkin(3) + t36;
t12 = (-pkin(4) * t28 - pkin(3)) * t25 + t36;
t10 = -t26 * t15 + t24;
t9 = pkin(7) * t40 + t11;
t8 = t24 + (-pkin(7) * t25 - t15) * t26;
t1 = [t34 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t27 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t34, t25 ^ 2 / 0.2e1, t37 * t25, -t39 * t25, -t14 * t40, (-t10 * t26 + t11 * t28) * t25, t11 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t17 ^ 2 / 0.2e1, -t17 * t16, t17 * qJD(5), -t16 * qJD(5), qJD(5) ^ 2 / 0.2e1, t12 * t16 + (-t30 * t9 + t32 * t8) * qJD(5), t12 * t17 - (t30 * t8 + t32 * t9) * qJD(5);];
T_reg = t1;
