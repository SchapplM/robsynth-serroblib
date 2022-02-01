% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:27
% EndTime: 2022-01-20 11:49:27
% DurationCPUTime: 0.07s
% Computational Cost: add. (212->29), mult. (303->69), div. (0->0), fcn. (169->6), ass. (0->29)
t27 = qJD(1) + qJD(2);
t25 = t27 ^ 2;
t44 = t25 / 0.2e1;
t43 = cos(qJ(4));
t29 = sin(qJ(3));
t42 = t27 * t29;
t31 = cos(qJ(3));
t41 = t27 * t31;
t39 = pkin(1) * qJD(1);
t36 = sin(qJ(2)) * t39;
t22 = t27 * pkin(7) + t36;
t34 = pkin(8) * t27 + t22;
t17 = qJD(3) * pkin(3) - t34 * t29;
t18 = t34 * t31;
t28 = sin(qJ(4));
t40 = t28 * t17 + t43 * t18;
t38 = qJD(3) * t29;
t37 = qJD(3) * t31;
t35 = cos(qJ(2)) * t39;
t33 = t43 * t17 - t28 * t18;
t21 = -t35 + (-pkin(3) * t31 - pkin(2)) * t27;
t26 = qJD(3) + qJD(4);
t23 = -t27 * pkin(2) - t35;
t20 = (t28 * t31 + t43 * t29) * t27;
t19 = t28 * t42 - t43 * t41;
t13 = t19 * pkin(4) + qJD(5) + t21;
t12 = -t19 * qJ(5) + t40;
t11 = t26 * pkin(4) - t20 * qJ(5) + t33;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t44, t27 * t35, -t27 * t36, t29 ^ 2 * t44, t29 * t25 * t31, t27 * t38, t27 * t37, qJD(3) ^ 2 / 0.2e1, -t22 * t38 - t23 * t41, -t22 * t37 + t23 * t42, t20 ^ 2 / 0.2e1, -t20 * t19, t20 * t26, -t19 * t26, t26 ^ 2 / 0.2e1, t21 * t19 + t33 * t26, t21 * t20 - t40 * t26, t11 * t26 + t13 * t19, -t12 * t26 + t13 * t20, -t11 * t20 - t12 * t19, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg = t1;
