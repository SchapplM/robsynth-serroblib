% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:55
% EndTime: 2022-01-20 10:19:55
% DurationCPUTime: 0.12s
% Computational Cost: add. (177->30), mult. (315->76), div. (0->0), fcn. (145->6), ass. (0->35)
t24 = qJD(1) + qJD(2);
t23 = t24 ^ 2;
t20 = t23 / 0.2e1;
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t31 = cos(qJ(2));
t38 = pkin(1) * qJD(1);
t34 = t31 * t38;
t12 = t24 * pkin(2) + t34;
t26 = sin(pkin(8));
t27 = cos(pkin(8));
t29 = sin(qJ(2));
t35 = t29 * t38;
t10 = t26 * t12 + t27 * t35;
t8 = t24 * pkin(7) + t10;
t5 = t28 * qJD(3) + t30 * t8;
t40 = t24 * t28;
t39 = t24 * t30;
t37 = qJ(5) * t24;
t36 = qJD(4) * t24;
t9 = t27 * t12 - t26 * t35;
t32 = qJD(1) ^ 2;
t25 = qJD(4) ^ 2 / 0.2e1;
t22 = t30 * qJD(3);
t19 = t30 * t36;
t18 = t28 * t36;
t17 = t30 ^ 2 * t20;
t16 = t28 ^ 2 * t20;
t13 = t28 * t23 * t30;
t7 = -t24 * pkin(3) - t9;
t4 = -t28 * t8 + t22;
t3 = qJD(5) + (-pkin(4) * t30 - pkin(3)) * t24 - t9;
t2 = t30 * t37 + t5;
t1 = qJD(4) * pkin(4) + t22 + (-t8 - t37) * t28;
t6 = [0, 0, 0, 0, 0, t32 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t24 * t34, -t24 * t35, 0, (t29 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t32, 0, 0, 0, 0, 0, t20, t9 * t24, -t10 * t24, 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t16, t13, t18, t17, t19, t25, t4 * qJD(4) - t7 * t39, -t5 * qJD(4) + t7 * t40, (-t28 * t4 + t30 * t5) * t24, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t16, t13, t18, t17, t19, t25, t1 * qJD(4) - t3 * t39, -t2 * qJD(4) + t3 * t40, (-t1 * t28 + t2 * t30) * t24, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t6;
