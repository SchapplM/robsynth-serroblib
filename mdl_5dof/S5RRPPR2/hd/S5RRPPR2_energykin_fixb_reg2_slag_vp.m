% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:46
% EndTime: 2022-01-20 10:05:46
% DurationCPUTime: 0.13s
% Computational Cost: add. (222->31), mult. (373->83), div. (0->0), fcn. (190->8), ass. (0->35)
t21 = qJD(1) + qJD(2);
t19 = t21 ^ 2;
t17 = t19 / 0.2e1;
t22 = sin(pkin(9));
t42 = t22 ^ 2 * t19;
t41 = t21 * t22;
t24 = cos(pkin(9));
t40 = t24 * t21;
t29 = cos(qJ(2));
t39 = pkin(1) * qJD(1);
t34 = t29 * t39;
t12 = t21 * pkin(2) + t34;
t23 = sin(pkin(8));
t25 = cos(pkin(8));
t27 = sin(qJ(2));
t35 = t27 * t39;
t10 = t23 * t12 + t25 * t35;
t8 = t21 * qJ(4) + t10;
t4 = -t24 * qJD(3) + t22 * t8;
t38 = t4 ^ 2 / 0.2e1;
t26 = sin(qJ(5));
t37 = t26 * t41;
t28 = cos(qJ(5));
t36 = t28 * t41;
t33 = t42 / 0.2e1;
t9 = t25 * t12 - t23 * t35;
t32 = qJD(4) - t9;
t30 = qJD(1) ^ 2;
t13 = -qJD(5) + t40;
t7 = -t21 * pkin(3) + t32;
t6 = t22 * qJD(3) + t24 * t8;
t3 = (-pkin(4) * t24 - pkin(7) * t22 - pkin(3)) * t21 + t32;
t2 = t26 * t3 + t28 * t6;
t1 = -t26 * t6 + t28 * t3;
t5 = [0, 0, 0, 0, 0, t30 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t21 * t34, -t21 * t35, 0, (t27 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t30, 0, 0, 0, 0, 0, t17, t9 * t21, -t10 * t21, 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t33, t22 * t19 * t24, 0, t24 ^ 2 * t17, 0, 0, -t7 * t40, t7 * t41, (t22 * t4 + t24 * t6) * t21, t6 ^ 2 / 0.2e1 + t38 + t7 ^ 2 / 0.2e1, t28 ^ 2 * t33, -t28 * t26 * t42, -t13 * t36, t26 ^ 2 * t33, t13 * t37, t13 ^ 2 / 0.2e1, -t1 * t13 + t4 * t37, t2 * t13 + t4 * t36, (-t1 * t28 - t2 * t26) * t41, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t38;];
T_reg = t5;
