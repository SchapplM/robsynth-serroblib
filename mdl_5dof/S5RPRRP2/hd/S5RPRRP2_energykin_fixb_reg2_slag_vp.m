% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRP2
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
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:09
% EndTime: 2022-01-23 09:28:09
% DurationCPUTime: 0.14s
% Computational Cost: add. (153->31), mult. (317->75), div. (0->0), fcn. (145->6), ass. (0->36)
t23 = qJD(1) + qJD(3);
t22 = t23 ^ 2;
t41 = t22 / 0.2e1;
t33 = qJD(1) ^ 2;
t40 = pkin(1) * t33;
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t28 = cos(pkin(8));
t12 = (pkin(1) * t28 + pkin(2)) * qJD(1);
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t27 = sin(pkin(8));
t35 = pkin(1) * qJD(1) * t27;
t10 = t30 * t12 + t32 * t35;
t8 = t23 * pkin(7) + t10;
t5 = t29 * qJD(2) + t31 * t8;
t39 = t23 * t29;
t38 = t23 * t31;
t37 = qJ(5) * t23;
t36 = qJD(4) * t23;
t9 = t32 * t12 - t30 * t35;
t26 = t33 / 0.2e1;
t25 = qJD(2) ^ 2 / 0.2e1;
t24 = qJD(4) ^ 2 / 0.2e1;
t21 = t31 * qJD(2);
t19 = t31 * t36;
t18 = t29 * t36;
t17 = t31 ^ 2 * t41;
t16 = t29 ^ 2 * t41;
t13 = t29 * t22 * t31;
t7 = -t23 * pkin(3) - t9;
t4 = -t29 * t8 + t21;
t3 = qJD(5) + (-pkin(4) * t31 - pkin(3)) * t23 - t9;
t2 = t31 * t37 + t5;
t1 = qJD(4) * pkin(4) + t21 + (-t8 - t37) * t29;
t6 = [0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t28 * t40, -t27 * t40, 0, t25 + (t27 ^ 2 / 0.2e1 + t28 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t33, 0, 0, 0, 0, 0, t41, t9 * t23, -t10 * t23, 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t25, t16, t13, t18, t17, t19, t24, t4 * qJD(4) - t7 * t38, -t5 * qJD(4) + t7 * t39, (-t29 * t4 + t31 * t5) * t23, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t16, t13, t18, t17, t19, t24, t1 * qJD(4) - t3 * t38, -t2 * qJD(4) + t3 * t39, (-t1 * t29 + t2 * t31) * t23, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t6;
