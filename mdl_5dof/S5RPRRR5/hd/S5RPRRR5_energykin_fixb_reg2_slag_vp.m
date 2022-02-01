% Calculate inertial parameters regressor of fixed base kinetic energy for
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
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:57
% EndTime: 2022-01-20 09:48:57
% DurationCPUTime: 0.16s
% Computational Cost: add. (227->37), mult. (432->96), div. (0->0), fcn. (229->8), ass. (0->36)
t25 = qJD(1) + qJD(3);
t23 = t25 ^ 2;
t43 = t23 / 0.2e1;
t35 = qJD(1) ^ 2;
t42 = pkin(1) * t35;
t29 = cos(pkin(9));
t17 = (pkin(1) * t29 + pkin(2)) * qJD(1);
t32 = sin(qJ(3));
t34 = cos(qJ(3));
t28 = sin(pkin(9));
t37 = pkin(1) * qJD(1) * t28;
t12 = t32 * t17 + t34 * t37;
t10 = t25 * pkin(7) + t12;
t31 = sin(qJ(4));
t33 = cos(qJ(4));
t6 = t31 * qJD(2) + t33 * t10;
t41 = cos(qJ(5));
t40 = t25 * t31;
t39 = t25 * t33;
t38 = qJD(4) * t25;
t11 = t34 * t17 - t32 * t37;
t30 = sin(qJ(5));
t27 = t35 / 0.2e1;
t26 = qJD(2) ^ 2 / 0.2e1;
t24 = qJD(4) + qJD(5);
t22 = t33 * qJD(2);
t15 = (t30 * t33 + t31 * t41) * t25;
t13 = t30 * t40 - t39 * t41;
t9 = -t25 * pkin(3) - t11;
t7 = (-pkin(4) * t33 - pkin(3)) * t25 - t11;
t5 = -t31 * t10 + t22;
t4 = pkin(8) * t39 + t6;
t3 = qJD(4) * pkin(4) + t22 + (-pkin(8) * t25 - t10) * t31;
t2 = t30 * t3 + t4 * t41;
t1 = t3 * t41 - t30 * t4;
t8 = [0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29 * t42, -t28 * t42, 0, t26 + (t28 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t35, 0, 0, 0, 0, 0, t43, t11 * t25, -t12 * t25, 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t26, t31 ^ 2 * t43, t31 * t23 * t33, t31 * t38, t33 ^ 2 * t43, t33 * t38, qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) - t39 * t9, -t6 * qJD(4) + t40 * t9, (-t31 * t5 + t33 * t6) * t25, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t24, t13 ^ 2 / 0.2e1, -t13 * t24, t24 ^ 2 / 0.2e1, t1 * t24 + t7 * t13, t7 * t15 - t2 * t24, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
