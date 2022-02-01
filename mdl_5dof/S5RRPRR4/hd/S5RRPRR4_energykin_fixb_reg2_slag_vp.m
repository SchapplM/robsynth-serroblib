% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRPRR4
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
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:27
% EndTime: 2022-01-20 10:48:27
% DurationCPUTime: 0.17s
% Computational Cost: add. (257->36), mult. (430->97), div. (0->0), fcn. (229->8), ass. (0->35)
t26 = qJD(1) + qJD(2);
t24 = t26 ^ 2;
t21 = t24 / 0.2e1;
t33 = cos(qJ(2));
t39 = pkin(1) * qJD(1);
t36 = t33 * t39;
t17 = t26 * pkin(2) + t36;
t27 = sin(pkin(9));
t28 = cos(pkin(9));
t31 = sin(qJ(2));
t37 = t31 * t39;
t12 = t27 * t17 + t28 * t37;
t10 = t26 * pkin(7) + t12;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t6 = t30 * qJD(3) + t32 * t10;
t42 = cos(qJ(5));
t41 = t26 * t30;
t40 = t26 * t32;
t38 = qJD(4) * t26;
t11 = t28 * t17 - t27 * t37;
t34 = qJD(1) ^ 2;
t29 = sin(qJ(5));
t25 = qJD(4) + qJD(5);
t23 = t32 * qJD(3);
t15 = (t29 * t32 + t42 * t30) * t26;
t13 = t29 * t41 - t42 * t40;
t9 = -t26 * pkin(3) - t11;
t7 = (-pkin(4) * t32 - pkin(3)) * t26 - t11;
t5 = -t30 * t10 + t23;
t4 = pkin(8) * t40 + t6;
t3 = qJD(4) * pkin(4) + t23 + (-pkin(8) * t26 - t10) * t30;
t2 = t29 * t3 + t42 * t4;
t1 = -t29 * t4 + t42 * t3;
t8 = [0, 0, 0, 0, 0, t34 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t26 * t36, -t26 * t37, 0, (t31 ^ 2 / 0.2e1 + t33 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t34, 0, 0, 0, 0, 0, t21, t11 * t26, -t12 * t26, 0, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t30 ^ 2 * t21, t30 * t24 * t32, t30 * t38, t32 ^ 2 * t21, t32 * t38, qJD(4) ^ 2 / 0.2e1, t5 * qJD(4) - t40 * t9, -t6 * qJD(4) + t9 * t41, (-t30 * t5 + t32 * t6) * t26, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t15 ^ 2 / 0.2e1, -t15 * t13, t15 * t25, t13 ^ 2 / 0.2e1, -t13 * t25, t25 ^ 2 / 0.2e1, t1 * t25 + t7 * t13, t7 * t15 - t2 * t25, -t1 * t15 - t2 * t13, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t8;
