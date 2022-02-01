% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:05
% EndTime: 2022-01-20 12:02:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (309->32), mult. (412->93), div. (0->0), fcn. (209->8), ass. (0->38)
t26 = sin(qJ(4));
t23 = t26 ^ 2;
t44 = t23 / 0.2e1;
t29 = cos(qJ(4));
t24 = t29 ^ 2;
t43 = t24 / 0.2e1;
t42 = cos(qJ(5));
t22 = qJD(1) + qJD(2);
t20 = qJD(3) + t22;
t41 = t20 * t26;
t40 = t20 * t29;
t31 = cos(qJ(2));
t39 = pkin(1) * qJD(1);
t35 = t31 * t39;
t15 = t22 * pkin(2) + t35;
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t28 = sin(qJ(2));
t36 = t28 * t39;
t10 = t27 * t15 + t30 * t36;
t38 = qJD(4) * t26;
t37 = qJD(4) * t29;
t8 = t20 * pkin(8) + t10;
t34 = pkin(9) * t20 + t8;
t9 = t30 * t15 - t27 * t36;
t32 = qJD(1) ^ 2;
t25 = sin(qJ(5));
t21 = qJD(4) + qJD(5);
t19 = t20 ^ 2;
t13 = (t25 * t29 + t42 * t26) * t20;
t11 = t25 * t41 - t42 * t40;
t7 = -t20 * pkin(3) - t9;
t5 = (-pkin(4) * t29 - pkin(3)) * t20 - t9;
t4 = t34 * t29;
t3 = qJD(4) * pkin(4) - t34 * t26;
t2 = t25 * t3 + t42 * t4;
t1 = -t25 * t4 + t42 * t3;
t6 = [0, 0, 0, 0, 0, t32 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 ^ 2 / 0.2e1, t22 * t35, -t22 * t36, 0, (t28 ^ 2 / 0.2e1 + t31 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t32, 0, 0, 0, 0, 0, t19 / 0.2e1, t9 * t20, -t10 * t20, 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t19 * t44, t26 * t19 * t29, t20 * t38, t19 * t43, t20 * t37, qJD(4) ^ 2 / 0.2e1, -t8 * t38 - t7 * t40, -t8 * t37 + t7 * t41, (t23 + t24) * t8 * t20, t7 ^ 2 / 0.2e1 + (t43 + t44) * t8 ^ 2, t13 ^ 2 / 0.2e1, -t13 * t11, t13 * t21, t11 ^ 2 / 0.2e1, -t11 * t21, t21 ^ 2 / 0.2e1, t1 * t21 + t5 * t11, t5 * t13 - t2 * t21, -t1 * t13 - t2 * t11, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
