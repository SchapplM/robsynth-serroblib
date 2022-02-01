% Calculate inertial parameters regressor of fixed base kinetic energy for
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
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:23
% EndTime: 2022-01-20 11:49:23
% DurationCPUTime: 0.15s
% Computational Cost: add. (306->38), mult. (470->93), div. (0->0), fcn. (251->6), ass. (0->41)
t31 = sin(qJ(3));
t28 = t31 ^ 2;
t47 = t28 / 0.2e1;
t33 = cos(qJ(3));
t29 = t33 ^ 2;
t46 = t29 / 0.2e1;
t27 = qJD(1) + qJD(2);
t32 = sin(qJ(2));
t42 = pkin(1) * qJD(1);
t39 = t32 * t42;
t21 = t27 * pkin(7) + t39;
t37 = pkin(8) * t27 + t21;
t10 = qJD(3) * pkin(3) - t37 * t31;
t13 = t37 * t33;
t30 = sin(qJ(4));
t45 = cos(qJ(4));
t4 = t30 * t10 + t45 * t13;
t44 = t27 * t31;
t43 = t27 * t33;
t41 = qJD(3) * t31;
t40 = qJD(3) * t33;
t34 = cos(qJ(2));
t38 = t34 * t42;
t3 = t45 * t10 - t30 * t13;
t19 = -t38 + (-pkin(3) * t33 - pkin(2)) * t27;
t35 = qJD(1) ^ 2;
t26 = qJD(3) + qJD(4);
t25 = t27 ^ 2;
t24 = t26 ^ 2 / 0.2e1;
t22 = -t27 * pkin(2) - t38;
t18 = (t30 * t33 + t45 * t31) * t27;
t16 = t30 * t44 - t45 * t43;
t15 = t18 ^ 2 / 0.2e1;
t14 = t16 ^ 2 / 0.2e1;
t12 = t18 * t26;
t11 = t16 * t26;
t6 = t18 * t16;
t5 = t16 * pkin(4) + qJD(5) + t19;
t2 = -t16 * qJ(5) + t4;
t1 = t26 * pkin(4) - t18 * qJ(5) + t3;
t7 = [0, 0, 0, 0, 0, t35 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 / 0.2e1, t27 * t38, -t27 * t39, 0, (t32 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t35, t25 * t47, t31 * t25 * t33, t27 * t41, t25 * t46, t27 * t40, qJD(3) ^ 2 / 0.2e1, -t21 * t41 - t22 * t43, -t21 * t40 + t22 * t44, (t28 + t29) * t27 * t21, t22 ^ 2 / 0.2e1 + (t46 + t47) * t21 ^ 2, t15, -t6, t12, t14, -t11, t24, t19 * t16 + t3 * t26, t19 * t18 - t4 * t26, -t4 * t16 - t3 * t18, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t19 ^ 2 / 0.2e1, t15, -t6, t12, t14, -t11, t24, t1 * t26 + t5 * t16, t5 * t18 - t2 * t26, -t1 * t18 - t2 * t16, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t7;
