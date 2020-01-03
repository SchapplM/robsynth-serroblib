% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:53
% EndTime: 2019-12-31 18:12:53
% DurationCPUTime: 0.14s
% Computational Cost: add. (190->42), mult. (552->85), div. (0->0), fcn. (331->4), ass. (0->35)
t31 = qJD(1) ^ 2;
t42 = t31 / 0.2e1;
t41 = cos(qJ(3));
t30 = sin(qJ(3));
t29 = cos(pkin(7));
t36 = qJD(1) * t29;
t28 = sin(pkin(7));
t37 = qJD(1) * t28;
t16 = t30 * t37 - t41 * t36;
t18 = (-t41 * t28 - t29 * t30) * qJD(1);
t40 = t16 * t18;
t39 = pkin(3) + qJ(5);
t38 = pkin(6) + qJ(2);
t20 = t38 * t37;
t21 = t38 * t36;
t8 = -t30 * t20 + t41 * t21;
t12 = qJD(3) * t16;
t35 = t18 * qJD(3);
t9 = t16 ^ 2 / 0.2e1;
t10 = t18 ^ 2 / 0.2e1;
t7 = -t41 * t20 - t30 * t21;
t6 = -qJD(3) * qJ(4) - t8;
t34 = qJD(4) - t7;
t22 = qJD(2) + (-pkin(2) * t29 - pkin(1)) * qJD(1);
t33 = t18 * qJ(4) + t22;
t27 = qJD(3) ^ 2 / 0.2e1;
t26 = t29 ^ 2;
t25 = t28 ^ 2;
t24 = -qJD(1) * pkin(1) + qJD(2);
t5 = -qJD(3) * pkin(3) + t34;
t4 = t16 * pkin(3) + t33;
t3 = -t16 * pkin(4) + qJD(5) - t6;
t2 = -t18 * pkin(4) - t39 * qJD(3) + t34;
t1 = t39 * t16 + t33;
t11 = [0, 0, 0, 0, 0, t42, 0, 0, 0, 0, t25 * t42, t28 * t31 * t29, 0, t26 * t42, 0, 0, -t24 * t36, t24 * t37, (t25 + t26) * t31 * qJ(2), t24 ^ 2 / 0.2e1 + (t26 / 0.2e1 + t25 / 0.2e1) * qJ(2) ^ 2 * t31, t10, t40, -t35, t9, -t12, t27, t7 * qJD(3) + t22 * t16, -t8 * qJD(3) - t22 * t18, -t8 * t16 + t7 * t18, t8 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1, t27, t35, t12, t10, t40, t9, t6 * t16 - t5 * t18, t5 * qJD(3) - t4 * t16, -t6 * qJD(3) + t4 * t18, t4 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1, t27, t12, -t35, t9, -t40, t10, -t3 * t16 - t2 * t18, t3 * qJD(3) + t1 * t18, -t2 * qJD(3) + t1 * t16, t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t11;
