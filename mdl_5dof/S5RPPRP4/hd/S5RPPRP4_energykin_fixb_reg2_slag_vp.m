% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:19
% EndTime: 2019-12-31 17:52:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (132->34), mult. (264->69), div. (0->0), fcn. (97->4), ass. (0->30)
t28 = qJD(1) ^ 2;
t22 = t28 / 0.2e1;
t26 = sin(qJ(4));
t27 = cos(qJ(4));
t12 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t23 = sin(pkin(7));
t24 = cos(pkin(7));
t33 = qJ(2) * qJD(1);
t10 = t23 * t12 + t24 * t33;
t8 = -qJD(1) * pkin(6) + t10;
t4 = t26 * qJD(3) + t27 * t8;
t35 = qJD(1) * t26;
t34 = qJD(1) * t27;
t32 = qJ(5) * qJD(1);
t31 = qJD(1) * qJD(4);
t30 = t26 * t31;
t29 = t27 * t31;
t9 = t24 * t12 - t23 * t33;
t7 = qJD(1) * pkin(3) - t9;
t21 = qJD(4) ^ 2 / 0.2e1;
t20 = t27 * qJD(3);
t18 = -qJD(1) * pkin(1) + qJD(2);
t17 = t27 ^ 2 * t22;
t16 = t26 ^ 2 * t22;
t13 = t26 * t28 * t27;
t5 = pkin(4) * t34 + qJD(5) + t7;
t3 = -t26 * t8 + t20;
t2 = -t27 * t32 + t4;
t1 = qJD(4) * pkin(4) + t20 + (-t8 + t32) * t26;
t6 = [0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, -t18 * qJD(1), 0, t28 * qJ(2), qJ(2) ^ 2 * t22 + t18 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t22, -t9 * qJD(1), t10 * qJD(1), 0, t10 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t16, t13, -t30, t17, -t29, t21, t3 * qJD(4) + t7 * t34, -t4 * qJD(4) - t7 * t35, (t26 * t3 - t27 * t4) * qJD(1), t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1, t16, t13, -t30, t17, -t29, t21, t1 * qJD(4) + t5 * t34, -t2 * qJD(4) - t5 * t35, (t1 * t26 - t2 * t27) * qJD(1), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1;];
T_reg = t6;
