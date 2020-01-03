% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR8_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:22
% EndTime: 2019-12-31 17:08:22
% DurationCPUTime: 0.11s
% Computational Cost: add. (98->34), mult. (281->84), div. (0->0), fcn. (121->4), ass. (0->33)
t27 = qJD(1) ^ 2;
t38 = t27 / 0.2e1;
t37 = pkin(2) + pkin(3);
t26 = cos(qJ(2));
t36 = t26 * t27;
t34 = qJD(1) * t26;
t11 = pkin(5) * t34 + qJD(2) * qJ(3);
t24 = sin(qJ(2));
t35 = qJD(1) * t24;
t33 = pkin(5) * t35 + qJD(3);
t32 = qJD(1) * qJD(2);
t31 = t24 * t36;
t14 = t24 * t32;
t30 = t26 * t32;
t29 = qJ(3) * t24 + pkin(1);
t25 = cos(qJ(4));
t23 = sin(qJ(4));
t22 = t26 ^ 2;
t21 = t24 ^ 2;
t19 = qJD(2) ^ 2 / 0.2e1;
t17 = qJD(2) - qJD(4);
t13 = t22 * t38;
t12 = t21 * t38;
t10 = -qJD(2) * pkin(2) + t33;
t9 = (-pkin(2) * t26 - t29) * qJD(1);
t8 = -pkin(6) * t34 + t11;
t7 = (-t23 * t26 + t24 * t25) * qJD(1);
t5 = (-t23 * t24 - t25 * t26) * qJD(1);
t4 = -pkin(6) * t35 - t37 * qJD(2) + t33;
t3 = (t37 * t26 + t29) * qJD(1);
t2 = t23 * t4 + t25 * t8;
t1 = -t23 * t8 + t25 * t4;
t6 = [0, 0, 0, 0, 0, t38, 0, 0, 0, 0, t12, t31, t14, t13, t30, t19, pkin(1) * t36 - pkin(5) * t14, -t27 * pkin(1) * t24 - pkin(5) * t30, (t21 + t22) * t27 * pkin(5), (pkin(1) ^ 2 / 0.2e1 + (t22 / 0.2e1 + t21 / 0.2e1) * pkin(5) ^ 2) * t27, t12, t14, -t31, t19, -t30, t13, -t10 * qJD(2) - t9 * t34, (t10 * t24 + t11 * t26) * qJD(1), t11 * qJD(2) - t9 * t35, t11 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t7 ^ 2 / 0.2e1, t7 * t5, -t7 * t17, t5 ^ 2 / 0.2e1, -t5 * t17, t17 ^ 2 / 0.2e1, -t1 * t17 - t3 * t5, t2 * t17 + t3 * t7, -t1 * t7 + t2 * t5, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg = t6;
