% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:38
% EndTime: 2019-12-31 16:31:38
% DurationCPUTime: 0.08s
% Computational Cost: add. (44->13), mult. (100->48), div. (0->0), fcn. (35->4), ass. (0->19)
t6 = qJD(2) + qJD(3);
t5 = t6 ^ 2;
t19 = t5 / 0.2e1;
t11 = cos(qJ(3));
t17 = pkin(2) * qJD(2);
t14 = t11 * t17;
t4 = -t6 * pkin(3) - t14;
t18 = t4 * t6;
t16 = qJD(4) * t6;
t9 = sin(qJ(3));
t15 = t9 * t17;
t12 = qJD(2) ^ 2;
t10 = cos(qJ(4));
t8 = sin(qJ(4));
t7 = qJD(1) ^ 2 / 0.2e1;
t3 = t6 * pkin(6) + t15;
t2 = t8 * qJD(1) + t10 * t3;
t1 = t10 * qJD(1) - t8 * t3;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, t12 / 0.2e1, 0, 0, 0, t7, 0, 0, 0, 0, 0, t19, t6 * t14, -t6 * t15, 0, t7 + (t9 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t12, t8 ^ 2 * t19, t8 * t5 * t10, t8 * t16, t10 ^ 2 * t19, t10 * t16, qJD(4) ^ 2 / 0.2e1, qJD(4) * t1 - t10 * t18, -qJD(4) * t2 + t8 * t18, (-t1 * t8 + t10 * t2) * t6, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg = t13;
