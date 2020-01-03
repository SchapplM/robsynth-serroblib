% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PPRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:26
% EndTime: 2019-12-31 16:17:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (24->13), mult. (85->46), div. (0->0), fcn. (35->4), ass. (0->16)
t10 = qJD(3) ^ 2;
t15 = t10 / 0.2e1;
t9 = cos(qJ(3));
t4 = -qJD(3) * pkin(3) - t9 * qJD(2);
t14 = qJD(3) * t4;
t13 = qJD(2) * qJD(3);
t12 = qJD(3) * qJD(4);
t11 = qJD(2) ^ 2;
t8 = cos(qJ(4));
t7 = sin(qJ(3));
t6 = sin(qJ(4));
t5 = qJD(1) ^ 2 / 0.2e1;
t3 = qJD(3) * pkin(5) + t7 * qJD(2);
t2 = -t6 * qJD(1) + t8 * t3;
t1 = -t8 * qJD(1) - t6 * t3;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 + t11 / 0.2e1, 0, 0, 0, 0, 0, t15, t9 * t13, -t7 * t13, 0, t5 + (t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1) * t11, t6 ^ 2 * t15, t6 * t10 * t8, t6 * t12, t8 ^ 2 * t15, t8 * t12, qJD(4) ^ 2 / 0.2e1, t1 * qJD(4) - t8 * t14, -t2 * qJD(4) + t6 * t14, (-t1 * t6 + t2 * t8) * qJD(3), t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1;];
T_reg = t16;
