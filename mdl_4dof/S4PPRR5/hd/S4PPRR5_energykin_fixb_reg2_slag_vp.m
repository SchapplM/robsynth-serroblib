% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PPRR5_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:51
% EndTime: 2019-12-31 16:19:51
% DurationCPUTime: 0.07s
% Computational Cost: add. (30->12), mult. (90->42), div. (0->0), fcn. (40->4), ass. (0->17)
t15 = qJD(3) ^ 2;
t19 = t15 / 0.2e1;
t12 = sin(qJ(3));
t14 = cos(qJ(3));
t5 = t14 * qJD(1) + t12 * qJD(2);
t4 = -t12 * qJD(1) + t14 * qJD(2);
t2 = -qJD(3) * pkin(3) - t4;
t18 = qJD(3) * t2;
t3 = qJD(3) * pkin(5) + t5;
t17 = qJD(4) * t3;
t16 = qJD(3) * qJD(4);
t13 = cos(qJ(4));
t11 = sin(qJ(4));
t10 = t13 ^ 2;
t9 = t11 ^ 2;
t8 = qJD(1) ^ 2 / 0.2e1;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 + qJD(2) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t19, t4 * qJD(3), -t5 * qJD(3), 0, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, t9 * t19, t11 * t15 * t13, t11 * t16, t10 * t19, t13 * t16, qJD(4) ^ 2 / 0.2e1, -t11 * t17 - t13 * t18, t11 * t18 - t13 * t17, (t10 + t9) * t3 * qJD(3), t2 ^ 2 / 0.2e1 + (t10 / 0.2e1 + t9 / 0.2e1) * t3 ^ 2;];
T_reg = t1;
