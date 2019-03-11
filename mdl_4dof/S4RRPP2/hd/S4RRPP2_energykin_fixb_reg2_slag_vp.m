% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPP2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:05
% EndTime: 2019-03-08 18:34:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->13), mult. (66->28), div. (0->0), fcn. (12->2), ass. (0->15)
t15 = pkin(1) * qJD(1);
t8 = sin(qJ(2));
t14 = t8 * t15;
t6 = qJD(1) + qJD(2);
t4 = t6 * qJ(3) + t14;
t16 = t4 * t6;
t5 = t6 ^ 2 / 0.2e1;
t9 = cos(qJ(2));
t13 = t9 * t15;
t12 = qJD(3) - t13;
t10 = qJD(1) ^ 2;
t3 = -t6 * pkin(2) + t12;
t2 = t4 ^ 2 / 0.2e1;
t1 = (-pkin(2) - pkin(3)) * t6 + t12;
t7 = [0, 0, 0, 0, 0, t10 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6 * t13, -t6 * t14, 0 (t8 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t10, 0, 0, 0, t5, 0, 0, -t3 * t6, 0, t16, t2 + t3 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t5, -t1 * t6, t16, 0, t2 + t1 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t7;
