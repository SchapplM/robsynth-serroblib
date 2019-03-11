% Calculate inertial parameters regressor of fixed base kinetic energy for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% T_reg [1x(2*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S2RR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_energykin_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_energykin_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_energykin_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:59
% EndTime: 2019-03-08 18:00:59
% DurationCPUTime: 0.06s
% Computational Cost: add. (4->4), mult. (45->22), div. (0->0), fcn. (12->2), ass. (0->10)
t5 = qJD(1) ^ 2;
t10 = t5 / 0.2e1;
t9 = qJD(1) * qJD(2);
t3 = sin(qJ(2));
t8 = t3 * t9;
t4 = cos(qJ(2));
t7 = t4 * t9;
t2 = t4 ^ 2;
t1 = t3 ^ 2;
t6 = [0, 0, 0, 0, 0, t10, 0, 0, 0, 0, t1 * t10, t3 * t5 * t4, t8, t2 * t10, t7, qJD(2) ^ 2 / 0.2e1, -pkin(1) * t8, -pkin(1) * t7 (t1 + t2) * t5 * pkin(1) (t2 / 0.2e1 + t1 / 0.2e1) * pkin(1) ^ 2 * t5;];
T_reg  = t6;
