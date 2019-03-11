% Calculate minimal parameter regressor of fixed base kinetic energy for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% 
% Output:
% T_reg [1x9]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S3RPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_energykin_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_energykin_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_energykin_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:01
% EndTime: 2019-03-08 18:05:01
% DurationCPUTime: 0.02s
% Computational Cost: add. (11->7), mult. (25->17), div. (0->0), fcn. (0->0), ass. (0->6)
t20 = qJD(1) ^ 2;
t21 = t20 / 0.2e1;
t19 = -qJD(1) * pkin(1) + qJD(2);
t18 = qJD(1) * qJ(2) + qJD(3);
t17 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t1 = [t21, 0, 0, t19 * qJD(1), t20 * qJ(2), qJ(2) ^ 2 * t21 + t19 ^ 2 / 0.2e1, t18 * qJD(1), -t17 * qJD(1), t17 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1;];
T_reg  = t1;
