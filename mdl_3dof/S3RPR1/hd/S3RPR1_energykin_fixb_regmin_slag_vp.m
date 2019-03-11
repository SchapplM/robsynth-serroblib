% Calculate minimal parameter regressor of fixed base kinetic energy for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% T_reg [1x9]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S3RPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_energykin_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_energykin_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_energykin_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:57
% EndTime: 2019-03-08 18:05:57
% DurationCPUTime: 0.02s
% Computational Cost: add. (14->9), mult. (27->19), div. (0->0), fcn. (4->2), ass. (0->9)
t30 = qJD(1) ^ 2;
t32 = t30 / 0.2e1;
t31 = qJ(2) * qJD(1);
t29 = cos(qJ(3));
t28 = sin(qJ(3));
t27 = -qJD(1) + qJD(3);
t26 = -qJD(1) * pkin(1) + qJD(2);
t25 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t1 = [t32, 0, 0, -t26 * qJD(1), t30 * qJ(2), qJ(2) ^ 2 * t32 + t26 ^ 2 / 0.2e1, t27 ^ 2 / 0.2e1 (t29 * t25 - t28 * t31) * t27 -(t28 * t25 + t29 * t31) * t27;];
T_reg  = t1;
