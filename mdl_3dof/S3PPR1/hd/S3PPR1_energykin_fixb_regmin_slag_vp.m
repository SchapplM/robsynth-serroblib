% Calculate minimal parameter regressor of fixed base kinetic energy for
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% 
% Output:
% T_reg [1x5]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S3PPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_energykin_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR1_energykin_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_energykin_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:02:00
% EndTime: 2019-03-08 18:02:00
% DurationCPUTime: 0.01s
% Computational Cost: add. (4->4), mult. (14->11), div. (0->0), fcn. (4->2), ass. (0->4)
t12 = cos(qJ(3));
t11 = sin(qJ(3));
t10 = qJD(1) ^ 2 / 0.2e1;
t1 = [t10, t10 + qJD(2) ^ 2 / 0.2e1, qJD(3) ^ 2 / 0.2e1 (-t11 * qJD(1) + t12 * qJD(2)) * qJD(3) -(t12 * qJD(1) + t11 * qJD(2)) * qJD(3);];
T_reg  = t1;
