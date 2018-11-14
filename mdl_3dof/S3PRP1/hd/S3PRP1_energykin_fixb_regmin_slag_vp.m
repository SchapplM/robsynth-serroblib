% Calculate minimal parameter regressor of fixed base kinetic energy for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% 
% Output:
% T_reg [1x7]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:04
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S3PRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_energykin_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_energykin_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_energykin_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:04:15
% EndTime: 2018-11-14 10:04:15
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->6), mult. (22->17), div. (0->0), fcn. (6->2), ass. (0->6)
t18 = qJD(1) * qJD(2);
t17 = cos(qJ(2));
t16 = sin(qJ(2));
t15 = qJD(2) * qJ(3) + t16 * qJD(1);
t14 = -qJD(2) * pkin(2) - t17 * qJD(1) + qJD(3);
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t17 * t18, -t16 * t18, -t14 * qJD(2), t15 * qJD(2), t15 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1;];
T_reg  = t1;
