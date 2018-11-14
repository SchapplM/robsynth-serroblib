% Calculate minimal parameter regressor of fixed base kinetic energy for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% 
% Output:
% T_reg [1x7]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S3PRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_energykin_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_energykin_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_energykin_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:57
% EndTime: 2018-11-14 10:11:57
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->6), mult. (22->18), div. (0->0), fcn. (10->4), ass. (0->9)
t21 = sin(qJ(2));
t25 = qJD(1) * t21;
t24 = qJD(1) * qJD(2);
t23 = cos(qJ(2));
t22 = cos(qJ(3));
t20 = sin(qJ(3));
t19 = qJD(2) + qJD(3);
t18 = qJD(2) * pkin(2) + t23 * qJD(1);
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t23 * t24, -t21 * t24, t19 ^ 2 / 0.2e1 (t22 * t18 - t20 * t25) * t19 -(t20 * t18 + t22 * t25) * t19;];
T_reg  = t1;
