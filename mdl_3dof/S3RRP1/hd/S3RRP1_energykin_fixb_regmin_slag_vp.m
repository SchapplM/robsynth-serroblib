% Calculate minimal parameter regressor of fixed base kinetic energy for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% 
% Output:
% T_reg [1x9]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S3RRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_energykin_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_energykin_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_energykin_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:08
% EndTime: 2018-11-14 10:15:08
% DurationCPUTime: 0.02s
% Computational Cost: add. (18->7), mult. (28->17), div. (0->0), fcn. (6->2), ass. (0->7)
t29 = pkin(1) * qJD(1);
t28 = sin(qJ(2)) * t29;
t27 = cos(qJ(2)) * t29;
t24 = qJD(1) + qJD(2);
t23 = t24 * qJ(3) + t28;
t22 = -t24 * pkin(2) + qJD(3) - t27;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t24 ^ 2 / 0.2e1, t24 * t27, -t24 * t28, -t22 * t24, t23 * t24, t23 ^ 2 / 0.2e1 + t22 ^ 2 / 0.2e1;];
T_reg  = t1;
