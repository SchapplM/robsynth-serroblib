% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% 
% Output:
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRPP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_energykin_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:01:24
% EndTime: 2018-11-14 14:01:24
% DurationCPUTime: 0.02s
% Computational Cost: add. (20->11), mult. (40->22), div. (0->0), fcn. (10->2), ass. (0->10)
t32 = qJD(1) * qJD(2);
t30 = cos(qJ(2));
t31 = -t30 * qJD(1) + qJD(3);
t29 = sin(qJ(2));
t28 = qJD(2) * qJ(3) + t29 * qJD(1);
t27 = -qJD(2) * pkin(2) + t31;
t26 = t28 * qJD(2);
t25 = t28 ^ 2 / 0.2e1;
t24 = (-pkin(2) - pkin(3)) * qJD(2) + t31;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t30 * t32, -t29 * t32, -t27 * qJD(2), t26, t25 + t27 ^ 2 / 0.2e1, -t24 * qJD(2), t26, t25 + t24 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t1;
