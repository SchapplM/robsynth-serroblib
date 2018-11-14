% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta1]';
% 
% Output:
% T_reg [1x6]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:38
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:38:14
% EndTime: 2018-11-14 13:38:14
% DurationCPUTime: 0.02s
% Computational Cost: add. (6->5), mult. (20->13), div. (0->0), fcn. (4->2), ass. (0->5)
t28 = qJD(1) ^ 2 / 0.2e1;
t31 = t28 + qJD(2) ^ 2 / 0.2e1;
t30 = cos(qJ(4));
t29 = sin(qJ(4));
t1 = [t28, t31, qJD(3) ^ 2 / 0.2e1 + t31, qJD(4) ^ 2 / 0.2e1 (-t29 * qJD(2) + t30 * qJD(3)) * qJD(4) -(t30 * qJD(2) + t29 * qJD(3)) * qJD(4);];
T_reg  = t1;
