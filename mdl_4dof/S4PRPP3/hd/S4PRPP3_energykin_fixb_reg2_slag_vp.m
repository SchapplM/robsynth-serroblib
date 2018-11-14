% Calculate inertial parameters regressor of fixed base kinetic energy for
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
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRPP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_energykin_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:01:24
% EndTime: 2018-11-14 14:01:24
% DurationCPUTime: 0.06s
% Computational Cost: add. (21->12), mult. (52->26), div. (0->0), fcn. (12->2), ass. (0->12)
t11 = qJD(1) * qJD(2);
t8 = cos(qJ(2));
t10 = -t8 * qJD(1) + qJD(3);
t9 = qJD(1) ^ 2;
t7 = sin(qJ(2));
t6 = qJD(2) ^ 2 / 0.2e1;
t5 = qJD(2) * qJ(3) + t7 * qJD(1);
t4 = -qJD(2) * pkin(2) + t10;
t3 = t5 * qJD(2);
t2 = t5 ^ 2 / 0.2e1;
t1 = (-pkin(2) - pkin(3)) * qJD(2) + t10;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t9 / 0.2e1, 0, 0, 0, 0, 0, t6, t8 * t11, -t7 * t11, 0 (t7 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1) * t9, 0, 0, 0, t6, 0, 0, -t4 * qJD(2), 0, t3, t2 + t4 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t6, -t1 * qJD(2), t3, 0, t2 + t1 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t12;
