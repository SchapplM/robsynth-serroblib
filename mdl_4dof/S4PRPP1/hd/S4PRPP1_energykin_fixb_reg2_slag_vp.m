% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:42
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRPP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:41:09
% EndTime: 2018-11-14 13:41:09
% DurationCPUTime: 0.05s
% Computational Cost: add. (13->9), mult. (37->18), div. (0->0), fcn. (0->0), ass. (0->7)
t6 = qJD(2) ^ 2;
t4 = t6 / 0.2e1;
t5 = qJD(1) ^ 2 / 0.2e1;
t3 = -qJD(2) * pkin(2) + qJD(3);
t2 = qJD(2) * qJ(3) + qJD(4);
t1 = qJD(3) + (-pkin(2) - qJ(4)) * qJD(2);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, t4, 0, 0, 0, t5, t4, 0, 0, 0, 0, 0, 0, t3 * qJD(2), t6 * qJ(3), t5 + qJ(3) ^ 2 * t4 + t3 ^ 2 / 0.2e1, t4, 0, 0, 0, 0, 0, 0, t2 * qJD(2), -t1 * qJD(2), t5 + t1 ^ 2 / 0.2e1 + t2 ^ 2 / 0.2e1;];
T_reg  = t7;
