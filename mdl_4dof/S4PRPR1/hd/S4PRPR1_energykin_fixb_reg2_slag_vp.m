% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:22
% EndTime: 2018-11-14 13:42:22
% DurationCPUTime: 0.06s
% Computational Cost: add. (23->12), mult. (51->24), div. (0->0), fcn. (8->2), ass. (0->12)
t10 = qJD(2) ^ 2;
t6 = t10 / 0.2e1;
t11 = qJ(3) * qJD(2);
t9 = cos(qJ(4));
t8 = sin(qJ(4));
t7 = qJD(1) ^ 2 / 0.2e1;
t5 = -qJD(2) + qJD(4);
t4 = -qJD(2) * pkin(2) + qJD(3);
t3 = qJD(3) + (-pkin(2) - pkin(3)) * qJD(2);
t2 = t9 * t11 + t8 * t3;
t1 = -t8 * t11 + t9 * t3;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, t6, 0, 0, 0, t7, 0, 0, 0, t6, 0, 0, -t4 * qJD(2), 0, t10 * qJ(3), qJ(3) ^ 2 * t6 + t7 + t4 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t5 ^ 2 / 0.2e1, t1 * t5, -t2 * t5, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7;];
T_reg  = t12;
