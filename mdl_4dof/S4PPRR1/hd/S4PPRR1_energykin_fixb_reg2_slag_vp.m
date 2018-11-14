% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PPRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:18
% EndTime: 2018-11-14 13:40:18
% DurationCPUTime: 0.06s
% Computational Cost: add. (18->11), mult. (52->28), div. (0->0), fcn. (20->4), ass. (0->13)
t7 = sin(qJ(3));
t12 = qJD(2) * t7;
t11 = qJD(2) * qJD(3);
t10 = qJD(2) ^ 2;
t9 = cos(qJ(3));
t8 = cos(qJ(4));
t6 = sin(qJ(4));
t5 = qJD(1) ^ 2 / 0.2e1;
t4 = qJD(3) + qJD(4);
t3 = qJD(3) * pkin(3) + t9 * qJD(2);
t2 = t8 * t12 + t6 * t3;
t1 = -t6 * t12 + t8 * t3;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 + t10 / 0.2e1, 0, 0, 0, 0, 0, qJD(3) ^ 2 / 0.2e1, t9 * t11, -t7 * t11, 0, t5 + (t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1) * t10, 0, 0, 0, 0, 0, t4 ^ 2 / 0.2e1, t1 * t4, -t2 * t4, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5;];
T_reg  = t13;
