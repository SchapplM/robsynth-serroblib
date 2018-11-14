% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_energykin_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:44:29
% EndTime: 2018-11-14 13:44:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (27->11), mult. (66->30), div. (0->0), fcn. (20->4), ass. (0->15)
t15 = pkin(2) * qJD(2);
t5 = qJD(2) + qJD(3);
t8 = sin(qJ(3));
t14 = t8 * t15;
t10 = cos(qJ(3));
t13 = t10 * t15;
t11 = qJD(2) ^ 2;
t9 = cos(qJ(4));
t7 = sin(qJ(4));
t6 = qJD(1) ^ 2 / 0.2e1;
t4 = qJD(4) + t5;
t3 = t5 * pkin(3) + t13;
t2 = t14 * t9 + t7 * t3;
t1 = -t14 * t7 + t9 * t3;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, t11 / 0.2e1, 0, 0, 0, t6, 0, 0, 0, 0, 0, t5 ^ 2 / 0.2e1, t5 * t13, -t5 * t14, 0, t6 + (t8 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1) * pkin(2) ^ 2 * t11, 0, 0, 0, 0, 0, t4 ^ 2 / 0.2e1, t1 * t4, -t2 * t4, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t6;];
T_reg  = t12;
