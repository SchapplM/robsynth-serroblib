% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:29
% EndTime: 2018-11-14 14:12:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->12), mult. (76->31), div. (0->0), fcn. (36->4), ass. (0->16)
t11 = cos(qJ(3));
t10 = sin(qJ(2));
t15 = qJD(1) * t10;
t12 = cos(qJ(2));
t6 = qJD(2) * pkin(2) + t12 * qJD(1);
t9 = sin(qJ(3));
t4 = t11 * t15 + t9 * t6;
t8 = qJD(2) + qJD(3);
t16 = t4 * t8;
t14 = qJD(1) * qJD(2);
t3 = t11 * t6 - t9 * t15;
t13 = qJD(1) ^ 2;
t7 = t8 ^ 2 / 0.2e1;
t2 = t4 ^ 2 / 0.2e1;
t1 = t8 * pkin(3) + t3;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t13 / 0.2e1, 0, 0, 0, 0, 0, qJD(2) ^ 2 / 0.2e1, t12 * t14, -t10 * t14, 0 (t10 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1) * t13, 0, 0, 0, 0, 0, t7, t3 * t8, -t16, 0, t2 + t3 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t7, t1 * t8, -t16, 0, t2 + t1 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t5;
