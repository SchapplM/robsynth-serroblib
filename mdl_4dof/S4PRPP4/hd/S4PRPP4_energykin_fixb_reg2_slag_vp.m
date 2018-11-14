% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRPP4_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP4_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:29
% EndTime: 2018-11-14 14:09:29
% DurationCPUTime: 0.07s
% Computational Cost: add. (30->14), mult. (80->33), div. (0->0), fcn. (36->4), ass. (0->15)
t10 = sin(pkin(5));
t11 = cos(pkin(5));
t12 = sin(qJ(2));
t16 = qJD(1) * t12;
t13 = cos(qJ(2));
t6 = qJD(2) * pkin(2) + t13 * qJD(1);
t4 = t10 * t6 + t11 * t16;
t15 = qJD(1) * qJD(2);
t3 = -t10 * t16 + t11 * t6;
t14 = qJD(1) ^ 2;
t9 = qJD(2) ^ 2 / 0.2e1;
t8 = qJD(3) ^ 2 / 0.2e1;
t2 = qJD(2) * qJ(4) + t4;
t1 = -qJD(2) * pkin(3) + qJD(4) - t3;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t14 / 0.2e1, 0, 0, 0, 0, 0, t9, t13 * t15, -t12 * t15, 0 (t12 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1) * t14, 0, 0, 0, 0, 0, t9, t3 * qJD(2), -t4 * qJD(2), 0, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t8, 0, 0, 0, t9, 0, 0, -t1 * qJD(2), 0, t2 * qJD(2), t2 ^ 2 / 0.2e1 + t8 + t1 ^ 2 / 0.2e1;];
T_reg  = t5;
