% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR3_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:11:25
% EndTime: 2018-11-14 14:11:25
% DurationCPUTime: 0.07s
% Computational Cost: add. (43->15), mult. (108->38), div. (0->0), fcn. (60->6), ass. (0->19)
t14 = sin(qJ(2));
t19 = qJD(1) * t14;
t18 = qJD(1) * qJD(2);
t11 = sin(pkin(6));
t12 = cos(pkin(6));
t16 = cos(qJ(2));
t7 = qJD(2) * pkin(2) + t16 * qJD(1);
t4 = -t11 * t19 + t12 * t7;
t17 = qJD(1) ^ 2;
t15 = cos(qJ(4));
t13 = sin(qJ(4));
t10 = qJD(2) ^ 2 / 0.2e1;
t9 = qJD(3) ^ 2 / 0.2e1;
t8 = qJD(2) + qJD(4);
t5 = t11 * t7 + t12 * t19;
t3 = qJD(2) * pkin(3) + t4;
t2 = t13 * t3 + t15 * t5;
t1 = -t13 * t5 + t15 * t3;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t17 / 0.2e1, 0, 0, 0, 0, 0, t10, t16 * t18, -t14 * t18, 0 (t14 ^ 2 / 0.2e1 + t16 ^ 2 / 0.2e1) * t17, 0, 0, 0, 0, 0, t10, t4 * qJD(2), -t5 * qJD(2), 0, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t9, 0, 0, 0, 0, 0, t8 ^ 2 / 0.2e1, t1 * t8, -t2 * t8, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t9;];
T_reg  = t6;
