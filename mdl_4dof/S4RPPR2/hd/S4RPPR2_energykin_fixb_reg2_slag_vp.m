% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_energykin_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:31
% EndTime: 2018-11-14 13:47:31
% DurationCPUTime: 0.07s
% Computational Cost: add. (57->18), mult. (107->36), div. (0->0), fcn. (32->4), ass. (0->18)
t16 = qJD(1) ^ 2;
t11 = t16 / 0.2e1;
t17 = qJ(2) * qJD(1);
t15 = cos(qJ(4));
t14 = sin(qJ(4));
t13 = cos(pkin(6));
t12 = sin(pkin(6));
t10 = qJD(3) ^ 2 / 0.2e1;
t9 = -qJD(1) + qJD(4);
t8 = -qJD(1) * pkin(1) + qJD(2);
t7 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t6 = t13 * t7;
t5 = t12 * t7 + t13 * t17;
t4 = -t12 * t17 + t6;
t3 = t6 + (-qJ(2) * t12 - pkin(3)) * qJD(1);
t2 = t14 * t3 + t15 * t5;
t1 = -t14 * t5 + t15 * t3;
t18 = [0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, -t8 * qJD(1), 0, t16 * qJ(2), qJ(2) ^ 2 * t11 + t8 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t11, -t4 * qJD(1), t5 * qJD(1), 0, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1 + t10, 0, 0, 0, 0, 0, t9 ^ 2 / 0.2e1, t1 * t9, -t2 * t9, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t10;];
T_reg  = t18;
