% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% T_reg [1x(4*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:34
% EndTime: 2019-07-18 18:16:34
% DurationCPUTime: 0.07s
% Computational Cost: add. (60->15), mult. (86->36), div. (0->0), fcn. (24->4), ass. (0->18)
t18 = pkin(1) * qJD(1);
t8 = qJD(1) + qJD(2);
t10 = sin(qJ(2));
t17 = t10 * t18;
t12 = cos(qJ(2));
t16 = t12 * t18;
t15 = qJD(3) - t16;
t13 = qJD(1) ^ 2;
t11 = cos(qJ(4));
t9 = sin(qJ(4));
t7 = qJD(4) - t8;
t6 = t8 ^ 2 / 0.2e1;
t5 = t8 * qJ(3) + t17;
t4 = -t8 * pkin(2) + t15;
t3 = (-pkin(2) - pkin(3)) * t8 + t15;
t2 = t11 * t5 + t9 * t3;
t1 = t11 * t3 - t9 * t5;
t14 = [0, 0, 0, 0, 0, t13 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t8 * t16, -t8 * t17, 0, (t10 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t13, 0, 0, 0, t6, 0, 0, -t4 * t8, 0, t5 * t8, t5 ^ 2 / 0.2e1 + t4 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t7 ^ 2 / 0.2e1, t1 * t7, -t2 * t7, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t14;
