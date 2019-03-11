% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP10_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energykin_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:31
% EndTime: 2019-03-09 03:32:31
% DurationCPUTime: 0.15s
% Computational Cost: add. (267->52), mult. (520->100), div. (0->0), fcn. (226->4), ass. (0->44)
t39 = qJD(1) ^ 2;
t31 = t39 / 0.2e1;
t36 = sin(qJ(5));
t52 = cos(qJ(5));
t22 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t38 = cos(qJ(3));
t7 = qJD(4) + (pkin(4) * qJD(1) - t22) * t38 + (-pkin(3) - pkin(8)) * qJD(3);
t37 = sin(qJ(3));
t47 = qJD(1) * t37;
t49 = pkin(3) * t47 + qJD(1) * qJ(2);
t9 = (pkin(8) * t37 - qJ(4) * t38) * qJD(1) + t49;
t4 = t36 * t7 + t52 * t9;
t16 = t36 * qJD(3) - t52 * t47;
t18 = t52 * qJD(3) + t36 * t47;
t51 = t18 * t16;
t45 = t38 * qJD(1);
t24 = qJD(5) + t45;
t50 = t24 * t16;
t14 = -qJD(3) * qJ(4) - t37 * t22;
t48 = t39 * qJ(2);
t46 = qJD(3) * t22;
t44 = t16 ^ 2 / 0.2e1;
t43 = qJD(1) * qJD(3);
t42 = t38 * t39 * t37;
t41 = t37 * t43;
t40 = t38 * t43;
t10 = -pkin(4) * t47 - t14;
t3 = -t36 * t9 + t52 * t7;
t35 = t38 ^ 2;
t34 = t37 ^ 2;
t30 = qJD(3) ^ 2 / 0.2e1;
t29 = qJ(2) ^ 2 * t31;
t28 = -pkin(1) * qJD(1) + qJD(2);
t26 = t35 * t31;
t25 = t34 * t31;
t20 = t24 ^ 2 / 0.2e1;
t15 = t18 ^ 2 / 0.2e1;
t13 = -qJ(4) * t45 + t49;
t12 = -qJD(3) * pkin(3) - t38 * t22 + qJD(4);
t11 = t18 * t24;
t5 = t16 * pkin(5) - t18 * qJ(6) + t10;
t2 = t24 * qJ(6) + t4;
t1 = -t24 * pkin(5) + qJD(6) - t3;
t6 = [0, 0, 0, 0, 0, t31, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, t28 * qJD(1), t48, t29 + t28 ^ 2 / 0.2e1, t26, -t42, t40, t25, -t41, t30, t37 * t48 + t38 * t46, -t37 * t46 + t38 * t48 (-t34 - t35) * t22 * qJD(1), t29 + (t34 / 0.2e1 + t35 / 0.2e1) * t22 ^ 2, t30, -t40, t41, t26, -t42, t25 (t12 * t38 + t14 * t37) * qJD(1), t12 * qJD(3) - t13 * t47, -t14 * qJD(3) - t13 * t45, t13 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1 + t12 ^ 2 / 0.2e1, t15, -t51, t11, t44, -t50, t20, t10 * t16 + t3 * t24, t10 * t18 - t4 * t24, -t4 * t16 - t3 * t18, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1, t15, t11, t51, t20, t50, t44, -t1 * t24 + t5 * t16, t1 * t18 - t2 * t16, -t5 * t18 + t2 * t24, t2 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
