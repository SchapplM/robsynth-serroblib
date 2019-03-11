% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_energykin_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:23
% EndTime: 2019-03-08 19:05:23
% DurationCPUTime: 0.20s
% Computational Cost: add. (510->52), mult. (1254->132), div. (0->0), fcn. (1035->14), ass. (0->48)
t37 = cos(pkin(6)) * qJD(1) + qJD(2);
t43 = sin(pkin(7));
t46 = cos(pkin(7));
t45 = cos(pkin(13));
t44 = sin(pkin(6));
t63 = qJD(1) * t44;
t58 = t45 * t63;
t68 = t37 * t43 + t46 * t58;
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t42 = sin(pkin(13));
t59 = t42 * t63;
t20 = -t50 * t59 + t68 * t52;
t53 = qJD(3) ^ 2;
t67 = t53 / 0.2e1;
t21 = t68 * t50 + t52 * t59;
t19 = qJD(3) * pkin(9) + t21;
t26 = t46 * t37 - t43 * t58;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t12 = t51 * t19 + t49 * t26;
t10 = qJD(4) * pkin(10) + t12;
t15 = (-pkin(4) * t51 - pkin(10) * t49 - pkin(3)) * qJD(3) - t20;
t48 = sin(qJ(5));
t66 = cos(qJ(5));
t6 = t66 * t10 + t48 * t15;
t65 = cos(qJ(6));
t62 = qJD(3) * t49;
t61 = t51 * qJD(3);
t60 = qJD(3) * qJD(4);
t5 = -t48 * t10 + t66 * t15;
t11 = -t49 * t19 + t51 * t26;
t38 = -qJD(5) + t61;
t9 = -qJD(4) * pkin(4) - t11;
t54 = qJD(1) ^ 2;
t47 = sin(qJ(6));
t35 = -qJD(6) + t38;
t32 = t48 * qJD(4) + t66 * t62;
t30 = -t66 * qJD(4) + t48 * t62;
t24 = -t47 * t30 + t65 * t32;
t22 = t65 * t30 + t47 * t32;
t18 = -qJD(3) * pkin(3) - t20;
t7 = t30 * pkin(5) + t9;
t4 = -t30 * pkin(11) + t6;
t3 = -t38 * pkin(5) - t32 * pkin(11) + t5;
t2 = t47 * t3 + t65 * t4;
t1 = t65 * t3 - t47 * t4;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t54 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 ^ 2 / 0.2e1 + (t42 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1) * t54 * t44 ^ 2, 0, 0, 0, 0, 0, t67, t20 * qJD(3), -t21 * qJD(3), 0, t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t49 ^ 2 * t67, t51 * t53 * t49, t49 * t60, t51 ^ 2 * t67, t51 * t60, qJD(4) ^ 2 / 0.2e1, t11 * qJD(4) - t18 * t61, -t12 * qJD(4) + t18 * t62 (-t11 * t49 + t12 * t51) * qJD(3), t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t30, -t32 * t38, t30 ^ 2 / 0.2e1, t30 * t38, t38 ^ 2 / 0.2e1, t9 * t30 - t5 * t38, t9 * t32 + t6 * t38, -t6 * t30 - t5 * t32, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, -t24 * t35, t22 ^ 2 / 0.2e1, t22 * t35, t35 ^ 2 / 0.2e1, -t1 * t35 + t7 * t22, t2 * t35 + t7 * t24, -t1 * t24 - t2 * t22, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg  = t8;
