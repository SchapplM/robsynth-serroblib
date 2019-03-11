% Calculate inertial parameters regressor of fixed base kinetic energy for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% T_reg [1x(6*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_energykin_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_energykin_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_energykin_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:39
% EndTime: 2019-03-08 19:11:39
% DurationCPUTime: 0.22s
% Computational Cost: add. (759->55), mult. (1932->140), div. (0->0), fcn. (1689->16), ass. (0->55)
t43 = cos(pkin(6)) * qJD(1) + qJD(2);
t48 = sin(pkin(7));
t52 = cos(pkin(7));
t50 = cos(pkin(14));
t49 = sin(pkin(6));
t69 = qJD(1) * t49;
t66 = t50 * t69;
t77 = t43 * t48 + t52 * t66;
t56 = sin(qJ(3));
t58 = cos(qJ(3));
t46 = sin(pkin(14));
t67 = t46 * t69;
t24 = -t56 * t67 + t58 * t77;
t20 = qJD(3) * pkin(3) + t24;
t29 = t43 * t52 - t48 * t66;
t47 = sin(pkin(8));
t51 = cos(pkin(8));
t76 = t20 * t51 + t29 * t47;
t25 = t56 * t77 + t58 * t67;
t68 = qJD(3) * t47;
t19 = pkin(10) * t68 + t25;
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t11 = -t19 * t55 + t57 * t76;
t12 = t57 * t19 + t55 * t76;
t42 = qJD(3) * t51 + qJD(4);
t10 = pkin(11) * t42 + t12;
t27 = t51 * t29;
t14 = t27 + (-t20 + (-pkin(4) * t57 - pkin(11) * t55) * qJD(3)) * t47;
t54 = sin(qJ(5));
t75 = cos(qJ(5));
t6 = t10 * t75 + t14 * t54;
t74 = cos(qJ(6));
t59 = qJD(3) ^ 2;
t70 = t47 ^ 2 * t59;
t65 = t55 * t68;
t64 = t57 * t68;
t63 = t70 / 0.2e1;
t30 = -t42 * t75 + t54 * t65;
t5 = -t10 * t54 + t14 * t75;
t9 = -t42 * pkin(4) - t11;
t60 = qJD(1) ^ 2;
t53 = sin(qJ(6));
t40 = -qJD(5) + t64;
t32 = t42 * t54 + t65 * t75;
t28 = qJD(6) + t30;
t23 = t32 * t74 - t40 * t53;
t21 = t32 * t53 + t40 * t74;
t15 = -t47 * t20 + t27;
t7 = t30 * pkin(5) - t32 * pkin(12) + t9;
t4 = -pkin(12) * t40 + t6;
t3 = pkin(5) * t40 - t5;
t2 = t4 * t74 + t53 * t7;
t1 = -t4 * t53 + t7 * t74;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t60 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 ^ 2 / 0.2e1 + (t46 ^ 2 / 0.2e1 + t50 ^ 2 / 0.2e1) * t60 * t49 ^ 2, 0, 0, 0, 0, 0, t59 / 0.2e1, t24 * qJD(3), -t25 * qJD(3), 0, t25 ^ 2 / 0.2e1 + t24 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1, t55 ^ 2 * t63, t57 * t55 * t70, t42 * t65, t57 ^ 2 * t63, t42 * t64, t42 ^ 2 / 0.2e1, t11 * t42 - t15 * t64, -t12 * t42 + t15 * t65 (-t11 * t55 + t12 * t57) * t68, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t30, -t32 * t40, t30 ^ 2 / 0.2e1, t30 * t40, t40 ^ 2 / 0.2e1, t30 * t9 - t40 * t5, t32 * t9 + t40 * t6, -t30 * t6 - t32 * t5, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t21, t23 * t28, t21 ^ 2 / 0.2e1, -t21 * t28, t28 ^ 2 / 0.2e1, t1 * t28 + t21 * t3, -t2 * t28 + t23 * t3, -t1 * t23 - t2 * t21, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1;];
T_reg  = t8;
