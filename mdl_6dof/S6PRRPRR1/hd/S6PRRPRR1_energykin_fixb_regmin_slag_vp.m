% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:31
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:28:21
% EndTime: 2021-01-16 03:28:21
% DurationCPUTime: 0.12s
% Computational Cost: add. (338->49), mult. (828->105), div. (0->0), fcn. (638->12), ass. (0->46)
t63 = qJD(2) ^ 2;
t76 = t63 / 0.2e1;
t58 = sin(qJ(2));
t72 = qJD(1) * sin(pkin(6));
t45 = qJD(2) * pkin(8) + t58 * t72;
t61 = cos(qJ(3));
t71 = qJD(1) * cos(pkin(6));
t49 = t61 * t71;
t57 = sin(qJ(3));
t68 = qJ(4) * qJD(2);
t38 = qJD(3) * pkin(3) + t49 + (-t45 - t68) * t57;
t74 = t61 * t45 + t57 * t71;
t39 = t61 * t68 + t74;
t52 = sin(pkin(12));
t73 = cos(pkin(12));
t27 = t73 * t38 - t52 * t39;
t43 = (t52 * t61 + t73 * t57) * qJD(2);
t25 = qJD(3) * pkin(4) - t43 * pkin(9) + t27;
t28 = t52 * t38 + t73 * t39;
t69 = qJD(2) * t61;
t70 = qJD(2) * t57;
t42 = t52 * t70 - t73 * t69;
t26 = -t42 * pkin(9) + t28;
t56 = sin(qJ(5));
t60 = cos(qJ(5));
t75 = t56 * t25 + t60 * t26;
t67 = qJD(2) * qJD(3);
t62 = cos(qJ(2));
t66 = t62 * t72;
t65 = qJD(2) * t72;
t32 = t60 * t42 + t56 * t43;
t64 = t60 * t25 - t56 * t26;
t41 = -t66 + qJD(4) + (-pkin(3) * t61 - pkin(2)) * qJD(2);
t34 = t42 * pkin(4) + t41;
t59 = cos(qJ(6));
t55 = sin(qJ(6));
t51 = qJD(3) + qJD(5);
t46 = -qJD(2) * pkin(2) - t66;
t33 = -t56 * t42 + t60 * t43;
t31 = qJD(6) + t32;
t30 = t59 * t33 + t55 * t51;
t29 = t55 * t33 - t59 * t51;
t22 = t32 * pkin(5) - t33 * pkin(10) + t34;
t21 = t51 * pkin(10) + t75;
t20 = -t51 * pkin(5) - t64;
t1 = [qJD(1) ^ 2 / 0.2e1, t76, t62 * t65, -t58 * t65, t57 ^ 2 * t76, t57 * t63 * t61, t57 * t67, t61 * t67, qJD(3) ^ 2 / 0.2e1, (-t57 * t45 + t49) * qJD(3) - t46 * t69, -t74 * qJD(3) + t46 * t70, t27 * qJD(3) + t41 * t42, -t28 * qJD(3) + t41 * t43, -t27 * t43 - t28 * t42, t28 ^ 2 / 0.2e1 + t27 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1, t33 ^ 2 / 0.2e1, -t33 * t32, t33 * t51, -t32 * t51, t51 ^ 2 / 0.2e1, t34 * t32 + t64 * t51, t34 * t33 - t75 * t51, t30 ^ 2 / 0.2e1, -t30 * t29, t30 * t31, -t29 * t31, t31 ^ 2 / 0.2e1, (-t55 * t21 + t59 * t22) * t31 + t20 * t29, -(t59 * t21 + t55 * t22) * t31 + t20 * t30;];
T_reg = t1;
