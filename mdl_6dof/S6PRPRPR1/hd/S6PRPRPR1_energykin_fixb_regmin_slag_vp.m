% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:06:17
% EndTime: 2021-01-16 01:06:17
% DurationCPUTime: 0.10s
% Computational Cost: add. (240->42), mult. (557->93), div. (0->0), fcn. (407->12), ass. (0->41)
t67 = qJD(2) ^ 2;
t76 = t67 / 0.2e1;
t66 = cos(qJ(2));
t74 = qJD(1) * sin(pkin(6));
t48 = qJD(2) * pkin(2) + t66 * t74;
t57 = sin(pkin(11));
t60 = cos(pkin(11));
t63 = sin(qJ(2));
t69 = t63 * t74;
t41 = t57 * t48 + t60 * t69;
t39 = qJD(2) * pkin(8) + t41;
t54 = cos(pkin(6)) * qJD(1) + qJD(3);
t65 = cos(qJ(4));
t52 = t65 * t54;
t62 = sin(qJ(4));
t71 = qJ(5) * qJD(2);
t34 = qJD(4) * pkin(4) + t52 + (-t39 - t71) * t62;
t75 = t65 * t39 + t62 * t54;
t35 = t65 * t71 + t75;
t56 = sin(pkin(12));
t59 = cos(pkin(12));
t30 = t56 * t34 + t59 * t35;
t73 = qJD(2) * t62;
t72 = qJD(2) * t65;
t70 = qJD(2) * qJD(4);
t68 = qJD(2) * t74;
t40 = t60 * t48 - t57 * t69;
t45 = t56 * t73 - t59 * t72;
t29 = t59 * t34 - t56 * t35;
t36 = qJD(5) + (-pkin(4) * t65 - pkin(3)) * qJD(2) - t40;
t64 = cos(qJ(6));
t61 = sin(qJ(6));
t46 = (t56 * t65 + t59 * t62) * qJD(2);
t44 = qJD(6) + t45;
t43 = t61 * qJD(4) + t64 * t46;
t42 = -t64 * qJD(4) + t61 * t46;
t38 = -qJD(2) * pkin(3) - t40;
t31 = t45 * pkin(5) - t46 * pkin(9) + t36;
t28 = qJD(4) * pkin(9) + t30;
t27 = -qJD(4) * pkin(5) - t29;
t1 = [qJD(1) ^ 2 / 0.2e1, t76, t66 * t68, -t63 * t68, t41 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1 + t54 ^ 2 / 0.2e1, t62 ^ 2 * t76, t62 * t67 * t65, t62 * t70, t65 * t70, qJD(4) ^ 2 / 0.2e1, (-t62 * t39 + t52) * qJD(4) - t38 * t72, -t75 * qJD(4) + t38 * t73, t29 * qJD(4) + t36 * t45, -t30 * qJD(4) + t36 * t46, -t29 * t46 - t30 * t45, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t36 ^ 2 / 0.2e1, t43 ^ 2 / 0.2e1, -t43 * t42, t43 * t44, -t42 * t44, t44 ^ 2 / 0.2e1, (-t61 * t28 + t64 * t31) * t44 + t27 * t42, -(t64 * t28 + t61 * t31) * t44 + t27 * t43;];
T_reg = t1;
