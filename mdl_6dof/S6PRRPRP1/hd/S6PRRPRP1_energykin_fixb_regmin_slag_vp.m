% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:38:47
% EndTime: 2021-01-16 02:38:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (338->46), mult. (795->98), div. (0->0), fcn. (580->10), ass. (0->41)
t61 = qJD(2) ^ 2;
t74 = t61 / 0.2e1;
t73 = cos(qJ(5));
t58 = sin(qJ(2));
t70 = qJD(1) * sin(pkin(6));
t46 = qJD(2) * pkin(8) + t58 * t70;
t59 = cos(qJ(3));
t69 = qJD(1) * cos(pkin(6));
t50 = t59 * t69;
t57 = sin(qJ(3));
t66 = qJ(4) * qJD(2);
t37 = qJD(3) * pkin(3) + t50 + (-t46 - t66) * t57;
t71 = t59 * t46 + t57 * t69;
t38 = t59 * t66 + t71;
t52 = sin(pkin(11));
t54 = cos(pkin(11));
t30 = t52 * t37 + t54 * t38;
t28 = qJD(3) * pkin(9) + t30;
t60 = cos(qJ(2));
t64 = t60 * t70;
t41 = -t64 + qJD(4) + (-pkin(3) * t59 - pkin(2)) * qJD(2);
t67 = qJD(2) * t59;
t68 = qJD(2) * t57;
t43 = t52 * t68 - t54 * t67;
t44 = (t52 * t59 + t54 * t57) * qJD(2);
t33 = t43 * pkin(4) - t44 * pkin(9) + t41;
t56 = sin(qJ(5));
t72 = t73 * t28 + t56 * t33;
t65 = qJD(2) * qJD(3);
t63 = qJD(2) * t70;
t62 = -t56 * t28 + t73 * t33;
t29 = t54 * t37 - t52 * t38;
t27 = -qJD(3) * pkin(4) - t29;
t47 = -qJD(2) * pkin(2) - t64;
t42 = qJD(5) + t43;
t40 = t56 * qJD(3) + t73 * t44;
t39 = -t73 * qJD(3) + t56 * t44;
t25 = t39 * pkin(5) + qJD(6) + t27;
t24 = -t39 * qJ(6) + t72;
t23 = t42 * pkin(5) - t40 * qJ(6) + t62;
t1 = [qJD(1) ^ 2 / 0.2e1, t74, t60 * t63, -t58 * t63, t57 ^ 2 * t74, t57 * t61 * t59, t57 * t65, t59 * t65, qJD(3) ^ 2 / 0.2e1, (-t57 * t46 + t50) * qJD(3) - t47 * t67, -t71 * qJD(3) + t47 * t68, t29 * qJD(3) + t41 * t43, -t30 * qJD(3) + t41 * t44, -t29 * t44 - t30 * t43, t30 ^ 2 / 0.2e1 + t29 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1, t40 ^ 2 / 0.2e1, -t40 * t39, t40 * t42, -t39 * t42, t42 ^ 2 / 0.2e1, t27 * t39 + t62 * t42, t27 * t40 - t72 * t42, t23 * t42 + t25 * t39, -t24 * t42 + t25 * t40, -t23 * t40 - t24 * t39, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1;];
T_reg = t1;
