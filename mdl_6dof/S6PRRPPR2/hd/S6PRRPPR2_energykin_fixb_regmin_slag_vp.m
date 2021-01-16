% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:16
% EndTime: 2021-01-16 02:21:16
% DurationCPUTime: 0.11s
% Computational Cost: add. (256->47), mult. (616->98), div. (0->0), fcn. (435->10), ass. (0->42)
t54 = qJD(2) ^ 2;
t67 = t54 / 0.2e1;
t66 = pkin(4) + pkin(9);
t50 = sin(qJ(2));
t64 = qJD(1) * sin(pkin(6));
t39 = qJD(2) * pkin(8) + t50 * t64;
t52 = cos(qJ(3));
t63 = qJD(1) * cos(pkin(6));
t43 = t52 * t63;
t49 = sin(qJ(3));
t60 = qJ(4) * qJD(2);
t29 = qJD(3) * pkin(3) + t43 + (-t39 - t60) * t49;
t65 = t52 * t39 + t49 * t63;
t30 = t52 * t60 + t65;
t44 = sin(pkin(11));
t46 = cos(pkin(11));
t24 = t44 * t29 + t46 * t30;
t62 = qJD(2) * t49;
t61 = qJD(2) * t52;
t59 = qJD(2) * qJD(3);
t53 = cos(qJ(2));
t58 = t53 * t64;
t57 = qJD(2) * t64;
t23 = t46 * t29 - t44 * t30;
t56 = qJD(5) - t23;
t21 = -qJD(3) * qJ(5) - t24;
t37 = (t44 * t52 + t46 * t49) * qJD(2);
t34 = -t58 + qJD(4) + (-pkin(3) * t52 - pkin(2)) * qJD(2);
t55 = -t37 * qJ(5) + t34;
t51 = cos(qJ(6));
t48 = sin(qJ(6));
t40 = -qJD(2) * pkin(2) - t58;
t36 = t44 * t62 - t46 * t61;
t35 = qJD(6) + t37;
t32 = t51 * qJD(3) + t48 * t36;
t31 = t48 * qJD(3) - t51 * t36;
t25 = t36 * pkin(4) + t55;
t22 = t66 * t36 + t55;
t20 = -qJD(3) * pkin(4) + t56;
t19 = -t36 * pkin(5) - t21;
t18 = t37 * pkin(5) - t66 * qJD(3) + t56;
t1 = [qJD(1) ^ 2 / 0.2e1, t67, t53 * t57, -t50 * t57, t49 ^ 2 * t67, t49 * t54 * t52, t49 * t59, t52 * t59, qJD(3) ^ 2 / 0.2e1, (-t49 * t39 + t43) * qJD(3) - t40 * t61, -t65 * qJD(3) + t40 * t62, t23 * qJD(3) + t34 * t36, -t24 * qJD(3) + t34 * t37, -t23 * t37 - t24 * t36, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t34 ^ 2 / 0.2e1, t20 * t37 + t21 * t36, t20 * qJD(3) - t25 * t36, -t21 * qJD(3) - t25 * t37, t25 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + t20 ^ 2 / 0.2e1, t32 ^ 2 / 0.2e1, -t32 * t31, t32 * t35, -t31 * t35, t35 ^ 2 / 0.2e1, (t51 * t18 - t48 * t22) * t35 + t19 * t31, -(t48 * t18 + t51 * t22) * t35 + t19 * t32;];
T_reg = t1;
