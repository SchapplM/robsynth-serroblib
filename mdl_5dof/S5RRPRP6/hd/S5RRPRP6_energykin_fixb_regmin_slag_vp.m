% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:29:56
% EndTime: 2021-01-15 20:29:56
% DurationCPUTime: 0.08s
% Computational Cost: add. (244->38), mult. (601->83), div. (0->0), fcn. (390->6), ass. (0->34)
t52 = qJD(1) ^ 2;
t63 = t52 / 0.2e1;
t62 = cos(qJ(4));
t51 = cos(qJ(2));
t61 = t51 * t52;
t60 = pkin(6) + qJ(3);
t47 = sin(pkin(8));
t48 = cos(pkin(8));
t57 = qJD(1) * t51;
t50 = sin(qJ(2));
t58 = qJD(1) * t50;
t38 = t47 * t58 - t48 * t57;
t39 = (t47 * t51 + t48 * t50) * qJD(1);
t44 = qJD(3) + (-pkin(2) * t51 - pkin(1)) * qJD(1);
t28 = t38 * pkin(3) - t39 * pkin(7) + t44;
t42 = qJD(2) * pkin(2) - t60 * t58;
t43 = t60 * t57;
t33 = t47 * t42 + t48 * t43;
t31 = qJD(2) * pkin(7) + t33;
t49 = sin(qJ(4));
t59 = t49 * t28 + t62 * t31;
t56 = qJD(1) * qJD(2);
t55 = t50 * t56;
t54 = t51 * t56;
t53 = t62 * t28 - t49 * t31;
t32 = t48 * t42 - t47 * t43;
t30 = -qJD(2) * pkin(3) - t32;
t37 = qJD(4) + t38;
t35 = t49 * qJD(2) + t62 * t39;
t34 = -t62 * qJD(2) + t49 * t39;
t25 = t34 * pkin(4) + qJD(5) + t30;
t24 = -t34 * qJ(5) + t59;
t23 = t37 * pkin(4) - t35 * qJ(5) + t53;
t1 = [t63, 0, 0, t50 ^ 2 * t63, t50 * t61, t55, t54, qJD(2) ^ 2 / 0.2e1, pkin(1) * t61 - pkin(6) * t55, -t52 * pkin(1) * t50 - pkin(6) * t54, t32 * qJD(2) + t44 * t38, -t33 * qJD(2) + t44 * t39, -t32 * t39 - t33 * t38, t33 ^ 2 / 0.2e1 + t32 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1, t35 ^ 2 / 0.2e1, -t35 * t34, t35 * t37, -t34 * t37, t37 ^ 2 / 0.2e1, t30 * t34 + t53 * t37, t30 * t35 - t59 * t37, t23 * t37 + t25 * t34, -t24 * t37 + t25 * t35, -t23 * t35 - t24 * t34, t24 ^ 2 / 0.2e1 + t23 ^ 2 / 0.2e1 + t25 ^ 2 / 0.2e1;];
T_reg = t1;
