% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:10
% EndTime: 2021-01-16 00:50:10
% DurationCPUTime: 0.12s
% Computational Cost: add. (246->38), mult. (606->88), div. (0->0), fcn. (489->12), ass. (0->41)
t32 = cos(pkin(6)) * qJD(1) + qJD(2);
t38 = sin(pkin(7));
t41 = cos(pkin(7));
t40 = cos(pkin(12));
t39 = sin(pkin(6));
t60 = qJD(1) * t39;
t54 = t40 * t60;
t67 = t32 * t38 + t41 * t54;
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t37 = sin(pkin(12));
t55 = t37 * t60;
t66 = -t44 * t55 + t67 * t46;
t47 = qJD(3) ^ 2;
t65 = t47 / 0.2e1;
t64 = cos(qJ(5));
t56 = t67 * t44 + t46 * t55;
t23 = qJD(3) * pkin(9) + t56;
t25 = t41 * t32 - t38 * t54;
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t61 = t45 * t23 + t43 * t25;
t16 = qJD(4) * pkin(10) + t61;
t19 = (-pkin(4) * t45 - pkin(10) * t43 - pkin(3)) * qJD(3) - t66;
t42 = sin(qJ(5));
t62 = t64 * t16 + t42 * t19;
t59 = qJD(3) * t43;
t58 = t45 * qJD(3);
t57 = qJD(3) * qJD(4);
t53 = -t42 * t16 + t64 * t19;
t52 = -t43 * t23 + t45 * t25;
t15 = -qJD(4) * pkin(4) - t52;
t48 = qJD(1) ^ 2;
t33 = -qJD(5) + t58;
t29 = t42 * qJD(4) + t64 * t59;
t28 = -t64 * qJD(4) + t42 * t59;
t22 = -qJD(3) * pkin(3) - t66;
t13 = t28 * pkin(5) + qJD(6) + t15;
t12 = -t28 * qJ(6) + t62;
t11 = -t33 * pkin(5) - t29 * qJ(6) + t53;
t1 = [t48 / 0.2e1, t32 ^ 2 / 0.2e1 + (t37 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1) * t48 * t39 ^ 2, t65, t66 * qJD(3), -t56 * qJD(3), t43 ^ 2 * t65, t43 * t47 * t45, t43 * t57, t45 * t57, qJD(4) ^ 2 / 0.2e1, t52 * qJD(4) - t22 * t58, -t61 * qJD(4) + t22 * t59, t29 ^ 2 / 0.2e1, -t29 * t28, -t29 * t33, t28 * t33, t33 ^ 2 / 0.2e1, t15 * t28 - t53 * t33, t15 * t29 + t62 * t33, -t11 * t33 + t13 * t28, t12 * t33 + t13 * t29, -t11 * t29 - t12 * t28, t12 ^ 2 / 0.2e1 + t11 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1;];
T_reg = t1;
