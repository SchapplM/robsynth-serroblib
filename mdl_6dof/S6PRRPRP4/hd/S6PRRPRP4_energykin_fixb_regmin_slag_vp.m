% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:12:46
% EndTime: 2021-01-16 03:12:46
% DurationCPUTime: 0.09s
% Computational Cost: add. (204->43), mult. (444->90), div. (0->0), fcn. (275->8), ass. (0->37)
t38 = qJD(2) ^ 2;
t53 = t38 / 0.2e1;
t52 = -pkin(9) - pkin(3);
t34 = sin(qJ(2));
t49 = qJD(1) * sin(pkin(6));
t24 = qJD(2) * pkin(8) + t34 * t49;
t33 = sin(qJ(3));
t36 = cos(qJ(3));
t48 = qJD(1) * cos(pkin(6));
t40 = -t33 * t24 + t36 * t48;
t39 = qJD(4) - t40;
t46 = t33 * qJD(2);
t13 = pkin(4) * t46 + t52 * qJD(3) + t39;
t42 = -qJ(4) * t33 - pkin(2);
t37 = cos(qJ(2));
t44 = t37 * t49;
t16 = -t44 + (t52 * t36 + t42) * qJD(2);
t32 = sin(qJ(5));
t35 = cos(qJ(5));
t51 = t32 * t13 + t35 * t16;
t50 = t36 * t24 + t33 * t48;
t47 = qJD(2) * t36;
t45 = qJD(2) * qJD(3);
t18 = -qJD(3) * qJ(4) - t50;
t43 = qJD(2) * t49;
t41 = t35 * t13 - t32 * t16;
t14 = pkin(4) * t47 - t18;
t27 = qJD(5) + t46;
t25 = -qJD(2) * pkin(2) - t44;
t23 = t35 * qJD(3) - t32 * t47;
t22 = t32 * qJD(3) + t35 * t47;
t19 = -t44 + (-pkin(3) * t36 + t42) * qJD(2);
t17 = -qJD(3) * pkin(3) + t39;
t10 = t22 * pkin(5) + qJD(6) + t14;
t9 = -t22 * qJ(6) + t51;
t8 = t27 * pkin(5) - t23 * qJ(6) + t41;
t1 = [qJD(1) ^ 2 / 0.2e1, t53, t37 * t43, -t34 * t43, t33 ^ 2 * t53, t33 * t38 * t36, t33 * t45, t36 * t45, qJD(3) ^ 2 / 0.2e1, t40 * qJD(3) - t25 * t47, -t50 * qJD(3) + t25 * t46, (t17 * t33 - t18 * t36) * qJD(2), t17 * qJD(3) + t19 * t47, -t18 * qJD(3) - t19 * t46, t19 ^ 2 / 0.2e1 + t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1, t23 ^ 2 / 0.2e1, -t23 * t22, t23 * t27, -t22 * t27, t27 ^ 2 / 0.2e1, t14 * t22 + t41 * t27, t14 * t23 - t51 * t27, t10 * t22 + t8 * t27, t10 * t23 - t9 * t27, -t9 * t22 - t8 * t23, t9 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1;];
T_reg = t1;
