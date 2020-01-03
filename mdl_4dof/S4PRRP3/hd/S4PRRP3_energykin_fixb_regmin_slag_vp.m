% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% T_reg [1x13]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:55
% EndTime: 2019-12-31 16:26:55
% DurationCPUTime: 0.03s
% Computational Cost: add. (19->13), mult. (66->38), div. (0->0), fcn. (26->2), ass. (0->14)
t51 = qJD(2) ^ 2;
t57 = t51 / 0.2e1;
t50 = cos(qJ(3));
t56 = t50 * t51;
t49 = sin(qJ(3));
t53 = qJD(2) * t50;
t55 = pkin(5) * t53 + t49 * qJD(1);
t54 = qJD(2) * t49;
t52 = qJD(2) * qJD(3);
t48 = t50 * qJD(1);
t45 = qJD(4) + (-pkin(3) * t50 - pkin(2)) * qJD(2);
t44 = qJ(4) * t53 + t55;
t43 = qJD(3) * pkin(3) + t48 + (-pkin(5) - qJ(4)) * t54;
t1 = [qJD(1) ^ 2 / 0.2e1, t57, 0, 0, t49 ^ 2 * t57, t49 * t56, t49 * t52, t50 * t52, qJD(3) ^ 2 / 0.2e1, pkin(2) * t56 + (-pkin(5) * t54 + t48) * qJD(3), -t51 * pkin(2) * t49 - t55 * qJD(3), (-t43 * t49 + t44 * t50) * qJD(2), t44 ^ 2 / 0.2e1 + t43 ^ 2 / 0.2e1 + t45 ^ 2 / 0.2e1;];
T_reg = t1;
