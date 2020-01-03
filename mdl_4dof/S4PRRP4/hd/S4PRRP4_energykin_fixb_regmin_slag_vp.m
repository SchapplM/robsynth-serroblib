% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRP4
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
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:57
% EndTime: 2019-12-31 16:27:58
% DurationCPUTime: 0.03s
% Computational Cost: add. (30->14), mult. (86->42), div. (0->0), fcn. (33->2), ass. (0->14)
t53 = qJD(2) ^ 2;
t60 = t53 / 0.2e1;
t52 = cos(qJ(3));
t59 = t52 * t53;
t51 = sin(qJ(3));
t56 = qJD(2) * t52;
t58 = pkin(5) * t56 + t51 * qJD(1);
t57 = qJD(2) * t51;
t55 = qJD(2) * qJD(3);
t54 = -pkin(5) * t57 + t52 * qJD(1);
t48 = qJD(3) * qJ(4) + t58;
t47 = (-pkin(3) * t52 - qJ(4) * t51 - pkin(2)) * qJD(2);
t46 = -qJD(3) * pkin(3) + qJD(4) - t54;
t1 = [qJD(1) ^ 2 / 0.2e1, t60, 0, 0, t51 ^ 2 * t60, t51 * t59, t51 * t55, t52 * t55, qJD(3) ^ 2 / 0.2e1, pkin(2) * t59 + t54 * qJD(3), -t53 * pkin(2) * t51 - t58 * qJD(3), -t46 * qJD(3) - t47 * t56, (t46 * t51 + t48 * t52) * qJD(2), t48 * qJD(3) - t47 * t57, t48 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1;];
T_reg = t1;
