% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:16
% EndTime: 2019-12-31 16:47:16
% DurationCPUTime: 0.03s
% Computational Cost: add. (43->18), mult. (95->45), div. (0->0), fcn. (25->2), ass. (0->14)
t53 = qJD(1) ^ 2;
t58 = t53 / 0.2e1;
t57 = t53 * qJ(2);
t51 = sin(qJ(3));
t52 = cos(qJ(3));
t47 = (pkin(3) * t51 - qJ(4) * t52 + qJ(2)) * qJD(1);
t56 = qJD(1) * t47;
t49 = qJD(2) + (-pkin(1) - pkin(5)) * qJD(1);
t55 = qJD(3) * t49;
t54 = qJD(1) * qJD(3);
t50 = -qJD(1) * pkin(1) + qJD(2);
t48 = qJD(3) * qJ(4) + t51 * t49;
t46 = -qJD(3) * pkin(3) - t52 * t49 + qJD(4);
t1 = [t58, 0, 0, t50 * qJD(1), t57, qJ(2) ^ 2 * t58 + t50 ^ 2 / 0.2e1, t52 ^ 2 * t58, -t52 * t53 * t51, t52 * t54, -t51 * t54, qJD(3) ^ 2 / 0.2e1, t51 * t57 + t52 * t55, -t51 * t55 + t52 * t57, -t46 * qJD(3) + t51 * t56, (t46 * t52 - t48 * t51) * qJD(1), t48 * qJD(3) - t52 * t56, t48 ^ 2 / 0.2e1 + t47 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1;];
T_reg = t1;
