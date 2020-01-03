% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRP6
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
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:13
% EndTime: 2019-12-31 16:46:13
% DurationCPUTime: 0.03s
% Computational Cost: add. (30->15), mult. (75->39), div. (0->0), fcn. (20->2), ass. (0->14)
t51 = qJD(1) ^ 2;
t56 = t51 / 0.2e1;
t55 = t51 * qJ(2);
t47 = qJD(2) + (-pkin(1) - pkin(5)) * qJD(1);
t54 = qJD(3) * t47;
t53 = qJD(1) * qJD(3);
t52 = -qJ(4) * qJD(1) + t47;
t50 = cos(qJ(3));
t49 = sin(qJ(3));
t48 = -qJD(1) * pkin(1) + qJD(2);
t46 = qJD(4) + (pkin(3) * t49 + qJ(2)) * qJD(1);
t45 = t52 * t49;
t44 = qJD(3) * pkin(3) + t52 * t50;
t1 = [t56, 0, 0, t48 * qJD(1), t55, qJ(2) ^ 2 * t56 + t48 ^ 2 / 0.2e1, t50 ^ 2 * t56, -t50 * t51 * t49, t50 * t53, -t49 * t53, qJD(3) ^ 2 / 0.2e1, t49 * t55 + t50 * t54, -t49 * t54 + t50 * t55, (-t44 * t50 - t45 * t49) * qJD(1), t45 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1 + t46 ^ 2 / 0.2e1;];
T_reg = t1;
