% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% T_reg [1x12]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:51
% EndTime: 2019-12-31 16:19:51
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->8), mult. (47->28), div. (0->0), fcn. (21->4), ass. (0->13)
t53 = qJD(3) ^ 2;
t59 = t53 / 0.2e1;
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t58 = t52 * qJD(1) + t50 * qJD(2);
t54 = -t50 * qJD(1) + t52 * qJD(2);
t57 = qJD(3) * (-qJD(3) * pkin(3) - t54);
t56 = qJD(4) * (qJD(3) * pkin(5) + t58);
t55 = qJD(3) * qJD(4);
t51 = cos(qJ(4));
t49 = sin(qJ(4));
t48 = qJD(1) ^ 2 / 0.2e1;
t1 = [t48, t48 + qJD(2) ^ 2 / 0.2e1, t59, t54 * qJD(3), -t58 * qJD(3), t49 ^ 2 * t59, t49 * t53 * t51, t49 * t55, t51 * t55, qJD(4) ^ 2 / 0.2e1, -t49 * t56 - t51 * t57, t49 * t57 - t51 * t56;];
T_reg = t1;
