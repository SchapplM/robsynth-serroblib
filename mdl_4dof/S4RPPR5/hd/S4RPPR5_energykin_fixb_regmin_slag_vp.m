% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:48
% EndTime: 2019-12-31 16:39:48
% DurationCPUTime: 0.06s
% Computational Cost: add. (41->19), mult. (88->44), div. (0->0), fcn. (27->4), ass. (0->16)
t60 = qJD(1) ^ 2;
t65 = t60 / 0.2e1;
t53 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t57 = cos(pkin(6));
t64 = t57 * t53;
t56 = sin(pkin(6));
t62 = qJ(2) * qJD(1);
t51 = t56 * t53 + t57 * t62;
t63 = (-t64 + (qJ(2) * t56 + pkin(3)) * qJD(1)) * qJD(1);
t61 = qJD(1) * qJD(4);
t59 = cos(qJ(4));
t58 = sin(qJ(4));
t55 = -qJD(1) * pkin(1) + qJD(2);
t50 = -t56 * t62 + t64;
t49 = -qJD(1) * pkin(5) + t51;
t1 = [t65, 0, 0, -t55 * qJD(1), t60 * qJ(2), qJ(2) ^ 2 * t65 + t55 ^ 2 / 0.2e1, -t50 * qJD(1), t51 * qJD(1), t51 ^ 2 / 0.2e1 + t50 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t58 ^ 2 * t65, t58 * t60 * t59, -t58 * t61, -t59 * t61, qJD(4) ^ 2 / 0.2e1, t59 * t63 + (t59 * qJD(3) - t58 * t49) * qJD(4), -t58 * t63 - (t58 * qJD(3) + t59 * t49) * qJD(4);];
T_reg = t1;
