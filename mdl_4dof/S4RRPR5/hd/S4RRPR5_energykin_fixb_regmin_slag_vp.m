% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:32
% EndTime: 2019-12-31 17:03:32
% DurationCPUTime: 0.03s
% Computational Cost: add. (42->12), mult. (63->32), div. (0->0), fcn. (19->4), ass. (0->16)
t52 = qJD(1) + qJD(2);
t51 = t52 ^ 2;
t64 = t51 / 0.2e1;
t62 = pkin(1) * qJD(1);
t59 = sin(qJ(2)) * t62;
t49 = t52 * qJ(3) + t59;
t63 = t49 * t52;
t53 = sin(qJ(4));
t61 = qJD(4) * t53;
t55 = cos(qJ(4));
t60 = qJD(4) * t55;
t58 = cos(qJ(2)) * t62;
t57 = qJD(3) - t58;
t48 = -t52 * pkin(2) + t57;
t47 = (-pkin(2) - pkin(6)) * t52 + t57;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t64, t52 * t58, -t52 * t59, t48 * t52, t63, t49 ^ 2 / 0.2e1 + t48 ^ 2 / 0.2e1, t55 ^ 2 * t64, -t55 * t51 * t53, t52 * t60, -t52 * t61, qJD(4) ^ 2 / 0.2e1, t47 * t60 + t53 * t63, -t47 * t61 + t55 * t63;];
T_reg = t1;
