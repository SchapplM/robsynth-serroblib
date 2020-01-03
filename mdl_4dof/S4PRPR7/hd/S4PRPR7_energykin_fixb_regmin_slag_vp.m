% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% T_reg [1x14]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:54
% EndTime: 2019-12-31 16:25:54
% DurationCPUTime: 0.03s
% Computational Cost: add. (21->11), mult. (53->32), div. (0->0), fcn. (19->4), ass. (0->14)
t61 = qJD(2) ^ 2;
t67 = t61 / 0.2e1;
t60 = cos(qJ(2));
t62 = -t60 * qJD(1) + qJD(3);
t66 = qJD(4) * ((-pkin(2) - pkin(5)) * qJD(2) + t62);
t58 = sin(qJ(2));
t55 = qJD(2) * qJ(3) + t58 * qJD(1);
t65 = t55 * qJD(2);
t64 = qJD(1) * qJD(2);
t63 = qJD(2) * qJD(4);
t59 = cos(qJ(4));
t57 = sin(qJ(4));
t54 = -qJD(2) * pkin(2) + t62;
t1 = [qJD(1) ^ 2 / 0.2e1, t67, t60 * t64, -t58 * t64, t54 * qJD(2), t65, t55 ^ 2 / 0.2e1 + t54 ^ 2 / 0.2e1, t59 ^ 2 * t67, -t59 * t61 * t57, t59 * t63, -t57 * t63, qJD(4) ^ 2 / 0.2e1, t57 * t65 + t59 * t66, -t57 * t66 + t59 * t65;];
T_reg = t1;
