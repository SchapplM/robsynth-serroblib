% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPPR7
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
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:40
% EndTime: 2019-12-31 16:41:40
% DurationCPUTime: 0.04s
% Computational Cost: add. (44->18), mult. (113->50), div. (0->0), fcn. (48->4), ass. (0->19)
t63 = qJD(1) ^ 2;
t66 = t63 / 0.2e1;
t59 = sin(pkin(6));
t65 = qJD(1) * t59;
t54 = qJD(1) * qJ(2) + qJD(3);
t53 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t64 = -pkin(5) * qJD(1) + t53;
t62 = cos(qJ(4));
t61 = sin(qJ(4));
t60 = cos(pkin(6));
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = -qJD(1) * pkin(1) + qJD(2);
t51 = pkin(3) * t65 + t54;
t50 = (-t59 * t61 + t60 * t62) * qJD(1);
t49 = (t59 * t62 + t60 * t61) * qJD(1);
t48 = t64 * t60;
t47 = t64 * t59;
t1 = [t66, 0, 0, t55 * qJD(1), t63 * qJ(2), qJ(2) ^ 2 * t66 + t55 ^ 2 / 0.2e1, t54 * t65, t54 * t60 * qJD(1), (-t56 - t57) * t53 * qJD(1), t54 ^ 2 / 0.2e1 + (t56 / 0.2e1 + t57 / 0.2e1) * t53 ^ 2, t50 ^ 2 / 0.2e1, -t50 * t49, t50 * qJD(4), -t49 * qJD(4), qJD(4) ^ 2 / 0.2e1, t51 * t49 + (-t61 * t47 + t62 * t48) * qJD(4), t51 * t50 - (t62 * t47 + t61 * t48) * qJD(4);];
T_reg = t1;
