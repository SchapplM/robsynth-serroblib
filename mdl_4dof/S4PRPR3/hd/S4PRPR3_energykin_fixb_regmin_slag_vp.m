% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:54
% EndTime: 2019-12-31 16:20:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (34->20), mult. (102->47), div. (0->0), fcn. (58->4), ass. (0->17)
t63 = sin(pkin(7));
t64 = cos(pkin(7));
t67 = qJ(3) * qJD(2);
t56 = t63 * qJD(1) + t64 * t67;
t69 = qJD(2) * t63;
t68 = qJD(2) * t64;
t66 = cos(qJ(4));
t65 = sin(qJ(4));
t62 = t64 * qJD(1);
t60 = -qJD(2) * pkin(2) + qJD(3);
t57 = qJD(3) + (-pkin(3) * t64 - pkin(2)) * qJD(2);
t55 = -t63 * t67 + t62;
t54 = (t63 * t66 + t64 * t65) * qJD(2);
t53 = t65 * t69 - t66 * t68;
t52 = pkin(5) * t68 + t56;
t51 = t62 + (-pkin(5) - qJ(3)) * t69;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, -t60 * t68, t60 * t69, (-t55 * t63 + t56 * t64) * qJD(2), t56 ^ 2 / 0.2e1 + t55 ^ 2 / 0.2e1 + t60 ^ 2 / 0.2e1, t54 ^ 2 / 0.2e1, -t54 * t53, t54 * qJD(4), -t53 * qJD(4), qJD(4) ^ 2 / 0.2e1, t57 * t53 + (t66 * t51 - t65 * t52) * qJD(4), t57 * t54 - (t65 * t51 + t66 * t52) * qJD(4);];
T_reg = t1;
