% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:39
% EndTime: 2019-12-31 16:32:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (32->18), mult. (102->49), div. (0->0), fcn. (57->4), ass. (0->19)
t67 = qJD(2) ^ 2;
t73 = t67 / 0.2e1;
t66 = cos(qJ(3));
t72 = t66 * t67;
t64 = sin(qJ(3));
t69 = qJD(2) * t66;
t71 = pkin(5) * t69 + t64 * qJD(1);
t70 = qJD(2) * t64;
t68 = qJD(2) * qJD(3);
t65 = cos(qJ(4));
t63 = sin(qJ(4));
t62 = qJD(3) + qJD(4);
t61 = t66 * qJD(1);
t57 = (-pkin(3) * t66 - pkin(2)) * qJD(2);
t56 = (t63 * t66 + t64 * t65) * qJD(2);
t55 = t63 * t70 - t65 * t69;
t54 = pkin(6) * t69 + t71;
t53 = qJD(3) * pkin(3) + t61 + (-pkin(6) - pkin(5)) * t70;
t1 = [qJD(1) ^ 2 / 0.2e1, t73, 0, 0, t64 ^ 2 * t73, t64 * t72, t64 * t68, t66 * t68, qJD(3) ^ 2 / 0.2e1, pkin(2) * t72 + (-pkin(5) * t70 + t61) * qJD(3), -t67 * pkin(2) * t64 - t71 * qJD(3), t56 ^ 2 / 0.2e1, -t56 * t55, t56 * t62, -t55 * t62, t62 ^ 2 / 0.2e1, t57 * t55 + (t65 * t53 - t63 * t54) * t62, t57 * t56 - (t63 * t53 + t65 * t54) * t62;];
T_reg = t1;
