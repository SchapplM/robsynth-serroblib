% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (43->20), mult. (111->51), div. (0->0), fcn. (51->4), ass. (0->19)
t63 = qJD(1) ^ 2;
t68 = t63 / 0.2e1;
t67 = t63 * qJ(2);
t56 = qJD(2) + (-pkin(1) - pkin(5)) * qJD(1);
t66 = qJD(3) * t56;
t65 = qJD(1) * qJD(3);
t64 = -pkin(6) * qJD(1) + t56;
t62 = cos(qJ(3));
t61 = cos(qJ(4));
t60 = sin(qJ(3));
t59 = sin(qJ(4));
t58 = qJD(3) + qJD(4);
t57 = -qJD(1) * pkin(1) + qJD(2);
t55 = (pkin(3) * t60 + qJ(2)) * qJD(1);
t54 = (-t59 * t60 + t61 * t62) * qJD(1);
t53 = (t59 * t62 + t60 * t61) * qJD(1);
t52 = t64 * t60;
t51 = qJD(3) * pkin(3) + t62 * t64;
t1 = [t68, 0, 0, t57 * qJD(1), t67, qJ(2) ^ 2 * t68 + t57 ^ 2 / 0.2e1, t62 ^ 2 * t68, -t62 * t63 * t60, t62 * t65, -t60 * t65, qJD(3) ^ 2 / 0.2e1, t60 * t67 + t62 * t66, -t60 * t66 + t62 * t67, t54 ^ 2 / 0.2e1, -t54 * t53, t54 * t58, -t53 * t58, t58 ^ 2 / 0.2e1, t55 * t53 + (t51 * t61 - t52 * t59) * t58, t55 * t54 - (t51 * t59 + t52 * t61) * t58;];
T_reg = t1;
