% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% T_reg [1x13]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:47
% EndTime: 2019-12-31 16:42:47
% DurationCPUTime: 0.04s
% Computational Cost: add. (30->17), mult. (96->47), div. (0->0), fcn. (37->4), ass. (0->17)
t67 = qJD(1) ^ 2;
t74 = t67 / 0.2e1;
t63 = sin(pkin(6));
t59 = (pkin(1) * t63 + pkin(5)) * qJD(1);
t65 = sin(qJ(3));
t66 = cos(qJ(3));
t73 = t65 * qJD(2) + t66 * t59;
t64 = cos(pkin(6));
t69 = -pkin(1) * t64 - pkin(2);
t72 = t67 * t69;
t71 = qJ(4) * qJD(1);
t70 = qJD(1) * qJD(3);
t62 = t66 * qJD(2);
t57 = qJD(4) + (-pkin(3) * t66 + t69) * qJD(1);
t56 = t66 * t71 + t73;
t55 = qJD(3) * pkin(3) + t62 + (-t59 - t71) * t65;
t1 = [t74, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t63 ^ 2 / 0.2e1 + t64 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t67, t65 ^ 2 * t74, t65 * t67 * t66, t65 * t70, t66 * t70, qJD(3) ^ 2 / 0.2e1, -t66 * t72 + (-t65 * t59 + t62) * qJD(3), -t73 * qJD(3) + t65 * t72, (-t55 * t65 + t56 * t66) * qJD(1), t56 ^ 2 / 0.2e1 + t55 ^ 2 / 0.2e1 + t57 ^ 2 / 0.2e1;];
T_reg = t1;
