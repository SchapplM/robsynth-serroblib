% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:41
% EndTime: 2019-12-31 18:16:41
% DurationCPUTime: 0.05s
% Computational Cost: add. (90->30), mult. (171->65), div. (0->0), fcn. (50->2), ass. (0->21)
t75 = qJD(1) ^ 2;
t84 = t75 / 0.2e1;
t83 = -pkin(3) - pkin(4);
t70 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t73 = sin(qJ(3));
t68 = qJD(3) * qJ(4) + t73 * t70;
t82 = t75 * qJ(2);
t81 = qJD(1) * t73;
t74 = cos(qJ(3));
t80 = qJD(1) * t74;
t79 = qJD(3) * t70;
t78 = qJ(5) * qJD(1);
t77 = qJD(1) * qJD(3);
t76 = qJ(4) * t74 - qJ(2);
t71 = -qJD(1) * pkin(1) + qJD(2);
t67 = (pkin(3) * t73 - t76) * qJD(1);
t66 = -qJD(3) * pkin(3) - t74 * t70 + qJD(4);
t65 = t73 * t78 + t68;
t64 = qJD(5) + (t83 * t73 + t76) * qJD(1);
t63 = qJD(4) + (-t70 - t78) * t74 + t83 * qJD(3);
t1 = [t84, 0, 0, t71 * qJD(1), t82, qJ(2) ^ 2 * t84 + t71 ^ 2 / 0.2e1, t74 ^ 2 * t84, -t74 * t75 * t73, t74 * t77, -t73 * t77, qJD(3) ^ 2 / 0.2e1, t73 * t82 + t74 * t79, -t73 * t79 + t74 * t82, -t66 * qJD(3) + t67 * t81, (t66 * t74 - t68 * t73) * qJD(1), t68 * qJD(3) - t67 * t80, t68 ^ 2 / 0.2e1 + t67 ^ 2 / 0.2e1 + t66 ^ 2 / 0.2e1, -t63 * qJD(3) - t64 * t81, t65 * qJD(3) + t64 * t80, (-t63 * t74 + t65 * t73) * qJD(1), t65 ^ 2 / 0.2e1 + t63 ^ 2 / 0.2e1 + t64 ^ 2 / 0.2e1;];
T_reg = t1;
