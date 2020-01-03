% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:53
% EndTime: 2019-12-31 18:05:53
% DurationCPUTime: 0.05s
% Computational Cost: add. (59->25), mult. (123->60), div. (0->0), fcn. (47->4), ass. (0->23)
t88 = qJD(1) ^ 2;
t95 = t88 / 0.2e1;
t94 = -pkin(1) - qJ(3);
t87 = cos(qJ(4));
t93 = qJD(1) * t87;
t80 = qJD(1) * qJ(2) + qJD(3);
t76 = -qJD(1) * pkin(6) + t80;
t92 = qJD(4) * t76;
t77 = -t94 * qJD(1) - qJD(2);
t91 = t77 * qJD(1);
t85 = sin(qJ(4));
t90 = t85 * qJD(1);
t89 = qJD(1) * qJD(4);
t86 = cos(qJ(5));
t84 = sin(qJ(5));
t81 = -qJD(1) * pkin(1) + qJD(2);
t79 = qJD(5) + t90;
t75 = t84 * qJD(4) + t86 * t93;
t74 = -t86 * qJD(4) + t84 * t93;
t73 = -qJD(4) * pkin(4) - t87 * t76;
t72 = qJD(4) * pkin(7) + t85 * t76;
t71 = -qJD(2) + (pkin(4) * t85 - pkin(7) * t87 - t94) * qJD(1);
t1 = [t95, 0, 0, t81 * qJD(1), t88 * qJ(2), qJ(2) ^ 2 * t95 + t81 ^ 2 / 0.2e1, t80 * qJD(1), t91, t77 ^ 2 / 0.2e1 + t80 ^ 2 / 0.2e1, t87 ^ 2 * t95, -t87 * t88 * t85, t87 * t89, -t85 * t89, qJD(4) ^ 2 / 0.2e1, t77 * t90 + t87 * t92, -t85 * t92 + t87 * t91, t75 ^ 2 / 0.2e1, -t75 * t74, t75 * t79, -t74 * t79, t79 ^ 2 / 0.2e1, (t86 * t71 - t84 * t72) * t79 + t73 * t74, -(t84 * t71 + t86 * t72) * t79 + t73 * t75;];
T_reg = t1;
