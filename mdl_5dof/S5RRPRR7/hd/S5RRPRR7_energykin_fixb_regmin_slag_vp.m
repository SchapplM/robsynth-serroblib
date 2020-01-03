% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:40
% EndTime: 2019-12-31 20:15:40
% DurationCPUTime: 0.05s
% Computational Cost: add. (110->25), mult. (146->59), div. (0->0), fcn. (67->6), ass. (0->25)
t83 = qJD(1) + qJD(2);
t81 = t83 ^ 2;
t98 = t81 / 0.2e1;
t96 = pkin(1) * qJD(1);
t93 = sin(qJ(2)) * t96;
t79 = t83 * qJ(3) + t93;
t97 = t79 * t83;
t85 = sin(qJ(4));
t95 = qJD(4) * t85;
t88 = cos(qJ(4));
t94 = qJD(4) * t88;
t92 = cos(qJ(2)) * t96;
t90 = qJD(3) - t92;
t77 = (-pkin(2) - pkin(7)) * t83 + t90;
t91 = -pkin(8) * t83 + t77;
t87 = cos(qJ(5));
t84 = sin(qJ(5));
t82 = qJD(4) + qJD(5);
t78 = -t83 * pkin(2) + t90;
t76 = t93 + (pkin(4) * t85 + qJ(3)) * t83;
t75 = (-t84 * t85 + t87 * t88) * t83;
t74 = (t84 * t88 + t85 * t87) * t83;
t73 = t91 * t85;
t72 = qJD(4) * pkin(4) + t91 * t88;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t98, t83 * t92, -t83 * t93, t78 * t83, t97, t79 ^ 2 / 0.2e1 + t78 ^ 2 / 0.2e1, t88 ^ 2 * t98, -t88 * t81 * t85, t83 * t94, -t83 * t95, qJD(4) ^ 2 / 0.2e1, t77 * t94 + t85 * t97, -t77 * t95 + t88 * t97, t75 ^ 2 / 0.2e1, -t75 * t74, t75 * t82, -t74 * t82, t82 ^ 2 / 0.2e1, t76 * t74 + (t87 * t72 - t84 * t73) * t82, t76 * t75 - (t84 * t72 + t87 * t73) * t82;];
T_reg = t1;
