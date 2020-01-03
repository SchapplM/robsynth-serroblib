% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:04
% EndTime: 2019-12-31 18:26:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (88->20), mult. (138->50), div. (0->0), fcn. (51->6), ass. (0->23)
t82 = -qJD(1) + qJD(3);
t81 = t82 ^ 2;
t95 = t81 / 0.2e1;
t89 = qJD(1) ^ 2;
t94 = t89 / 0.2e1;
t79 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t86 = sin(qJ(3));
t88 = cos(qJ(3));
t91 = qJ(2) * qJD(1);
t90 = t88 * t79 - t86 * t91;
t75 = t82 * pkin(3) + t90;
t77 = t86 * t79 + t88 * t91;
t83 = sin(pkin(8));
t84 = cos(pkin(8));
t72 = t84 * t75 - t83 * t77;
t93 = (-t82 * pkin(4) - t72) * t82;
t73 = t83 * t75 + t84 * t77;
t92 = qJD(5) * t82;
t87 = cos(qJ(5));
t85 = sin(qJ(5));
t80 = -qJD(1) * pkin(1) + qJD(2);
t71 = t82 * pkin(7) + t73;
t1 = [t94, 0, 0, -t80 * qJD(1), t89 * qJ(2), qJ(2) ^ 2 * t94 + t80 ^ 2 / 0.2e1, t95, t90 * t82, -t77 * t82, t73 ^ 2 / 0.2e1 + t72 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t85 ^ 2 * t95, t85 * t81 * t87, t85 * t92, t87 * t92, qJD(5) ^ 2 / 0.2e1, -t87 * t93 + (t87 * qJD(4) - t85 * t71) * qJD(5), t85 * t93 - (t85 * qJD(4) + t87 * t71) * qJD(5);];
T_reg = t1;
