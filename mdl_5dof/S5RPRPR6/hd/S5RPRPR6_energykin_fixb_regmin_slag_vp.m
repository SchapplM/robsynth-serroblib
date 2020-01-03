% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:48
% EndTime: 2019-12-31 18:17:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (67->20), mult. (131->47), div. (0->0), fcn. (53->6), ass. (0->22)
t69 = qJD(1) + qJD(3);
t68 = t69 ^ 2;
t85 = t68 / 0.2e1;
t72 = cos(pkin(8));
t66 = (pkin(1) * t72 + pkin(2)) * qJD(1);
t74 = sin(qJ(3));
t76 = cos(qJ(3));
t71 = sin(pkin(8));
t82 = pkin(1) * qJD(1) * t71;
t79 = t74 * t66 + t76 * t82;
t64 = t69 * qJ(4) + t79;
t84 = t64 * t69;
t83 = qJD(5) * t69;
t81 = t76 * t66 - t74 * t82;
t80 = qJD(4) - t81;
t77 = qJD(1) ^ 2;
t75 = cos(qJ(5));
t73 = sin(qJ(5));
t70 = qJD(2) ^ 2 / 0.2e1;
t63 = -t69 * pkin(3) + t80;
t62 = (-pkin(3) - pkin(7)) * t69 + t80;
t1 = [t77 / 0.2e1, 0, 0, t70 + (t71 ^ 2 / 0.2e1 + t72 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t77, t85, t81 * t69, -t79 * t69, t63 * t69, t84, t70 + t64 ^ 2 / 0.2e1 + t63 ^ 2 / 0.2e1, t75 ^ 2 * t85, -t75 * t68 * t73, t75 * t83, -t73 * t83, qJD(5) ^ 2 / 0.2e1, t73 * t84 + (-t73 * qJD(2) + t75 * t62) * qJD(5), t75 * t84 - (t75 * qJD(2) + t73 * t62) * qJD(5);];
T_reg = t1;
