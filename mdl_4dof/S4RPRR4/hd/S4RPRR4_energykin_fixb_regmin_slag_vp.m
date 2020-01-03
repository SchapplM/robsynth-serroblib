% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:35
% EndTime: 2019-12-31 16:50:35
% DurationCPUTime: 0.05s
% Computational Cost: add. (44->21), mult. (130->58), div. (0->0), fcn. (65->6), ass. (0->23)
t87 = qJD(1) ^ 2;
t95 = t87 / 0.2e1;
t81 = sin(pkin(7));
t75 = (pkin(1) * t81 + pkin(5)) * qJD(1);
t84 = sin(qJ(3));
t86 = cos(qJ(3));
t94 = t84 * qJD(2) + t86 * t75;
t93 = qJD(1) * t84;
t92 = t86 * qJD(1);
t91 = qJD(1) * qJD(3);
t82 = cos(pkin(7));
t90 = -pkin(1) * t82 - pkin(2);
t89 = t86 * qJD(2) - t84 * t75;
t85 = cos(qJ(4));
t83 = sin(qJ(4));
t77 = -qJD(4) + t92;
t76 = t90 * qJD(1);
t74 = t83 * qJD(3) + t85 * t93;
t73 = -t85 * qJD(3) + t83 * t93;
t71 = (-pkin(3) * t86 - pkin(6) * t84 + t90) * qJD(1);
t70 = qJD(3) * pkin(6) + t94;
t69 = -qJD(3) * pkin(3) - t89;
t1 = [t95, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t81 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t87, t84 ^ 2 * t95, t84 * t87 * t86, t84 * t91, t86 * t91, qJD(3) ^ 2 / 0.2e1, t89 * qJD(3) - t76 * t92, -t94 * qJD(3) + t76 * t93, t74 ^ 2 / 0.2e1, -t74 * t73, -t74 * t77, t73 * t77, t77 ^ 2 / 0.2e1, -(-t83 * t70 + t85 * t71) * t77 + t69 * t73, (t85 * t70 + t83 * t71) * t77 + t69 * t74;];
T_reg = t1;
