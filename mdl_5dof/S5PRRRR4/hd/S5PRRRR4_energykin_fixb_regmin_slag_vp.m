% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:58
% EndTime: 2019-12-05 17:07:58
% DurationCPUTime: 0.05s
% Computational Cost: add. (82->23), mult. (130->58), div. (0->0), fcn. (69->6), ass. (0->24)
t84 = qJD(2) + qJD(3);
t82 = t84 ^ 2;
t98 = t82 / 0.2e1;
t86 = sin(qJ(4));
t97 = t84 * t86;
t89 = cos(qJ(4));
t96 = t84 * t89;
t94 = pkin(2) * qJD(2);
t92 = sin(qJ(3)) * t94;
t77 = t84 * pkin(7) + t92;
t95 = t86 * qJD(1) + t89 * t77;
t93 = qJD(4) * t84;
t91 = cos(qJ(3)) * t94;
t88 = cos(qJ(5));
t85 = sin(qJ(5));
t83 = qJD(4) + qJD(5);
t81 = t89 * qJD(1);
t78 = -t84 * pkin(3) - t91;
t75 = -t91 + (-pkin(4) * t89 - pkin(3)) * t84;
t74 = (t85 * t89 + t86 * t88) * t84;
t73 = t85 * t97 - t88 * t96;
t72 = pkin(8) * t96 + t95;
t71 = qJD(4) * pkin(4) + t81 + (-pkin(8) * t84 - t77) * t86;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t98, t84 * t91, -t84 * t92, t86 ^ 2 * t98, t86 * t82 * t89, t86 * t93, t89 * t93, qJD(4) ^ 2 / 0.2e1, -t78 * t96 + (-t86 * t77 + t81) * qJD(4), -t95 * qJD(4) + t78 * t97, t74 ^ 2 / 0.2e1, -t74 * t73, t74 * t83, -t73 * t83, t83 ^ 2 / 0.2e1, t75 * t73 + (t88 * t71 - t85 * t72) * t83, t75 * t74 - (t85 * t71 + t88 * t72) * t83;];
T_reg = t1;
