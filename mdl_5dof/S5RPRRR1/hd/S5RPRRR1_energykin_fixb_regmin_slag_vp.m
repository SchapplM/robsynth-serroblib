% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:09:50
% EndTime: 2019-12-05 18:09:50
% DurationCPUTime: 0.05s
% Computational Cost: add. (59->22), mult. (175->67), div. (0->0), fcn. (118->6), ass. (0->23)
t82 = qJD(1) ^ 2;
t89 = t82 / 0.2e1;
t78 = sin(qJ(3));
t88 = qJD(1) * t78;
t87 = qJ(2) * qJD(1);
t86 = qJ(2) * qJD(3);
t85 = qJD(1) * qJD(3);
t84 = t78 * t87;
t81 = cos(qJ(3));
t83 = t81 * t87;
t77 = sin(qJ(4));
t80 = cos(qJ(4));
t70 = -t80 * qJD(3) + t77 * t88;
t79 = cos(qJ(5));
t76 = sin(qJ(5));
t73 = t81 * qJD(1) - qJD(4);
t71 = t77 * qJD(3) + t80 * t88;
t69 = t77 * qJD(2) + t80 * t83;
t68 = -t80 * qJD(2) + t77 * t83;
t67 = qJD(5) + t70;
t66 = t79 * t71 - t76 * t73;
t65 = t76 * t71 + t79 * t73;
t1 = [t89, 0, 0, -qJD(2) * qJD(1), t82 * qJ(2), qJ(2) ^ 2 * t89 + qJD(2) ^ 2 / 0.2e1, t78 ^ 2 * t89, t78 * t82 * t81, t78 * t85, t81 * t85, qJD(3) ^ 2 / 0.2e1, (-qJD(2) * t81 - t78 * t86) * qJD(1), (qJD(2) * t78 - t81 * t86) * qJD(1), t71 ^ 2 / 0.2e1, -t71 * t70, -t71 * t73, t70 * t73, t73 ^ 2 / 0.2e1, t68 * t73 + t70 * t84, t69 * t73 + t71 * t84, t66 ^ 2 / 0.2e1, -t66 * t65, t66 * t67, -t65 * t67, t67 ^ 2 / 0.2e1, (-t76 * t69 + t79 * t84) * t67 + t68 * t65, -(t79 * t69 + t76 * t84) * t67 + t68 * t66;];
T_reg = t1;
