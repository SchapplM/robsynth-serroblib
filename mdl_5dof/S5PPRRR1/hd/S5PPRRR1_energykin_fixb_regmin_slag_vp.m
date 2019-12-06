% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:50
% EndTime: 2019-12-05 15:12:50
% DurationCPUTime: 0.05s
% Computational Cost: add. (51->16), mult. (121->48), div. (0->0), fcn. (81->8), ass. (0->21)
t74 = qJD(3) + qJD(4);
t73 = t74 ^ 2;
t89 = t73 / 0.2e1;
t75 = sin(pkin(9));
t76 = cos(pkin(9));
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t85 = (-t75 * t79 + t76 * t82) * qJD(1);
t70 = qJD(3) * pkin(3) + t85;
t71 = (t75 * t82 + t76 * t79) * qJD(1);
t78 = sin(qJ(4));
t81 = cos(qJ(4));
t84 = t81 * t70 - t78 * t71;
t88 = (-t74 * pkin(4) - t84) * t74;
t87 = t78 * t70 + t81 * t71;
t86 = qJD(5) * t74;
t83 = qJD(1) ^ 2;
t80 = cos(qJ(5));
t77 = sin(qJ(5));
t67 = t74 * pkin(7) + t87;
t1 = [t83 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t75 ^ 2 / 0.2e1 + t76 ^ 2 / 0.2e1) * t83, qJD(3) ^ 2 / 0.2e1, t85 * qJD(3), -t71 * qJD(3), t89, t84 * t74, -t87 * t74, t77 ^ 2 * t89, t77 * t73 * t80, t77 * t86, t80 * t86, qJD(5) ^ 2 / 0.2e1, (t80 * qJD(2) - t77 * t67) * qJD(5) - t80 * t88, -(t77 * qJD(2) + t80 * t67) * qJD(5) + t77 * t88;];
T_reg = t1;
