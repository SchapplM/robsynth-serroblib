% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:49
% EndTime: 2019-12-05 16:15:49
% DurationCPUTime: 0.05s
% Computational Cost: add. (88->24), mult. (136->55), div. (0->0), fcn. (73->6), ass. (0->22)
t81 = qJD(2) + qJD(3);
t82 = sin(pkin(9));
t93 = t81 * t82;
t83 = cos(pkin(9));
t92 = t81 * t83;
t91 = pkin(2) * qJD(2);
t90 = sin(qJ(3)) * t91;
t77 = t81 * qJ(4) + t90;
t71 = t82 * qJD(1) + t83 * t77;
t89 = cos(qJ(3)) * t91;
t88 = qJD(4) - t89;
t86 = cos(qJ(5));
t84 = sin(qJ(5));
t80 = t83 * qJD(1);
t76 = -t81 * pkin(3) + t88;
t74 = (t82 * t86 + t83 * t84) * t81;
t73 = t84 * t93 - t86 * t92;
t72 = (-pkin(4) * t83 - pkin(3)) * t81 + t88;
t70 = -t82 * t77 + t80;
t69 = pkin(7) * t92 + t71;
t68 = t80 + (-pkin(7) * t81 - t77) * t82;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t81 ^ 2 / 0.2e1, t81 * t89, -t81 * t90, -t76 * t92, t76 * t93, (-t70 * t82 + t71 * t83) * t81, t71 ^ 2 / 0.2e1 + t70 ^ 2 / 0.2e1 + t76 ^ 2 / 0.2e1, t74 ^ 2 / 0.2e1, -t74 * t73, t74 * qJD(5), -t73 * qJD(5), qJD(5) ^ 2 / 0.2e1, t72 * t73 + (t86 * t68 - t84 * t69) * qJD(5), t72 * t74 - (t84 * t68 + t86 * t69) * qJD(5);];
T_reg = t1;
