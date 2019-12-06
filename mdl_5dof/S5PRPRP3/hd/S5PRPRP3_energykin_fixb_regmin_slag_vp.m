% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x14]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:40
% EndTime: 2019-12-05 15:33:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (53->22), mult. (131->55), div. (0->0), fcn. (72->6), ass. (0->23)
t101 = qJD(2) ^ 2;
t108 = t101 / 0.2e1;
t98 = sin(qJ(2));
t106 = qJD(1) * t98;
t100 = cos(qJ(2));
t90 = qJD(2) * pkin(2) + t100 * qJD(1);
t95 = sin(pkin(8));
t96 = cos(pkin(8));
t88 = t96 * t106 + t95 * t90;
t86 = qJD(2) * pkin(6) + t88;
t97 = sin(qJ(4));
t99 = cos(qJ(4));
t107 = t97 * qJD(3) + t99 * t86;
t87 = -t95 * t106 + t96 * t90;
t105 = qJD(2) * (-qJD(2) * pkin(3) - t87);
t104 = qJ(5) * qJD(2);
t103 = qJD(1) * qJD(2);
t102 = qJD(2) * qJD(4);
t94 = t99 * qJD(3);
t83 = qJD(5) + (-pkin(4) * t99 - pkin(3)) * qJD(2) - t87;
t82 = t99 * t104 + t107;
t81 = qJD(4) * pkin(4) + t94 + (-t86 - t104) * t97;
t1 = [qJD(1) ^ 2 / 0.2e1, t108, t100 * t103, -t98 * t103, t88 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t97 ^ 2 * t108, t97 * t101 * t99, t97 * t102, t99 * t102, qJD(4) ^ 2 / 0.2e1, (-t97 * t86 + t94) * qJD(4) - t99 * t105, -t107 * qJD(4) + t97 * t105, (-t81 * t97 + t82 * t99) * qJD(2), t82 ^ 2 / 0.2e1 + t81 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1;];
T_reg = t1;
