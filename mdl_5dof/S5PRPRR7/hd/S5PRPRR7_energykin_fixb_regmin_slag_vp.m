% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:44
% EndTime: 2019-12-05 16:00:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (59->24), mult. (130->59), div. (0->0), fcn. (67->6), ass. (0->25)
t99 = qJD(2) ^ 2;
t107 = t99 / 0.2e1;
t98 = cos(qJ(2));
t100 = -t98 * qJD(1) + qJD(3);
t88 = (-pkin(2) - pkin(6)) * qJD(2) + t100;
t106 = qJD(4) * t88;
t95 = sin(qJ(2));
t104 = t95 * qJD(1);
t90 = qJD(2) * qJ(3) + t104;
t105 = t90 * qJD(2);
t103 = qJD(1) * qJD(2);
t102 = qJD(2) * qJD(4);
t101 = -pkin(7) * qJD(2) + t88;
t97 = cos(qJ(4));
t96 = cos(qJ(5));
t94 = sin(qJ(4));
t93 = sin(qJ(5));
t92 = qJD(4) + qJD(5);
t89 = -qJD(2) * pkin(2) + t100;
t87 = t104 + (pkin(4) * t94 + qJ(3)) * qJD(2);
t86 = (-t93 * t94 + t96 * t97) * qJD(2);
t85 = (t93 * t97 + t94 * t96) * qJD(2);
t84 = t101 * t94;
t83 = qJD(4) * pkin(4) + t101 * t97;
t1 = [qJD(1) ^ 2 / 0.2e1, t107, t98 * t103, -t95 * t103, t89 * qJD(2), t105, t90 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1, t97 ^ 2 * t107, -t97 * t99 * t94, t97 * t102, -t94 * t102, qJD(4) ^ 2 / 0.2e1, t94 * t105 + t97 * t106, t97 * t105 - t94 * t106, t86 ^ 2 / 0.2e1, -t86 * t85, t86 * t92, -t85 * t92, t92 ^ 2 / 0.2e1, (t96 * t83 - t93 * t84) * t92 + t87 * t85, -(t93 * t83 + t96 * t84) * t92 + t87 * t86;];
T_reg = t1;
