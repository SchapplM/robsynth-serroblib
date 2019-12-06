% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:55
% EndTime: 2019-12-05 17:39:55
% DurationCPUTime: 0.06s
% Computational Cost: add. (120->31), mult. (284->73), div. (0->0), fcn. (172->6), ass. (0->29)
t109 = qJD(1) ^ 2;
t114 = t109 / 0.2e1;
t106 = sin(qJ(4));
t108 = cos(qJ(4));
t103 = sin(pkin(8));
t95 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t111 = -pkin(6) * qJD(1) + t95;
t89 = t111 * t103;
t104 = cos(pkin(8));
t90 = t111 * t104;
t113 = t106 * t90 + t108 * t89;
t112 = qJD(1) * t103;
t97 = qJD(1) * qJ(2) + qJD(3);
t93 = pkin(3) * t112 + t97;
t110 = -t106 * t89 + t108 * t90;
t107 = cos(qJ(5));
t105 = sin(qJ(5));
t101 = qJD(4) + qJD(5);
t100 = t104 ^ 2;
t99 = t103 ^ 2;
t98 = -qJD(1) * pkin(1) + qJD(2);
t92 = (-t103 * t106 + t104 * t108) * qJD(1);
t91 = (t103 * t108 + t104 * t106) * qJD(1);
t84 = t91 * pkin(4) + t93;
t83 = -t105 * t91 + t107 * t92;
t82 = t105 * t92 + t107 * t91;
t81 = -t91 * pkin(7) + t113;
t80 = qJD(4) * pkin(4) - t92 * pkin(7) + t110;
t1 = [t114, 0, 0, t98 * qJD(1), t109 * qJ(2), qJ(2) ^ 2 * t114 + t98 ^ 2 / 0.2e1, t97 * t112, t97 * t104 * qJD(1), (-t100 - t99) * t95 * qJD(1), t97 ^ 2 / 0.2e1 + (t99 / 0.2e1 + t100 / 0.2e1) * t95 ^ 2, t92 ^ 2 / 0.2e1, -t92 * t91, t92 * qJD(4), -t91 * qJD(4), qJD(4) ^ 2 / 0.2e1, t110 * qJD(4) + t93 * t91, -t113 * qJD(4) + t93 * t92, t83 ^ 2 / 0.2e1, -t83 * t82, t83 * t101, -t82 * t101, t101 ^ 2 / 0.2e1, t84 * t82 + (-t105 * t81 + t107 * t80) * t101, t84 * t83 - (t105 * t80 + t107 * t81) * t101;];
T_reg = t1;
