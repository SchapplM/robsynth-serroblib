% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:53
% EndTime: 2019-12-05 18:03:53
% DurationCPUTime: 0.11s
% Computational Cost: add. (98->30), mult. (250->72), div. (0->0), fcn. (139->6), ass. (0->28)
t105 = qJD(1) ^ 2;
t115 = t105 / 0.2e1;
t114 = cos(qJ(4));
t102 = sin(qJ(4));
t103 = sin(qJ(3));
t100 = sin(pkin(8));
t94 = (pkin(1) * t100 + pkin(6)) * qJD(1);
t104 = cos(qJ(3));
t98 = t104 * qJD(2);
t88 = qJD(3) * pkin(3) + t98 + (-pkin(7) * qJD(1) - t94) * t103;
t110 = qJD(1) * t104;
t112 = t103 * qJD(2) + t104 * t94;
t89 = pkin(7) * t110 + t112;
t113 = t102 * t88 + t114 * t89;
t111 = qJD(1) * t103;
t109 = qJD(1) * qJD(3);
t101 = cos(pkin(8));
t108 = -pkin(1) * t101 - pkin(2);
t107 = -t102 * t89 + t114 * t88;
t92 = (-pkin(3) * t104 + t108) * qJD(1);
t99 = qJD(3) + qJD(4);
t95 = t108 * qJD(1);
t91 = (t102 * t104 + t103 * t114) * qJD(1);
t90 = t102 * t111 - t110 * t114;
t84 = t90 * pkin(4) + qJD(5) + t92;
t83 = -qJ(5) * t90 + t113;
t82 = pkin(4) * t99 - qJ(5) * t91 + t107;
t1 = [t115, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t100 ^ 2 / 0.2e1 + t101 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t105, t103 ^ 2 * t115, t103 * t105 * t104, t103 * t109, t104 * t109, qJD(3) ^ 2 / 0.2e1, -t95 * t110 + (-t103 * t94 + t98) * qJD(3), -qJD(3) * t112 + t111 * t95, t91 ^ 2 / 0.2e1, -t91 * t90, t91 * t99, -t90 * t99, t99 ^ 2 / 0.2e1, t107 * t99 + t92 * t90, -t113 * t99 + t92 * t91, -t82 * t91 - t83 * t90, t83 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1 + t84 ^ 2 / 0.2e1;];
T_reg = t1;
