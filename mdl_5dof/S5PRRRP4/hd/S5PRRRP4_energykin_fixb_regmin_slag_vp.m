% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:19
% EndTime: 2019-12-05 16:46:19
% DurationCPUTime: 0.05s
% Computational Cost: add. (99->21), mult. (161->56), div. (0->0), fcn. (87->6), ass. (0->24)
t97 = qJD(2) + qJD(3);
t96 = t97 ^ 2;
t112 = t96 / 0.2e1;
t98 = sin(qJ(4));
t111 = t97 * t98;
t102 = cos(qJ(3));
t100 = sin(qJ(2));
t107 = qJD(1) * t100;
t103 = cos(qJ(2));
t93 = qJD(2) * pkin(2) + t103 * qJD(1);
t99 = sin(qJ(3));
t110 = t102 * t107 + t99 * t93;
t101 = cos(qJ(4));
t109 = t101 * t97;
t108 = qJD(4) * t98;
t106 = qJD(4) * t101;
t105 = qJD(1) * qJD(2);
t104 = t102 * t93 - t99 * t107;
t91 = t97 * pkin(7) + t110;
t90 = -t97 * pkin(3) - t104;
t89 = qJD(4) * qJ(5) + t101 * t91;
t88 = -qJD(4) * pkin(4) + t98 * t91 + qJD(5);
t87 = (-pkin(4) * t101 - qJ(5) * t98 - pkin(3)) * t97 - t104;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t103 * t105, -t100 * t105, t112, t104 * t97, -t110 * t97, t98 ^ 2 * t112, t98 * t96 * t101, t97 * t108, t97 * t106, qJD(4) ^ 2 / 0.2e1, -t91 * t108 - t90 * t109, -t91 * t106 + t90 * t111, -t88 * qJD(4) - t87 * t109, (t101 * t89 + t88 * t98) * t97, t89 * qJD(4) - t87 * t111, t89 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1;];
T_reg = t1;
