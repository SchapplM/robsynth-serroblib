% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:20
% EndTime: 2020-01-03 11:36:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (113->26), mult. (212->64), div. (0->0), fcn. (108->8), ass. (0->28)
t100 = sin(pkin(9));
t99 = qJD(1) + qJD(3);
t97 = t99 ^ 2;
t118 = t97 * t100 ^ 2;
t105 = sin(qJ(3));
t107 = cos(qJ(3));
t101 = sin(pkin(8));
t114 = pkin(1) * qJD(1) * t101;
t103 = cos(pkin(8));
t91 = (pkin(1) * t103 + pkin(2)) * qJD(1);
t117 = t105 * t91 + t107 * t114;
t116 = t100 * t99;
t102 = cos(pkin(9));
t115 = t102 * t99;
t104 = sin(qJ(5));
t113 = t104 * t116;
t106 = cos(qJ(5));
t112 = t106 * t116;
t111 = -t105 * t114 + t107 * t91;
t110 = qJD(4) - t111;
t108 = qJD(1) ^ 2;
t92 = -qJD(5) + t115;
t89 = t99 * qJ(4) + t117;
t88 = -t99 * pkin(3) + t110;
t87 = t100 * qJD(2) + t102 * t89;
t85 = -t102 * qJD(2) + t100 * t89;
t84 = (-pkin(4) * t102 - pkin(7) * t100 - pkin(3)) * t99 + t110;
t1 = [t108 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t101 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t108, t97 / 0.2e1, t111 * t99, -t117 * t99, -t88 * t115, t88 * t116, (t100 * t85 + t102 * t87) * t99, t87 ^ 2 / 0.2e1 + t85 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1, t106 ^ 2 * t118 / 0.2e1, -t106 * t104 * t118, -t92 * t112, t92 * t113, t92 ^ 2 / 0.2e1, -(-t104 * t87 + t106 * t84) * t92 + t85 * t113, (t104 * t84 + t106 * t87) * t92 + t85 * t112;];
T_reg = t1;
