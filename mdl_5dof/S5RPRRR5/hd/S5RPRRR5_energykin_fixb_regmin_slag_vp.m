% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:27
% EndTime: 2019-12-05 18:16:27
% DurationCPUTime: 0.11s
% Computational Cost: add. (108->28), mult. (202->69), div. (0->0), fcn. (107->8), ass. (0->30)
t97 = qJD(1) + qJD(3);
t95 = t97 ^ 2;
t115 = t95 / 0.2e1;
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t98 = sin(pkin(9));
t109 = pkin(1) * qJD(1) * t98;
t99 = cos(pkin(9));
t89 = (pkin(1) * t99 + pkin(2)) * qJD(1);
t113 = t102 * t89 + t105 * t109;
t85 = t97 * pkin(7) + t113;
t114 = t101 * qJD(2) + t104 * t85;
t112 = t101 * t97;
t111 = t104 * t97;
t110 = qJD(4) * t97;
t108 = -t102 * t109 + t105 * t89;
t106 = qJD(1) ^ 2;
t103 = cos(qJ(5));
t100 = sin(qJ(5));
t96 = qJD(4) + qJD(5);
t94 = t104 * qJD(2);
t87 = (t100 * t104 + t101 * t103) * t97;
t86 = t100 * t112 - t103 * t111;
t84 = -t97 * pkin(3) - t108;
t82 = (-pkin(4) * t104 - pkin(3)) * t97 - t108;
t81 = pkin(8) * t111 + t114;
t80 = qJD(4) * pkin(4) + t94 + (-pkin(8) * t97 - t85) * t101;
t1 = [t106 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t98 ^ 2 / 0.2e1 + t99 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t106, t115, t108 * t97, -t113 * t97, t101 ^ 2 * t115, t101 * t95 * t104, t101 * t110, t104 * t110, qJD(4) ^ 2 / 0.2e1, -t84 * t111 + (-t101 * t85 + t94) * qJD(4), -t114 * qJD(4) + t84 * t112, t87 ^ 2 / 0.2e1, -t87 * t86, t87 * t96, -t86 * t96, t96 ^ 2 / 0.2e1, t82 * t86 + (-t100 * t81 + t103 * t80) * t96, t82 * t87 - (t100 * t80 + t103 * t81) * t96;];
T_reg = t1;
