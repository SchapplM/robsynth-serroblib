% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:02:05
% EndTime: 2020-01-03 12:02:05
% DurationCPUTime: 0.10s
% Computational Cost: add. (120->28), mult. (198->67), div. (0->0), fcn. (107->8), ass. (0->29)
t98 = qJD(1) + qJD(2);
t96 = t98 ^ 2;
t114 = t96 / 0.2e1;
t102 = sin(qJ(4));
t105 = cos(qJ(4));
t100 = cos(pkin(9));
t112 = pkin(1) * qJD(1);
t108 = sin(qJ(2)) * t112;
t107 = cos(qJ(2)) * t112;
t90 = pkin(2) * t98 + t107;
t99 = sin(pkin(9));
t86 = t100 * t108 + t99 * t90;
t84 = pkin(7) * t98 + t86;
t113 = t102 * qJD(3) + t105 * t84;
t111 = t102 * t98;
t110 = t105 * t98;
t109 = qJD(4) * t98;
t85 = t100 * t90 - t99 * t108;
t104 = cos(qJ(5));
t101 = sin(qJ(5));
t97 = qJD(4) + qJD(5);
t95 = t105 * qJD(3);
t88 = (t101 * t105 + t102 * t104) * t98;
t87 = t101 * t111 - t104 * t110;
t83 = -pkin(3) * t98 - t85;
t81 = (-pkin(4) * t105 - pkin(3)) * t98 - t85;
t80 = pkin(8) * t110 + t113;
t79 = qJD(4) * pkin(4) + t95 + (-pkin(8) * t98 - t84) * t102;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t114, t98 * t107, -t98 * t108, t86 ^ 2 / 0.2e1 + t85 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t102 ^ 2 * t114, t102 * t96 * t105, t102 * t109, t105 * t109, qJD(4) ^ 2 / 0.2e1, -t83 * t110 + (-t102 * t84 + t95) * qJD(4), -qJD(4) * t113 + t111 * t83, t88 ^ 2 / 0.2e1, -t88 * t87, t88 * t97, -t87 * t97, t97 ^ 2 / 0.2e1, t81 * t87 + (-t101 * t80 + t104 * t79) * t97, t81 * t88 - (t101 * t79 + t104 * t80) * t97;];
T_reg = t1;
