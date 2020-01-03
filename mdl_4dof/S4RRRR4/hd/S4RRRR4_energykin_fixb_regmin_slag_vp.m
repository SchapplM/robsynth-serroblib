% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:14
% EndTime: 2019-12-31 17:26:14
% DurationCPUTime: 0.06s
% Computational Cost: add. (95->27), mult. (245->66), div. (0->0), fcn. (159->6), ass. (0->30)
t114 = -pkin(6) - pkin(5);
t104 = qJD(1) ^ 2;
t113 = t104 / 0.2e1;
t102 = cos(qJ(3));
t100 = sin(qJ(2));
t110 = qJD(1) * t100;
t92 = qJD(2) * pkin(2) + t114 * t110;
t103 = cos(qJ(2));
t109 = qJD(1) * t103;
t93 = t114 * t109;
t99 = sin(qJ(3));
t112 = -t102 * t93 + t99 * t92;
t111 = t103 * t104;
t108 = qJD(1) * qJD(2);
t107 = t100 * t108;
t106 = t103 * t108;
t89 = -t102 * t109 + t99 * t110;
t94 = (-pkin(2) * t103 - pkin(1)) * qJD(1);
t105 = t102 * t92 + t99 * t93;
t101 = cos(qJ(4));
t98 = sin(qJ(4));
t97 = qJD(2) + qJD(3);
t90 = (t100 * t102 + t103 * t99) * qJD(1);
t88 = qJD(4) + t89;
t86 = t101 * t90 + t98 * t97;
t85 = -t101 * t97 + t98 * t90;
t84 = t97 * pkin(7) + t112;
t83 = -t97 * pkin(3) - t105;
t82 = t89 * pkin(3) - t90 * pkin(7) + t94;
t1 = [t113, 0, 0, t100 ^ 2 * t113, t100 * t111, t107, t106, qJD(2) ^ 2 / 0.2e1, pkin(1) * t111 - pkin(5) * t107, -t104 * pkin(1) * t100 - pkin(5) * t106, t90 ^ 2 / 0.2e1, -t90 * t89, t90 * t97, -t89 * t97, t97 ^ 2 / 0.2e1, t105 * t97 + t94 * t89, -t112 * t97 + t94 * t90, t86 ^ 2 / 0.2e1, -t86 * t85, t86 * t88, -t85 * t88, t88 ^ 2 / 0.2e1, (t101 * t82 - t98 * t84) * t88 + t83 * t85, -(t101 * t84 + t98 * t82) * t88 + t83 * t86;];
T_reg = t1;
