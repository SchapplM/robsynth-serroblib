% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:05
% EndTime: 2019-12-31 19:06:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (116->26), mult. (169->63), div. (0->0), fcn. (75->6), ass. (0->30)
t96 = -qJD(1) + qJD(3);
t94 = t96 ^ 2;
t113 = t94 / 0.2e1;
t103 = qJD(1) ^ 2;
t112 = t103 / 0.2e1;
t98 = sin(qJ(4));
t111 = t96 * t98;
t102 = cos(qJ(3));
t107 = qJ(2) * qJD(1);
t90 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t99 = sin(qJ(3));
t110 = t102 * t107 + t99 * t90;
t101 = cos(qJ(4));
t109 = t101 * t96;
t108 = qJD(4) * t98;
t106 = qJD(4) * t101;
t85 = t96 * pkin(7) + t110;
t105 = pkin(8) * t96 + t85;
t104 = t102 * t90 - t99 * t107;
t100 = cos(qJ(5));
t97 = sin(qJ(5));
t95 = qJD(4) + qJD(5);
t93 = -qJD(1) * pkin(1) + qJD(2);
t87 = (t100 * t98 + t101 * t97) * t96;
t86 = -t100 * t109 + t97 * t111;
t84 = -t96 * pkin(3) - t104;
t83 = (-pkin(4) * t101 - pkin(3)) * t96 - t104;
t82 = t105 * t101;
t81 = qJD(4) * pkin(4) - t105 * t98;
t1 = [t112, 0, 0, -t93 * qJD(1), t103 * qJ(2), qJ(2) ^ 2 * t112 + t93 ^ 2 / 0.2e1, t113, t104 * t96, -t110 * t96, t98 ^ 2 * t113, t98 * t94 * t101, t96 * t108, t96 * t106, qJD(4) ^ 2 / 0.2e1, -t85 * t108 - t84 * t109, -t85 * t106 + t84 * t111, t87 ^ 2 / 0.2e1, -t87 * t86, t87 * t95, -t86 * t95, t95 ^ 2 / 0.2e1, t83 * t86 + (t100 * t81 - t97 * t82) * t95, t83 * t87 - (t100 * t82 + t97 * t81) * t95;];
T_reg = t1;
