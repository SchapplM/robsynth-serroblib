% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:05
% EndTime: 2022-01-20 12:02:05
% DurationCPUTime: 0.08s
% Computational Cost: add. (155->24), mult. (190->62), div. (0->0), fcn. (101->8), ass. (0->30)
t96 = qJD(1) + qJD(2);
t94 = qJD(3) + t96;
t93 = t94 ^ 2;
t115 = t93 / 0.2e1;
t98 = sin(qJ(4));
t114 = t94 * t98;
t103 = cos(qJ(3));
t112 = pkin(1) * qJD(1);
t108 = sin(qJ(2)) * t112;
t107 = cos(qJ(2)) * t112;
t89 = t96 * pkin(2) + t107;
t99 = sin(qJ(3));
t113 = t103 * t108 + t99 * t89;
t102 = cos(qJ(4));
t111 = t102 * t94;
t110 = qJD(4) * t98;
t109 = qJD(4) * t102;
t85 = t94 * pkin(8) + t113;
t106 = pkin(9) * t94 + t85;
t105 = t103 * t89 - t99 * t108;
t101 = cos(qJ(5));
t97 = sin(qJ(5));
t95 = qJD(4) + qJD(5);
t87 = (t101 * t98 + t102 * t97) * t94;
t86 = -t101 * t111 + t97 * t114;
t84 = -t94 * pkin(3) - t105;
t83 = (-pkin(4) * t102 - pkin(3)) * t94 - t105;
t82 = t106 * t102;
t81 = qJD(4) * pkin(4) - t106 * t98;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t96 ^ 2 / 0.2e1, t96 * t107, -t96 * t108, t115, t105 * t94, -t113 * t94, t98 ^ 2 * t115, t98 * t93 * t102, t94 * t110, t94 * t109, qJD(4) ^ 2 / 0.2e1, -t85 * t110 - t84 * t111, -t85 * t109 + t84 * t114, t87 ^ 2 / 0.2e1, -t87 * t86, t87 * t95, -t86 * t95, t95 ^ 2 / 0.2e1, t83 * t86 + (t101 * t81 - t97 * t82) * t95, t83 * t87 - (t101 * t82 + t97 * t81) * t95;];
T_reg = t1;
