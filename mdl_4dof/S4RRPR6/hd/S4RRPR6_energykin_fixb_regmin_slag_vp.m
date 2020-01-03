% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:53
% EndTime: 2019-12-31 17:04:53
% DurationCPUTime: 0.06s
% Computational Cost: add. (89->26), mult. (257->63), div. (0->0), fcn. (166->6), ass. (0->27)
t111 = qJD(1) * (pkin(5) + qJ(3));
t104 = qJD(1) ^ 2;
t110 = t104 / 0.2e1;
t101 = sin(qJ(2));
t94 = qJD(2) * pkin(2) - t101 * t111;
t103 = cos(qJ(2));
t95 = t103 * t111;
t98 = sin(pkin(7));
t99 = cos(pkin(7));
t86 = t98 * t94 + t99 * t95;
t108 = t103 * t104;
t107 = qJD(1) * qJD(2);
t85 = t99 * t94 - t98 * t95;
t106 = t101 * t107;
t105 = t103 * t107;
t96 = qJD(3) + (-pkin(2) * t103 - pkin(1)) * qJD(1);
t102 = cos(qJ(4));
t100 = sin(qJ(4));
t97 = qJD(2) + qJD(4);
t92 = (t101 * t99 + t103 * t98) * qJD(1);
t91 = (-t101 * t98 + t103 * t99) * qJD(1);
t87 = -t91 * pkin(3) + t96;
t84 = t100 * t91 + t102 * t92;
t83 = t100 * t92 - t102 * t91;
t82 = t91 * pkin(6) + t86;
t81 = qJD(2) * pkin(3) - t92 * pkin(6) + t85;
t1 = [t110, 0, 0, t101 ^ 2 * t110, t101 * t108, t106, t105, qJD(2) ^ 2 / 0.2e1, pkin(1) * t108 - pkin(5) * t106, -t104 * pkin(1) * t101 - pkin(5) * t105, -t85 * t92 + t86 * t91, t86 ^ 2 / 0.2e1 + t85 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1, t84 ^ 2 / 0.2e1, -t84 * t83, t84 * t97, -t83 * t97, t97 ^ 2 / 0.2e1, t87 * t83 + (-t100 * t82 + t102 * t81) * t97, t87 * t84 - (t100 * t81 + t102 * t82) * t97;];
T_reg = t1;
