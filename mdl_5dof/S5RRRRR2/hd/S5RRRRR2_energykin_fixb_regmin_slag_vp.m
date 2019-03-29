% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:30
% EndTime: 2019-03-29 15:26:30
% DurationCPUTime: 0.06s
% Computational Cost: add. (119->24), mult. (225->71), div. (0->0), fcn. (153->8), ass. (0->30)
t93 = qJD(1) + qJD(2);
t91 = t93 ^ 2;
t110 = t91 / 0.2e1;
t109 = pkin(1) * qJD(1);
t100 = cos(qJ(3));
t108 = t100 * t93;
t101 = cos(qJ(2));
t107 = t101 * t93;
t96 = sin(qJ(3));
t106 = qJD(3) * t96;
t105 = qJD(3) * t100;
t97 = sin(qJ(2));
t104 = t97 * t109;
t103 = t101 * t109;
t102 = t100 * t104;
t95 = sin(qJ(4));
t99 = cos(qJ(4));
t85 = t95 * t96 * t93 - t99 * t108;
t98 = cos(qJ(5));
t94 = sin(qJ(5));
t92 = qJD(3) + qJD(4);
t88 = qJD(3) * pkin(2) - t96 * t104;
t87 = -pkin(2) * t108 - t103;
t86 = (t100 * t95 + t96 * t99) * t93;
t84 = qJD(5) + t85;
t83 = t99 * t102 + t95 * t88;
t82 = t95 * t102 - t99 * t88;
t81 = t98 * t86 + t94 * t92;
t80 = t94 * t86 - t98 * t92;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t110, t93 * t103, -t93 * t104, t96 ^ 2 * t110, t96 * t91 * t100, t93 * t106, t93 * t105, qJD(3) ^ 2 / 0.2e1 (t100 * t107 - t97 * t106) * t109 (-t97 * t105 - t96 * t107) * t109, t86 ^ 2 / 0.2e1, -t86 * t85, t86 * t92, -t85 * t92, t92 ^ 2 / 0.2e1, -t82 * t92 + t87 * t85, -t83 * t92 + t87 * t86, t81 ^ 2 / 0.2e1, -t81 * t80, t81 * t84, -t80 * t84, t84 ^ 2 / 0.2e1 (-t94 * t83 + t98 * t87) * t84 + t82 * t80 -(t98 * t83 + t94 * t87) * t84 + t82 * t81;];
T_reg  = t1;
