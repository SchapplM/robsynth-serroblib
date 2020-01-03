% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:58
% EndTime: 2019-12-31 20:51:58
% DurationCPUTime: 0.11s
% Computational Cost: add. (140->29), mult. (198->66), div. (0->0), fcn. (74->4), ass. (0->25)
t90 = qJD(1) + qJD(2);
t89 = t90 ^ 2;
t107 = t89 / 0.2e1;
t106 = pkin(3) + pkin(4);
t92 = sin(qJ(3));
t105 = t90 * t92;
t94 = cos(qJ(3));
t104 = t90 * t94;
t103 = pkin(1) * qJD(1);
t98 = sin(qJ(2)) * t103;
t87 = t90 * pkin(7) + t98;
t84 = qJD(3) * qJ(4) + t94 * t87;
t102 = qJ(5) * t90;
t101 = qJD(3) * t92;
t100 = qJD(3) * t94;
t99 = t92 * t87 + qJD(4);
t97 = cos(qJ(2)) * t103;
t96 = qJ(4) * t92 + pkin(2);
t88 = -t90 * pkin(2) - t97;
t83 = -qJD(3) * pkin(3) + t99;
t82 = -t94 * t102 + t84;
t81 = -t97 + (-pkin(3) * t94 - t96) * t90;
t80 = -t106 * qJD(3) - t92 * t102 + t99;
t79 = t97 + qJD(5) + (t106 * t94 + t96) * t90;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t107, t90 * t97, -t90 * t98, t92 ^ 2 * t107, t92 * t89 * t94, t90 * t101, t90 * t100, qJD(3) ^ 2 / 0.2e1, -t87 * t101 - t88 * t104, -t87 * t100 + t88 * t105, -t83 * qJD(3) - t81 * t104, (t83 * t92 + t84 * t94) * t90, t84 * qJD(3) - t81 * t105, t84 ^ 2 / 0.2e1 + t81 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1, -t80 * qJD(3) + t79 * t104, t82 * qJD(3) + t79 * t105, (-t80 * t92 - t82 * t94) * t90, t82 ^ 2 / 0.2e1 + t80 ^ 2 / 0.2e1 + t79 ^ 2 / 0.2e1;];
T_reg = t1;
