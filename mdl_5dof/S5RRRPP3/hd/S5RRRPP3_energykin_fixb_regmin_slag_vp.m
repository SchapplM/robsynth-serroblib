% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP3
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
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:40
% EndTime: 2019-12-31 20:53:40
% DurationCPUTime: 0.10s
% Computational Cost: add. (140->30), mult. (198->66), div. (0->0), fcn. (74->4), ass. (0->25)
t88 = qJD(1) + qJD(2);
t87 = t88 ^ 2;
t104 = t87 / 0.2e1;
t89 = sin(qJ(3));
t103 = t88 * t89;
t91 = cos(qJ(3));
t102 = t88 * t91;
t101 = -pkin(3) - qJ(5);
t100 = pkin(1) * qJD(1);
t99 = qJD(3) * t89;
t98 = qJD(3) * t91;
t95 = sin(qJ(2)) * t100;
t85 = t88 * pkin(7) + t95;
t97 = t89 * t85 + qJD(4);
t96 = qJD(3) * qJ(4);
t94 = cos(qJ(2)) * t100;
t93 = -qJ(4) * t89 - pkin(2);
t86 = -t88 * pkin(2) - t94;
t83 = -t91 * t85 - t96;
t82 = -qJD(3) * pkin(3) + t97;
t81 = -t94 + (-pkin(3) * t91 + t93) * t88;
t80 = t96 + qJD(5) + (pkin(4) * t88 + t85) * t91;
t79 = pkin(4) * t103 + t101 * qJD(3) + t97;
t78 = -t94 + (t101 * t91 + t93) * t88;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t104, t88 * t94, -t88 * t95, t89 ^ 2 * t104, t89 * t87 * t91, t88 * t99, t88 * t98, qJD(3) ^ 2 / 0.2e1, -t86 * t102 - t85 * t99, t86 * t103 - t85 * t98, (t82 * t89 - t83 * t91) * t88, t82 * qJD(3) + t81 * t102, -t83 * qJD(3) - t81 * t103, t81 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1, (t79 * t89 + t80 * t91) * t88, t80 * qJD(3) - t78 * t103, -t79 * qJD(3) - t78 * t102, t78 ^ 2 / 0.2e1 + t79 ^ 2 / 0.2e1 + t80 ^ 2 / 0.2e1;];
T_reg = t1;
