% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:36
% EndTime: 2019-12-31 18:24:37
% DurationCPUTime: 0.11s
% Computational Cost: add. (89->32), mult. (217->76), div. (0->0), fcn. (102->6), ass. (0->29)
t125 = -pkin(3) - pkin(7);
t114 = qJD(1) ^ 2;
t124 = t114 / 0.2e1;
t108 = sin(pkin(8));
t104 = (pkin(1) * t108 + pkin(6)) * qJD(1);
t111 = sin(qJ(3));
t113 = cos(qJ(3));
t123 = t111 * qJD(2) + t113 * t104;
t122 = qJD(1) * t113;
t121 = t111 * qJD(1);
t120 = qJD(1) * qJD(3);
t109 = cos(pkin(8));
t119 = -pkin(1) * t109 - pkin(2);
t118 = t113 * qJD(2) - t111 * t104;
t98 = -qJD(3) * qJ(4) - t123;
t117 = qJD(4) - t118;
t116 = -qJ(4) * t111 + t119;
t112 = cos(qJ(5));
t110 = sin(qJ(5));
t106 = qJD(5) + t121;
t105 = t119 * qJD(1);
t103 = t112 * qJD(3) - t110 * t122;
t102 = t110 * qJD(3) + t112 * t122;
t99 = (-pkin(3) * t113 + t116) * qJD(1);
t97 = -qJD(3) * pkin(3) + t117;
t96 = (t125 * t113 + t116) * qJD(1);
t95 = pkin(4) * t122 - t98;
t94 = pkin(4) * t121 + t125 * qJD(3) + t117;
t1 = [t124, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t108 ^ 2 / 0.2e1 + t109 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t114, t111 ^ 2 * t124, t111 * t114 * t113, t111 * t120, t113 * t120, qJD(3) ^ 2 / 0.2e1, t118 * qJD(3) - t105 * t122, -t123 * qJD(3) + t105 * t121, (t111 * t97 - t113 * t98) * qJD(1), t97 * qJD(3) + t99 * t122, -t98 * qJD(3) - t99 * t121, t99 ^ 2 / 0.2e1 + t98 ^ 2 / 0.2e1 + t97 ^ 2 / 0.2e1, t103 ^ 2 / 0.2e1, -t103 * t102, t103 * t106, -t102 * t106, t106 ^ 2 / 0.2e1, (-t110 * t96 + t112 * t94) * t106 + t95 * t102, -(t110 * t94 + t112 * t96) * t106 + t95 * t103;];
T_reg = t1;
