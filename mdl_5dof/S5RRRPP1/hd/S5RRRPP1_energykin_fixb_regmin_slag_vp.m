% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:57
% EndTime: 2019-12-31 20:49:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (199->28), mult. (293->66), div. (0->0), fcn. (151->6), ass. (0->28)
t107 = qJD(1) + qJD(2);
t106 = t107 ^ 2;
t122 = t106 / 0.2e1;
t108 = sin(pkin(8));
t109 = cos(pkin(8));
t110 = sin(qJ(3));
t121 = pkin(1) * qJD(1);
t116 = sin(qJ(2)) * t121;
t103 = t107 * pkin(7) + t116;
t114 = qJ(4) * t107 + t103;
t98 = qJD(3) * pkin(3) - t114 * t110;
t112 = cos(qJ(3));
t99 = t114 * t112;
t95 = t108 * t98 + t109 * t99;
t120 = t107 * t110;
t119 = t107 * t112;
t118 = qJD(3) * t110;
t117 = qJD(3) * t112;
t115 = cos(qJ(2)) * t121;
t94 = -t108 * t99 + t109 * t98;
t100 = -t115 + qJD(4) + (-pkin(3) * t112 - pkin(2)) * t107;
t104 = -t107 * pkin(2) - t115;
t102 = (t108 * t112 + t109 * t110) * t107;
t101 = t108 * t120 - t109 * t119;
t93 = qJD(3) * qJ(5) + t95;
t92 = -qJD(3) * pkin(4) + qJD(5) - t94;
t91 = t101 * pkin(4) - t102 * qJ(5) + t100;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t122, t107 * t115, -t107 * t116, t110 ^ 2 * t122, t110 * t106 * t112, t107 * t118, t107 * t117, qJD(3) ^ 2 / 0.2e1, -t103 * t118 - t104 * t119, -t103 * t117 + t104 * t120, -t95 * t101 - t94 * t102, t95 ^ 2 / 0.2e1 + t94 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1, -t92 * qJD(3) + t91 * t101, -t93 * t101 + t92 * t102, t93 * qJD(3) - t91 * t102, t93 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1 + t92 ^ 2 / 0.2e1;];
T_reg = t1;
