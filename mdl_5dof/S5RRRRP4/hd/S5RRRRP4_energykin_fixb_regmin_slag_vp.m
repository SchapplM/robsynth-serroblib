% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:12
% EndTime: 2019-12-31 21:51:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (206->29), mult. (291->69), div. (0->0), fcn. (157->6), ass. (0->29)
t100 = qJD(1) + qJD(2);
t98 = t100 ^ 2;
t117 = t98 / 0.2e1;
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t102 = sin(qJ(3));
t115 = pkin(1) * qJD(1);
t110 = sin(qJ(2)) * t115;
t95 = t100 * pkin(7) + t110;
t108 = pkin(8) * t100 + t95;
t90 = qJD(3) * pkin(3) - t108 * t102;
t105 = cos(qJ(3));
t91 = t108 * t105;
t116 = t101 * t90 + t104 * t91;
t114 = t100 * t102;
t113 = t100 * t105;
t112 = qJD(3) * t102;
t111 = qJD(3) * t105;
t109 = cos(qJ(2)) * t115;
t107 = -t101 * t91 + t104 * t90;
t94 = -t109 + (-pkin(3) * t105 - pkin(2)) * t100;
t99 = qJD(3) + qJD(4);
t96 = -t100 * pkin(2) - t109;
t93 = (t101 * t105 + t102 * t104) * t100;
t92 = t101 * t114 - t104 * t113;
t87 = t92 * pkin(4) - t93 * qJ(5) + t94;
t86 = t99 * qJ(5) + t116;
t85 = -t99 * pkin(4) + qJD(5) - t107;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t117, t100 * t109, -t100 * t110, t102 ^ 2 * t117, t102 * t98 * t105, t100 * t112, t100 * t111, qJD(3) ^ 2 / 0.2e1, -t95 * t112 - t96 * t113, -t95 * t111 + t96 * t114, t93 ^ 2 / 0.2e1, -t93 * t92, t93 * t99, -t92 * t99, t99 ^ 2 / 0.2e1, t107 * t99 + t94 * t92, -t116 * t99 + t94 * t93, -t85 * t99 + t87 * t92, t85 * t93 - t86 * t92, t86 * t99 - t87 * t93, t86 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1 + t85 ^ 2 / 0.2e1;];
T_reg = t1;
