% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR14V3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energykin_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:10:05
% EndTime: 2019-04-12 15:10:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (136->30), mult. (342->78), div. (0->0), fcn. (253->8), ass. (0->30)
t110 = sin(qJ(4));
t114 = cos(qJ(4));
t111 = sin(qJ(2));
t123 = qJD(1) * t111;
t126 = t114 * qJD(2) - t110 * t123;
t116 = qJD(2) ^ 2;
t125 = t116 / 0.2e1;
t117 = qJD(1) ^ 2;
t124 = t111 ^ 2 * t117;
t122 = qJD(1) * qJD(2);
t115 = cos(qJ(2));
t121 = t111 * t117 * t115;
t119 = t124 / 0.2e1;
t102 = t110 * qJD(2) + t114 * t123;
t104 = t115 * qJD(1) - qJD(4);
t109 = sin(qJ(5));
t113 = cos(qJ(5));
t93 = t109 * t102 + t113 * t104;
t112 = cos(qJ(6));
t108 = sin(qJ(6));
t100 = qJD(5) - t126;
t99 = t126 * qJ(3);
t98 = t102 * qJ(3);
t96 = t109 * qJD(3) + t113 * t99;
t95 = -t113 * qJD(3) + t109 * t99;
t94 = t113 * t102 - t109 * t104;
t92 = qJD(6) + t93;
t91 = t108 * t100 + t112 * t94;
t90 = -t112 * t100 + t108 * t94;
t1 = [t117 / 0.2e1, 0, 0, t119, t121, t111 * t122, t115 * t122, t125, 0, 0, qJ(3) * t121 - qJD(3) * qJD(2) (qJ(3) * qJD(2) * t115 + qJD(3) * t111) * qJD(1) (t116 + t124) * qJ(3), qJD(3) ^ 2 / 0.2e1 + (t125 + t119) * qJ(3) ^ 2, t102 ^ 2 / 0.2e1, t102 * t126, -t102 * t104, -t126 * t104, t104 ^ 2 / 0.2e1, -qJD(3) * t126 + t98 * t104, qJD(3) * t102 + t99 * t104, t94 ^ 2 / 0.2e1, -t94 * t93, t94 * t100, -t93 * t100, t100 ^ 2 / 0.2e1, -t95 * t100 + t98 * t93, -t96 * t100 + t98 * t94, t91 ^ 2 / 0.2e1, -t91 * t90, t91 * t92, -t90 * t92, t92 ^ 2 / 0.2e1 (-t108 * t96 + t112 * t98) * t92 + t95 * t90 -(t108 * t98 + t112 * t96) * t92 + t95 * t91;];
T_reg  = t1;
