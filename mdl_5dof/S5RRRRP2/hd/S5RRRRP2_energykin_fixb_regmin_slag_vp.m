% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP2
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
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:45
% EndTime: 2020-01-03 12:11:45
% DurationCPUTime: 0.10s
% Computational Cost: add. (154->27), mult. (230->65), div. (0->0), fcn. (125->6), ass. (0->29)
t98 = qJD(1) + qJD(2);
t96 = t98 ^ 2;
t115 = t96 / 0.2e1;
t114 = cos(qJ(4));
t100 = sin(qJ(3));
t112 = pkin(1) * qJD(1);
t107 = sin(qJ(2)) * t112;
t93 = pkin(7) * t98 + t107;
t105 = pkin(8) * t98 + t93;
t88 = qJD(3) * pkin(3) - t100 * t105;
t102 = cos(qJ(3));
t89 = t105 * t102;
t99 = sin(qJ(4));
t113 = t114 * t89 + t99 * t88;
t111 = t100 * t98;
t110 = t102 * t98;
t109 = qJD(3) * t100;
t108 = qJD(3) * t102;
t106 = cos(qJ(2)) * t112;
t104 = t114 * t88 - t89 * t99;
t92 = -t106 + (-pkin(3) * t102 - pkin(2)) * t98;
t97 = qJD(3) + qJD(4);
t94 = -pkin(2) * t98 - t106;
t91 = (t100 * t114 + t102 * t99) * t98;
t90 = -t110 * t114 + t111 * t99;
t84 = t90 * pkin(4) + qJD(5) + t92;
t83 = -qJ(5) * t90 + t113;
t82 = pkin(4) * t97 - qJ(5) * t91 + t104;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t115, t98 * t106, -t98 * t107, t100 ^ 2 * t115, t100 * t96 * t102, t98 * t109, t98 * t108, qJD(3) ^ 2 / 0.2e1, -t109 * t93 - t110 * t94, -t108 * t93 + t111 * t94, t91 ^ 2 / 0.2e1, -t91 * t90, t91 * t97, -t90 * t97, t97 ^ 2 / 0.2e1, t104 * t97 + t92 * t90, -t113 * t97 + t92 * t91, -t82 * t91 - t83 * t90, t83 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1 + t84 ^ 2 / 0.2e1;];
T_reg = t1;
