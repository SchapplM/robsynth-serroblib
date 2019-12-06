% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:17
% EndTime: 2019-12-05 16:49:17
% DurationCPUTime: 0.07s
% Computational Cost: add. (83->26), mult. (209->65), div. (0->0), fcn. (125->6), ass. (0->29)
t108 = qJD(2) ^ 2;
t119 = t108 / 0.2e1;
t118 = cos(qJ(4));
t103 = sin(qJ(4));
t104 = sin(qJ(3));
t105 = sin(qJ(2));
t99 = qJD(2) * pkin(6) + t105 * qJD(1);
t110 = pkin(7) * qJD(2) + t99;
t94 = qJD(3) * pkin(3) - t110 * t104;
t106 = cos(qJ(3));
t95 = t110 * t106;
t117 = t103 * t94 + t118 * t95;
t116 = qJD(3) * t99;
t115 = qJD(2) * t104;
t114 = qJD(2) * t106;
t107 = cos(qJ(2));
t113 = t107 * qJD(1);
t112 = qJD(1) * qJD(2);
t111 = qJD(2) * qJD(3);
t109 = -t103 * t95 + t118 * t94;
t98 = -t113 + (-pkin(3) * t106 - pkin(2)) * qJD(2);
t102 = qJD(3) + qJD(4);
t100 = -qJD(2) * pkin(2) - t113;
t97 = (t103 * t106 + t118 * t104) * qJD(2);
t96 = t103 * t115 - t118 * t114;
t90 = t96 * pkin(4) + qJD(5) + t98;
t89 = -t96 * qJ(5) + t117;
t88 = t102 * pkin(4) - t97 * qJ(5) + t109;
t1 = [qJD(1) ^ 2 / 0.2e1, t119, t107 * t112, -t105 * t112, t104 ^ 2 * t119, t104 * t108 * t106, t104 * t111, t106 * t111, qJD(3) ^ 2 / 0.2e1, -t100 * t114 - t104 * t116, t100 * t115 - t106 * t116, t97 ^ 2 / 0.2e1, -t97 * t96, t97 * t102, -t96 * t102, t102 ^ 2 / 0.2e1, t109 * t102 + t98 * t96, -t117 * t102 + t98 * t97, -t88 * t97 - t89 * t96, t89 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1 + t90 ^ 2 / 0.2e1;];
T_reg = t1;
