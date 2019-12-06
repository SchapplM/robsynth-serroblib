% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:32
% EndTime: 2019-12-05 17:51:32
% DurationCPUTime: 0.06s
% Computational Cost: add. (113->26), mult. (212->64), div. (0->0), fcn. (108->8), ass. (0->28)
t101 = sin(pkin(9));
t100 = qJD(1) + qJD(3);
t98 = t100 ^ 2;
t119 = t98 * t101 ^ 2;
t106 = sin(qJ(3));
t108 = cos(qJ(3));
t102 = sin(pkin(8));
t115 = pkin(1) * qJD(1) * t102;
t104 = cos(pkin(8));
t92 = (pkin(1) * t104 + pkin(2)) * qJD(1);
t118 = t106 * t92 + t108 * t115;
t117 = t100 * t101;
t103 = cos(pkin(9));
t116 = t103 * t100;
t105 = sin(qJ(5));
t114 = t105 * t117;
t107 = cos(qJ(5));
t113 = t107 * t117;
t112 = -t106 * t115 + t108 * t92;
t111 = qJD(4) - t112;
t109 = qJD(1) ^ 2;
t93 = -qJD(5) + t116;
t90 = t100 * qJ(4) + t118;
t89 = -t100 * pkin(3) + t111;
t88 = t101 * qJD(2) + t103 * t90;
t86 = -t103 * qJD(2) + t101 * t90;
t85 = (-pkin(4) * t103 - pkin(7) * t101 - pkin(3)) * t100 + t111;
t1 = [t109 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t102 ^ 2 / 0.2e1 + t104 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t109, t98 / 0.2e1, t112 * t100, -t118 * t100, -t89 * t116, t89 * t117, (t101 * t86 + t103 * t88) * t100, t88 ^ 2 / 0.2e1 + t86 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1, t107 ^ 2 * t119 / 0.2e1, -t107 * t105 * t119, -t93 * t113, t93 * t114, t93 ^ 2 / 0.2e1, -(-t105 * t88 + t107 * t85) * t93 + t86 * t114, (t105 * t85 + t107 * t88) * t93 + t86 * t113;];
T_reg = t1;
