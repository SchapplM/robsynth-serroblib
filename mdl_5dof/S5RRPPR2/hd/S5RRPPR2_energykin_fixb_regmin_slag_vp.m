% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:21
% EndTime: 2019-12-05 18:20:21
% DurationCPUTime: 0.06s
% Computational Cost: add. (128->26), mult. (208->62), div. (0->0), fcn. (108->8), ass. (0->27)
t101 = sin(pkin(9));
t100 = qJD(1) + qJD(2);
t98 = t100 ^ 2;
t117 = t98 * t101 ^ 2;
t102 = sin(pkin(8));
t104 = cos(pkin(8));
t116 = pkin(1) * qJD(1);
t113 = sin(qJ(2)) * t116;
t112 = cos(qJ(2)) * t116;
t92 = t100 * pkin(2) + t112;
t90 = t102 * t92 + t104 * t113;
t115 = t100 * t101;
t103 = cos(pkin(9));
t114 = t103 * t100;
t105 = sin(qJ(5));
t111 = t105 * t115;
t107 = cos(qJ(5));
t110 = t107 * t115;
t89 = -t102 * t113 + t104 * t92;
t109 = qJD(4) - t89;
t93 = -qJD(5) + t114;
t88 = t100 * qJ(4) + t90;
t87 = -t100 * pkin(3) + t109;
t86 = t101 * qJD(3) + t103 * t88;
t84 = -t103 * qJD(3) + t101 * t88;
t83 = (-pkin(4) * t103 - pkin(7) * t101 - pkin(3)) * t100 + t109;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t98 / 0.2e1, t100 * t112, -t100 * t113, t90 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, -t87 * t114, t87 * t115, (t101 * t84 + t103 * t86) * t100, t86 ^ 2 / 0.2e1 + t84 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1, t107 ^ 2 * t117 / 0.2e1, -t107 * t105 * t117, -t93 * t110, t93 * t111, t93 ^ 2 / 0.2e1, -(-t105 * t86 + t107 * t83) * t93 + t84 * t111, (t105 * t83 + t107 * t86) * t93 + t84 * t110;];
T_reg = t1;
