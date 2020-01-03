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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:57:35
% EndTime: 2020-01-03 11:57:35
% DurationCPUTime: 0.05s
% Computational Cost: add. (128->26), mult. (208->62), div. (0->0), fcn. (108->8), ass. (0->27)
t100 = sin(pkin(9));
t99 = qJD(1) + qJD(2);
t97 = t99 ^ 2;
t116 = t97 * t100 ^ 2;
t101 = sin(pkin(8));
t103 = cos(pkin(8));
t115 = pkin(1) * qJD(1);
t112 = sin(qJ(2)) * t115;
t109 = cos(qJ(2)) * t115;
t91 = t99 * pkin(2) + t109;
t89 = t101 * t91 + t103 * t112;
t114 = t100 * t99;
t102 = cos(pkin(9));
t113 = t102 * t99;
t104 = sin(qJ(5));
t111 = t104 * t114;
t106 = cos(qJ(5));
t110 = t106 * t114;
t88 = -t101 * t112 + t103 * t91;
t108 = qJD(4) - t88;
t92 = -qJD(5) + t113;
t87 = t99 * qJ(4) + t89;
t86 = -t99 * pkin(3) + t108;
t85 = t100 * qJD(3) + t102 * t87;
t83 = -t102 * qJD(3) + t100 * t87;
t82 = (-pkin(4) * t102 - pkin(7) * t100 - pkin(3)) * t99 + t108;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t97 / 0.2e1, t99 * t109, -t99 * t112, t89 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, -t86 * t113, t86 * t114, (t100 * t83 + t102 * t85) * t99, t85 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1 + t86 ^ 2 / 0.2e1, t106 ^ 2 * t116 / 0.2e1, -t106 * t104 * t116, -t92 * t110, t92 * t111, t92 ^ 2 / 0.2e1, -(-t104 * t85 + t106 * t82) * t92 + t83 * t111, (t104 * t82 + t106 * t85) * t92 + t83 * t110;];
T_reg = t1;
