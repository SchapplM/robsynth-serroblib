% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:25
% EndTime: 2019-03-08 18:26:25
% DurationCPUTime: 0.10s
% Computational Cost: add. (87->34), mult. (287->72), div. (0->0), fcn. (168->4), ass. (0->24)
t94 = sin(pkin(6));
t95 = sin(pkin(4));
t106 = t94 * t95;
t96 = cos(pkin(6));
t105 = t95 * t96;
t97 = cos(pkin(4));
t103 = qJD(1) * t97;
t101 = pkin(1) * t103;
t104 = qJD(1) * t95;
t99 = qJ(2) * t104;
t89 = t94 * t101 + t96 * t99;
t90 = t94 * t99;
t102 = qJD(3) + t90;
t100 = -pkin(1) * t96 - pkin(2);
t98 = -qJ(3) * t94 - pkin(1);
t93 = -pkin(1) * t104 + qJD(2);
t88 = t96 * t101 - t90;
t87 = -qJ(3) * t103 - t89;
t86 = qJD(2) + (-pkin(2) * t96 + t98) * t104;
t85 = t100 * t103 + t102;
t84 = qJD(2) + ((-pkin(2) - qJ(4)) * t96 + t98) * t104;
t83 = qJD(4) + (pkin(3) * t105 + qJ(3) * t97) * qJD(1) + t89;
t82 = (pkin(3) * t106 + (-qJ(4) + t100) * t97) * qJD(1) + t102;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (-t93 * t105 + t88 * t97) * qJD(1) (t93 * t106 - t89 * t97) * qJD(1) (-t88 * t94 + t89 * t96) * t104, t89 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1 + t93 ^ 2 / 0.2e1 (t85 * t94 - t87 * t96) * t104 (t86 * t105 + t85 * t97) * qJD(1) (-t86 * t106 - t87 * t97) * qJD(1), t86 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1 + t85 ^ 2 / 0.2e1 (t82 * t94 + t83 * t96) * t104 (-t84 * t106 + t83 * t97) * qJD(1) (-t84 * t105 - t82 * t97) * qJD(1), t84 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1;];
T_reg  = t1;
