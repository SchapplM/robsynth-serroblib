% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:55
% EndTime: 2019-12-31 17:36:55
% DurationCPUTime: 0.05s
% Computational Cost: add. (60->29), mult. (151->56), div. (0->0), fcn. (80->4), ass. (0->21)
t93 = sin(pkin(8));
t94 = cos(pkin(8));
t98 = qJ(3) * qJD(2);
t88 = t93 * qJD(1) + t94 * t98;
t100 = qJD(2) * t93;
t99 = qJD(2) * t94;
t97 = qJ(4) * t93 + pkin(2);
t87 = t94 * qJD(1) - t93 * t98;
t86 = qJD(4) - t87;
t96 = cos(qJ(5));
t95 = sin(qJ(5));
t91 = -qJD(2) * pkin(2) + qJD(3);
t85 = t88 ^ 2 / 0.2e1;
t84 = (t93 * t96 - t94 * t95) * qJD(2);
t83 = (t93 * t95 + t94 * t96) * qJD(2);
t82 = t88 * t99;
t81 = qJD(3) + (-pkin(3) * t94 - t97) * qJD(2);
t80 = -pkin(6) * t99 + t88;
t79 = -pkin(6) * t100 + t86;
t78 = -qJD(3) + ((pkin(3) + pkin(4)) * t94 + t97) * qJD(2);
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, -t91 * t99, t91 * t100, -t87 * t100 + t82, t85 + t87 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1, -t81 * t99, t86 * t100 + t82, -t81 * t100, t85 + t81 ^ 2 / 0.2e1 + t86 ^ 2 / 0.2e1, t84 ^ 2 / 0.2e1, -t84 * t83, t84 * qJD(5), -t83 * qJD(5), qJD(5) ^ 2 / 0.2e1, t78 * t83 + (t96 * t79 - t95 * t80) * qJD(5), t78 * t84 - (t95 * t79 + t96 * t80) * qJD(5);];
T_reg = t1;
