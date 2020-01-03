% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR2
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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:00
% EndTime: 2020-01-03 11:34:00
% DurationCPUTime: 0.06s
% Computational Cost: add. (120->29), mult. (223->66), div. (0->0), fcn. (120->8), ass. (0->28)
t94 = qJD(1) + qJD(3);
t95 = sin(pkin(9));
t110 = t94 * t95;
t97 = cos(pkin(9));
t109 = t94 * t97;
t100 = sin(qJ(3));
t102 = cos(qJ(3));
t96 = sin(pkin(8));
t107 = pkin(1) * qJD(1) * t96;
t98 = cos(pkin(8));
t89 = (pkin(1) * t98 + pkin(2)) * qJD(1);
t108 = t100 * t89 + t102 * t107;
t84 = t94 * qJ(4) + t108;
t80 = t95 * qJD(2) + t97 * t84;
t106 = -t100 * t107 + t102 * t89;
t105 = qJD(4) - t106;
t103 = qJD(1) ^ 2;
t101 = cos(qJ(5));
t99 = sin(qJ(5));
t93 = t97 * qJD(2);
t86 = (t101 * t95 + t97 * t99) * t94;
t85 = -t101 * t109 + t99 * t110;
t83 = -t94 * pkin(3) + t105;
t81 = (-pkin(4) * t97 - pkin(3)) * t94 + t105;
t79 = -t95 * t84 + t93;
t78 = pkin(7) * t109 + t80;
t77 = t93 + (-pkin(7) * t94 - t84) * t95;
t1 = [t103 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t96 ^ 2 / 0.2e1 + t98 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t103, t94 ^ 2 / 0.2e1, t106 * t94, -t108 * t94, -t83 * t109, t83 * t110, (-t79 * t95 + t80 * t97) * t94, t80 ^ 2 / 0.2e1 + t79 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1, t86 ^ 2 / 0.2e1, -t86 * t85, t86 * qJD(5), -t85 * qJD(5), qJD(5) ^ 2 / 0.2e1, t81 * t85 + (t101 * t77 - t99 * t78) * qJD(5), t81 * t86 - (t101 * t78 + t99 * t77) * qJD(5);];
T_reg = t1;
