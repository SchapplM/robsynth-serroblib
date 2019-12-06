% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR1
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:26
% EndTime: 2019-12-05 18:18:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (135->29), mult. (219->64), div. (0->0), fcn. (120->8), ass. (0->27)
t94 = qJD(1) + qJD(2);
t95 = sin(pkin(9));
t108 = t94 * t95;
t97 = cos(pkin(9));
t107 = t94 * t97;
t106 = pkin(1) * qJD(1);
t105 = sin(qJ(2)) * t106;
t104 = cos(qJ(2)) * t106;
t88 = t94 * pkin(2) + t104;
t96 = sin(pkin(8));
t98 = cos(pkin(8));
t84 = t98 * t105 + t96 * t88;
t82 = t94 * qJ(4) + t84;
t78 = t95 * qJD(3) + t97 * t82;
t83 = -t96 * t105 + t98 * t88;
t103 = qJD(4) - t83;
t101 = cos(qJ(5));
t99 = sin(qJ(5));
t93 = t97 * qJD(3);
t86 = (t101 * t95 + t97 * t99) * t94;
t85 = -t101 * t107 + t99 * t108;
t81 = -t94 * pkin(3) + t103;
t79 = (-pkin(4) * t97 - pkin(3)) * t94 + t103;
t77 = -t95 * t82 + t93;
t76 = pkin(7) * t107 + t78;
t75 = t93 + (-pkin(7) * t94 - t82) * t95;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t94 ^ 2 / 0.2e1, t94 * t104, -t94 * t105, t84 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, -t81 * t107, t81 * t108, (-t77 * t95 + t78 * t97) * t94, t78 ^ 2 / 0.2e1 + t77 ^ 2 / 0.2e1 + t81 ^ 2 / 0.2e1, t86 ^ 2 / 0.2e1, -t86 * t85, t86 * qJD(5), -t85 * qJD(5), qJD(5) ^ 2 / 0.2e1, t79 * t85 + (t101 * t75 - t99 * t76) * qJD(5), t79 * t86 - (t101 * t76 + t99 * t75) * qJD(5);];
T_reg = t1;
