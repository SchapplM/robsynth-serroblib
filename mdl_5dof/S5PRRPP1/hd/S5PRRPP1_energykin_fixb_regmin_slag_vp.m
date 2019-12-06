% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:51
% EndTime: 2019-12-05 16:06:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (107->27), mult. (253->62), div. (0->0), fcn. (143->4), ass. (0->23)
t101 = qJD(2) ^ 2;
t107 = t101 / 0.2e1;
t99 = sin(qJ(3));
t105 = qJD(2) * t99;
t100 = cos(qJ(3));
t96 = t100 * qJD(1);
t88 = qJD(3) * pkin(3) + t96 + (-pkin(6) - qJ(4)) * t105;
t103 = qJD(2) * t100;
t106 = pkin(6) * t103 + t99 * qJD(1);
t89 = qJ(4) * t103 + t106;
t97 = sin(pkin(8));
t98 = cos(pkin(8));
t85 = t97 * t88 + t98 * t89;
t104 = t100 * t101;
t102 = qJD(2) * qJD(3);
t84 = t98 * t88 - t97 * t89;
t92 = qJD(4) + (-pkin(3) * t100 - pkin(2)) * qJD(2);
t91 = (t100 * t97 + t98 * t99) * qJD(2);
t90 = -t98 * t103 + t97 * t105;
t83 = t90 * pkin(4) - t91 * qJ(5) + t92;
t82 = qJD(3) * qJ(5) + t85;
t81 = -qJD(3) * pkin(4) + qJD(5) - t84;
t1 = [qJD(1) ^ 2 / 0.2e1, t107, 0, 0, t99 ^ 2 * t107, t99 * t104, t99 * t102, t100 * t102, qJD(3) ^ 2 / 0.2e1, pkin(2) * t104 + (-pkin(6) * t105 + t96) * qJD(3), -t101 * pkin(2) * t99 - t106 * qJD(3), -t84 * t91 - t85 * t90, t85 ^ 2 / 0.2e1 + t84 ^ 2 / 0.2e1 + t92 ^ 2 / 0.2e1, -t81 * qJD(3) + t83 * t90, t81 * t91 - t82 * t90, t82 * qJD(3) - t83 * t91, t82 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1 + t81 ^ 2 / 0.2e1;];
T_reg = t1;
