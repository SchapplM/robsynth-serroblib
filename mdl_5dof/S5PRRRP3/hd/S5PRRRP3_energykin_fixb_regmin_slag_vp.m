% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP3
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:12
% EndTime: 2019-12-05 16:44:12
% DurationCPUTime: 0.05s
% Computational Cost: add. (77->26), mult. (200->61), div. (0->0), fcn. (118->4), ass. (0->24)
t92 = qJD(2) ^ 2;
t101 = t92 / 0.2e1;
t100 = cos(qJ(4));
t91 = cos(qJ(3));
t99 = t91 * t92;
t87 = t91 * qJD(1);
t90 = sin(qJ(3));
t96 = qJD(2) * t90;
t79 = qJD(3) * pkin(3) + t87 + (-pkin(7) - pkin(6)) * t96;
t95 = qJD(2) * t91;
t97 = pkin(6) * t95 + t90 * qJD(1);
t80 = pkin(7) * t95 + t97;
t89 = sin(qJ(4));
t98 = t100 * t80 + t89 * t79;
t94 = qJD(2) * qJD(3);
t93 = t100 * t79 - t89 * t80;
t83 = (-pkin(3) * t91 - pkin(2)) * qJD(2);
t88 = qJD(3) + qJD(4);
t82 = (t100 * t90 + t89 * t91) * qJD(2);
t81 = -t100 * t95 + t89 * t96;
t75 = t81 * pkin(4) + qJD(5) + t83;
t74 = -t81 * qJ(5) + t98;
t73 = t88 * pkin(4) - t82 * qJ(5) + t93;
t1 = [qJD(1) ^ 2 / 0.2e1, t101, 0, 0, t90 ^ 2 * t101, t90 * t99, t90 * t94, t91 * t94, qJD(3) ^ 2 / 0.2e1, pkin(2) * t99 + (-pkin(6) * t96 + t87) * qJD(3), -t92 * pkin(2) * t90 - t97 * qJD(3), t82 ^ 2 / 0.2e1, -t82 * t81, t82 * t88, -t81 * t88, t88 ^ 2 / 0.2e1, t83 * t81 + t93 * t88, t83 * t82 - t98 * t88, -t73 * t82 - t74 * t81, t74 ^ 2 / 0.2e1 + t73 ^ 2 / 0.2e1 + t75 ^ 2 / 0.2e1;];
T_reg = t1;
