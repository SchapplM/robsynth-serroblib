% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:29
% EndTime: 2019-12-05 18:58:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (155->24), mult. (190->62), div. (0->0), fcn. (101->8), ass. (0->30)
t93 = qJD(1) + qJD(2);
t91 = qJD(3) + t93;
t90 = t91 ^ 2;
t112 = t90 / 0.2e1;
t95 = sin(qJ(4));
t111 = t91 * t95;
t99 = cos(qJ(4));
t110 = t91 * t99;
t100 = cos(qJ(3));
t108 = pkin(1) * qJD(1);
t105 = sin(qJ(2)) * t108;
t104 = cos(qJ(2)) * t108;
t86 = t93 * pkin(2) + t104;
t96 = sin(qJ(3));
t109 = t100 * t105 + t96 * t86;
t107 = qJD(4) * t95;
t106 = qJD(4) * t99;
t82 = t91 * pkin(8) + t109;
t103 = pkin(9) * t91 + t82;
t102 = t100 * t86 - t96 * t105;
t98 = cos(qJ(5));
t94 = sin(qJ(5));
t92 = qJD(4) + qJD(5);
t84 = (t94 * t99 + t95 * t98) * t91;
t83 = -t98 * t110 + t94 * t111;
t81 = -t91 * pkin(3) - t102;
t80 = (-pkin(4) * t99 - pkin(3)) * t91 - t102;
t79 = t103 * t99;
t78 = qJD(4) * pkin(4) - t103 * t95;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t93 ^ 2 / 0.2e1, t93 * t104, -t93 * t105, t112, t102 * t91, -t109 * t91, t95 ^ 2 * t112, t95 * t90 * t99, t91 * t107, t91 * t106, qJD(4) ^ 2 / 0.2e1, -t82 * t107 - t81 * t110, -t82 * t106 + t81 * t111, t84 ^ 2 / 0.2e1, -t84 * t83, t84 * t92, -t83 * t92, t92 ^ 2 / 0.2e1, t80 * t83 + (t98 * t78 - t94 * t79) * t92, t80 * t84 - (t94 * t78 + t98 * t79) * t92;];
T_reg = t1;
