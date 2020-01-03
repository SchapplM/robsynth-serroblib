% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:42
% EndTime: 2019-12-31 19:49:42
% DurationCPUTime: 0.05s
% Computational Cost: add. (123->24), mult. (203->59), div. (0->0), fcn. (95->6), ass. (0->24)
t90 = qJD(1) + qJD(2);
t89 = t90 ^ 2;
t105 = t89 / 0.2e1;
t93 = sin(qJ(4));
t104 = t90 * t93;
t95 = cos(qJ(4));
t103 = t90 * t95;
t101 = pkin(1) * qJD(1);
t98 = cos(qJ(2)) * t101;
t85 = t90 * pkin(2) + t98;
t91 = sin(pkin(8));
t92 = cos(pkin(8));
t99 = sin(qJ(2)) * t101;
t83 = t91 * t85 + t92 * t99;
t81 = t90 * pkin(7) + t83;
t102 = t93 * qJD(3) + t95 * t81;
t100 = qJD(4) * t90;
t82 = t92 * t85 - t91 * t99;
t97 = t95 * qJD(3) - t93 * t81;
t80 = -t90 * pkin(3) - t82;
t78 = qJD(4) * qJ(5) + t102;
t77 = -qJD(4) * pkin(4) + qJD(5) - t97;
t76 = (-pkin(4) * t95 - qJ(5) * t93 - pkin(3)) * t90 - t82;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t105, t90 * t98, -t90 * t99, t83 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t93 ^ 2 * t105, t93 * t89 * t95, t93 * t100, t95 * t100, qJD(4) ^ 2 / 0.2e1, t97 * qJD(4) - t80 * t103, -t102 * qJD(4) + t80 * t104, -t77 * qJD(4) - t76 * t103, (t77 * t93 + t78 * t95) * t90, t78 * qJD(4) - t76 * t104, t78 ^ 2 / 0.2e1 + t76 ^ 2 / 0.2e1 + t77 ^ 2 / 0.2e1;];
T_reg = t1;
