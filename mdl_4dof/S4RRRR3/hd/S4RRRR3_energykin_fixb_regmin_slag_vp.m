% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:36
% EndTime: 2019-12-31 17:24:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (95->27), mult. (261->66), div. (0->0), fcn. (175->6), ass. (0->30)
t99 = qJD(1) ^ 2;
t110 = t99 / 0.2e1;
t109 = -pkin(6) - pkin(5);
t108 = cos(qJ(3));
t98 = cos(qJ(2));
t107 = t98 * t99;
t96 = sin(qJ(2));
t105 = qJD(1) * t96;
t88 = qJD(2) * pkin(2) + t109 * t105;
t104 = qJD(1) * t98;
t89 = t109 * t104;
t95 = sin(qJ(3));
t106 = -t108 * t89 + t95 * t88;
t103 = qJD(1) * qJD(2);
t93 = qJD(2) + qJD(3);
t102 = t96 * t103;
t101 = t98 * t103;
t100 = t108 * t88 + t95 * t89;
t90 = (-pkin(2) * t98 - pkin(1)) * qJD(1);
t97 = cos(qJ(4));
t94 = sin(qJ(4));
t92 = qJD(4) + t93;
t86 = (t108 * t96 + t95 * t98) * qJD(1);
t85 = -t108 * t104 + t95 * t105;
t81 = t85 * pkin(3) + t90;
t80 = -t94 * t85 + t97 * t86;
t79 = t97 * t85 + t94 * t86;
t78 = -t85 * pkin(7) + t106;
t77 = t93 * pkin(3) - t86 * pkin(7) + t100;
t1 = [t110, 0, 0, t96 ^ 2 * t110, t96 * t107, t102, t101, qJD(2) ^ 2 / 0.2e1, pkin(1) * t107 - pkin(5) * t102, -t99 * pkin(1) * t96 - pkin(5) * t101, t86 ^ 2 / 0.2e1, -t86 * t85, t86 * t93, -t85 * t93, t93 ^ 2 / 0.2e1, t100 * t93 + t90 * t85, -t106 * t93 + t90 * t86, t80 ^ 2 / 0.2e1, -t80 * t79, t80 * t92, -t79 * t92, t92 ^ 2 / 0.2e1, t81 * t79 + (t97 * t77 - t94 * t78) * t92, t81 * t80 - (t94 * t77 + t97 * t78) * t92;];
T_reg = t1;
