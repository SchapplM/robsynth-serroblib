% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:56
% EndTime: 2019-12-31 18:40:57
% DurationCPUTime: 0.10s
% Computational Cost: add. (108->24), mult. (207->61), div. (0->0), fcn. (95->6), ass. (0->25)
t89 = qJD(1) + qJD(3);
t88 = t89 ^ 2;
t106 = t88 / 0.2e1;
t92 = sin(qJ(4));
t105 = t89 * t92;
t94 = cos(qJ(4));
t104 = t89 * t94;
t90 = sin(pkin(8));
t100 = pkin(1) * qJD(1) * t90;
t91 = cos(pkin(8));
t84 = (pkin(1) * t91 + pkin(2)) * qJD(1);
t93 = sin(qJ(3));
t95 = cos(qJ(3));
t102 = t95 * t100 + t93 * t84;
t82 = pkin(7) * t89 + t102;
t103 = t92 * qJD(2) + t94 * t82;
t101 = qJD(4) * t89;
t99 = -t93 * t100 + t95 * t84;
t98 = qJD(2) * t94 - t82 * t92;
t96 = qJD(1) ^ 2;
t81 = -pkin(3) * t89 - t99;
t79 = qJD(4) * qJ(5) + t103;
t78 = -qJD(4) * pkin(4) + qJD(5) - t98;
t77 = (-pkin(4) * t94 - qJ(5) * t92 - pkin(3)) * t89 - t99;
t1 = [t96 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t90 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t96, t106, t99 * t89, -t102 * t89, t92 ^ 2 * t106, t92 * t88 * t94, t92 * t101, t94 * t101, qJD(4) ^ 2 / 0.2e1, qJD(4) * t98 - t104 * t81, -qJD(4) * t103 + t105 * t81, -qJD(4) * t78 - t104 * t77, (t78 * t92 + t79 * t94) * t89, qJD(4) * t79 - t105 * t77, t79 ^ 2 / 0.2e1 + t77 ^ 2 / 0.2e1 + t78 ^ 2 / 0.2e1;];
T_reg = t1;
