% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP2
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
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:28
% EndTime: 2020-01-03 11:45:28
% DurationCPUTime: 0.06s
% Computational Cost: add. (78->23), mult. (159->57), div. (0->0), fcn. (72->6), ass. (0->24)
t85 = qJD(1) + qJD(3);
t84 = t85 ^ 2;
t101 = t84 / 0.2e1;
t87 = cos(pkin(8));
t79 = (pkin(1) * t87 + pkin(2)) * qJD(1);
t89 = sin(qJ(3));
t91 = cos(qJ(3));
t86 = sin(pkin(8));
t95 = pkin(1) * qJD(1) * t86;
t94 = t91 * t79 - t89 * t95;
t100 = (-t85 * pkin(3) - t94) * t85;
t98 = t89 * t79 + t91 * t95;
t77 = t85 * pkin(7) + t98;
t88 = sin(qJ(4));
t90 = cos(qJ(4));
t99 = t88 * qJD(2) + t90 * t77;
t97 = qJ(5) * t85;
t96 = qJD(4) * t85;
t92 = qJD(1) ^ 2;
t83 = t90 * qJD(2);
t74 = qJD(5) + (-pkin(4) * t90 - pkin(3)) * t85 - t94;
t73 = t90 * t97 + t99;
t72 = qJD(4) * pkin(4) + t83 + (-t77 - t97) * t88;
t1 = [t92 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t86 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t92, t101, t94 * t85, -t98 * t85, t88 ^ 2 * t101, t88 * t84 * t90, t88 * t96, t90 * t96, qJD(4) ^ 2 / 0.2e1, -t90 * t100 + (-t88 * t77 + t83) * qJD(4), -t99 * qJD(4) + t88 * t100, (-t72 * t88 + t73 * t90) * t85, t73 ^ 2 / 0.2e1 + t72 ^ 2 / 0.2e1 + t74 ^ 2 / 0.2e1;];
T_reg = t1;
