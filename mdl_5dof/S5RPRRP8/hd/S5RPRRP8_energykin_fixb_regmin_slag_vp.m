% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:29
% EndTime: 2019-12-31 18:47:29
% DurationCPUTime: 0.05s
% Computational Cost: add. (117->24), mult. (166->57), div. (0->0), fcn. (55->4), ass. (0->24)
t85 = -qJD(1) + qJD(3);
t84 = t85 ^ 2;
t99 = t84 / 0.2e1;
t90 = qJD(1) ^ 2;
t98 = t90 / 0.2e1;
t86 = sin(qJ(4));
t97 = t85 * t86;
t88 = cos(qJ(4));
t96 = t85 * t88;
t80 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t87 = sin(qJ(3));
t89 = cos(qJ(3));
t92 = qJ(2) * qJD(1);
t95 = t87 * t80 + t89 * t92;
t94 = qJD(4) * t86;
t93 = qJD(4) * t88;
t91 = t89 * t80 - t87 * t92;
t83 = -qJD(1) * pkin(1) + qJD(2);
t78 = t85 * pkin(7) + t95;
t77 = -t85 * pkin(3) - t91;
t76 = qJD(4) * qJ(5) + t88 * t78;
t75 = -qJD(4) * pkin(4) + t86 * t78 + qJD(5);
t74 = (-pkin(4) * t88 - qJ(5) * t86 - pkin(3)) * t85 - t91;
t1 = [t98, 0, 0, -t83 * qJD(1), t90 * qJ(2), qJ(2) ^ 2 * t98 + t83 ^ 2 / 0.2e1, t99, t91 * t85, -t95 * t85, t86 ^ 2 * t99, t86 * t84 * t88, t85 * t94, t85 * t93, qJD(4) ^ 2 / 0.2e1, -t77 * t96 - t78 * t94, t77 * t97 - t78 * t93, -t75 * qJD(4) - t74 * t96, (t75 * t86 + t76 * t88) * t85, t76 * qJD(4) - t74 * t97, t76 ^ 2 / 0.2e1 + t74 ^ 2 / 0.2e1 + t75 ^ 2 / 0.2e1;];
T_reg = t1;
