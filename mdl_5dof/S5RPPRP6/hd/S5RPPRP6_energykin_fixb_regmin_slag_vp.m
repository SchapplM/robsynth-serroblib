% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:16
% EndTime: 2019-12-31 17:55:16
% DurationCPUTime: 0.06s
% Computational Cost: add. (122->28), mult. (258->66), div. (0->0), fcn. (127->4), ass. (0->24)
t95 = qJD(1) ^ 2;
t100 = t95 / 0.2e1;
t91 = sin(pkin(7));
t84 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t97 = -pkin(6) * qJD(1) + t84;
t78 = t97 * t91;
t92 = cos(pkin(7));
t79 = t97 * t92;
t93 = sin(qJ(4));
t94 = cos(qJ(4));
t99 = t94 * t78 + t93 * t79;
t98 = qJD(1) * t91;
t86 = qJD(1) * qJ(2) + qJD(3);
t82 = pkin(3) * t98 + t86;
t96 = -t93 * t78 + t94 * t79;
t89 = t92 ^ 2;
t88 = t91 ^ 2;
t87 = -qJD(1) * pkin(1) + qJD(2);
t81 = (-t91 * t93 + t92 * t94) * qJD(1);
t80 = (t91 * t94 + t92 * t93) * qJD(1);
t75 = t80 * pkin(4) - t81 * qJ(5) + t82;
t74 = qJD(4) * qJ(5) + t99;
t73 = -qJD(4) * pkin(4) + qJD(5) - t96;
t1 = [t100, 0, 0, t87 * qJD(1), t95 * qJ(2), qJ(2) ^ 2 * t100 + t87 ^ 2 / 0.2e1, t86 * t98, t86 * t92 * qJD(1), (-t88 - t89) * t84 * qJD(1), t86 ^ 2 / 0.2e1 + (t88 / 0.2e1 + t89 / 0.2e1) * t84 ^ 2, t81 ^ 2 / 0.2e1, -t81 * t80, t81 * qJD(4), -t80 * qJD(4), qJD(4) ^ 2 / 0.2e1, t96 * qJD(4) + t82 * t80, -t99 * qJD(4) + t82 * t81, -t73 * qJD(4) + t75 * t80, t73 * t81 - t74 * t80, t74 * qJD(4) - t75 * t81, t74 ^ 2 / 0.2e1 + t75 ^ 2 / 0.2e1 + t73 ^ 2 / 0.2e1;];
T_reg = t1;
