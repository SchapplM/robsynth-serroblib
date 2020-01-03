% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:43
% EndTime: 2019-12-31 16:40:43
% DurationCPUTime: 0.07s
% Computational Cost: add. (48->27), mult. (143->53), div. (0->0), fcn. (65->4), ass. (0->24)
t86 = qJD(1) ^ 2;
t93 = t86 / 0.2e1;
t92 = qJ(2) * t86;
t82 = sin(pkin(6));
t91 = qJD(1) * t82;
t83 = cos(pkin(6));
t90 = qJD(1) * t83;
t75 = qJ(2) * t91 + qJD(3);
t89 = qJ(2) ^ 2 * t93;
t88 = qJ(3) * t82 + pkin(1);
t85 = cos(qJ(4));
t84 = sin(qJ(4));
t81 = t83 ^ 2;
t80 = t82 ^ 2;
t79 = -qJD(1) * pkin(1) + qJD(2);
t78 = t81 * t92;
t76 = t81 * t89;
t74 = (-pkin(5) + qJ(2)) * t90;
t73 = -pkin(5) * t91 + t75;
t72 = (t82 * t85 - t83 * t84) * qJD(1);
t71 = (t82 * t84 + t83 * t85) * qJD(1);
t70 = qJD(2) + (-pkin(2) * t83 - t88) * qJD(1);
t69 = -qJD(2) + ((pkin(2) + pkin(3)) * t83 + t88) * qJD(1);
t1 = [t93, 0, 0, -t79 * t90, t79 * t91, t80 * t92 + t78, t76 + t80 * t89 + t79 ^ 2 / 0.2e1, -t70 * t90, t75 * t91 + t78, -t70 * t91, t76 + t70 ^ 2 / 0.2e1 + t75 ^ 2 / 0.2e1, t72 ^ 2 / 0.2e1, -t72 * t71, t72 * qJD(4), -t71 * qJD(4), qJD(4) ^ 2 / 0.2e1, t69 * t71 + (t73 * t85 - t74 * t84) * qJD(4), t69 * t72 - (t73 * t84 + t74 * t85) * qJD(4);];
T_reg = t1;
