% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:05
% EndTime: 2019-12-31 17:17:05
% DurationCPUTime: 0.05s
% Computational Cost: add. (91->24), mult. (233->59), div. (0->0), fcn. (130->4), ass. (0->25)
t84 = qJD(1) ^ 2;
t94 = t84 / 0.2e1;
t93 = -pkin(6) - pkin(5);
t83 = cos(qJ(2));
t92 = t83 * t84;
t81 = sin(qJ(2));
t90 = qJD(1) * t81;
t75 = qJD(2) * pkin(2) + t93 * t90;
t89 = qJD(1) * t83;
t76 = t93 * t89;
t80 = sin(qJ(3));
t82 = cos(qJ(3));
t91 = t80 * t75 - t82 * t76;
t88 = qJD(1) * qJD(2);
t87 = t81 * t88;
t86 = t83 * t88;
t77 = (-pkin(2) * t83 - pkin(1)) * qJD(1);
t85 = t82 * t75 + t80 * t76;
t79 = qJD(2) + qJD(3);
t73 = (t80 * t83 + t81 * t82) * qJD(1);
t72 = t80 * t90 - t82 * t89;
t70 = t79 * qJ(4) + t91;
t69 = -t79 * pkin(3) + qJD(4) - t85;
t68 = t72 * pkin(3) - t73 * qJ(4) + t77;
t1 = [t94, 0, 0, t81 ^ 2 * t94, t81 * t92, t87, t86, qJD(2) ^ 2 / 0.2e1, pkin(1) * t92 - pkin(5) * t87, -t84 * pkin(1) * t81 - pkin(5) * t86, t73 ^ 2 / 0.2e1, -t73 * t72, t73 * t79, -t72 * t79, t79 ^ 2 / 0.2e1, t77 * t72 + t85 * t79, t77 * t73 - t91 * t79, t68 * t72 - t69 * t79, t69 * t73 - t70 * t72, -t68 * t73 + t70 * t79, t70 ^ 2 / 0.2e1 + t68 ^ 2 / 0.2e1 + t69 ^ 2 / 0.2e1;];
T_reg = t1;
