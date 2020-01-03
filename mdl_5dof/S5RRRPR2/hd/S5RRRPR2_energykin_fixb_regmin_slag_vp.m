% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:36
% EndTime: 2020-01-03 12:07:36
% DurationCPUTime: 0.07s
% Computational Cost: add. (107->18), mult. (163->49), div. (0->0), fcn. (81->8), ass. (0->23)
t78 = qJD(1) + qJD(2);
t77 = qJD(3) + t78;
t76 = t77 ^ 2;
t93 = t76 / 0.2e1;
t91 = pkin(1) * qJD(1);
t88 = cos(qJ(2)) * t91;
t75 = t78 * pkin(2) + t88;
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t89 = sin(qJ(2)) * t91;
t87 = t85 * t75 - t82 * t89;
t71 = t77 * pkin(3) + t87;
t73 = t82 * t75 + t85 * t89;
t79 = sin(pkin(9));
t80 = cos(pkin(9));
t68 = t80 * t71 - t79 * t73;
t92 = (-t77 * pkin(4) - t68) * t77;
t69 = t79 * t71 + t80 * t73;
t90 = qJD(5) * t77;
t84 = cos(qJ(5));
t81 = sin(qJ(5));
t67 = t77 * pkin(8) + t69;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t78 ^ 2 / 0.2e1, t78 * t88, -t78 * t89, t93, t87 * t77, -t73 * t77, t69 ^ 2 / 0.2e1 + t68 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t81 ^ 2 * t93, t81 * t76 * t84, t81 * t90, t84 * t90, qJD(5) ^ 2 / 0.2e1, -t84 * t92 + (t84 * qJD(4) - t81 * t67) * qJD(5), t81 * t92 - (t81 * qJD(4) + t84 * t67) * qJD(5);];
T_reg = t1;
