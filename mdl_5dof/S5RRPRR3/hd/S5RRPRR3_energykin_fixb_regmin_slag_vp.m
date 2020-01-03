% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:34
% EndTime: 2020-01-03 12:00:34
% DurationCPUTime: 0.09s
% Computational Cost: add. (101->18), mult. (163->49), div. (0->0), fcn. (81->8), ass. (0->23)
t78 = qJD(1) + qJD(2);
t77 = qJD(4) + t78;
t76 = t77 ^ 2;
t94 = t76 / 0.2e1;
t91 = pkin(1) * qJD(1);
t88 = cos(qJ(2)) * t91;
t75 = t78 * pkin(2) + t88;
t79 = sin(pkin(9));
t80 = cos(pkin(9));
t89 = sin(qJ(2)) * t91;
t72 = t80 * t75 - t79 * t89;
t70 = t78 * pkin(3) + t72;
t73 = t79 * t75 + t80 * t89;
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t87 = t85 * t70 - t82 * t73;
t93 = (-t77 * pkin(4) - t87) * t77;
t92 = t82 * t70 + t85 * t73;
t90 = qJD(5) * t77;
t84 = cos(qJ(5));
t81 = sin(qJ(5));
t68 = pkin(8) * t77 + t92;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t78 ^ 2 / 0.2e1, t78 * t88, -t78 * t89, t73 ^ 2 / 0.2e1 + t72 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t94, t87 * t77, -t92 * t77, t81 ^ 2 * t94, t81 * t76 * t84, t81 * t90, t84 * t90, qJD(5) ^ 2 / 0.2e1, -t84 * t93 + (t84 * qJD(3) - t81 * t68) * qJD(5), t81 * t93 - (t81 * qJD(3) + t84 * t68) * qJD(5);];
T_reg = t1;
