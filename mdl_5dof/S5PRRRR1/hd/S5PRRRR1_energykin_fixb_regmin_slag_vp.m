% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:48
% EndTime: 2019-07-18 13:28:48
% DurationCPUTime: 0.06s
% Computational Cost: add. (68->23), mult. (203->71), div. (0->0), fcn. (153->8), ass. (0->28)
t92 = qJD(2) ^ 2;
t100 = t92 / 0.2e1;
t87 = sin(qJ(2));
t99 = qJD(1) * t87;
t90 = cos(qJ(3));
t98 = qJD(2) * t90;
t91 = cos(qJ(2));
t97 = qJD(2) * t91;
t96 = qJD(3) * t87;
t95 = qJD(1) * qJD(2);
t94 = qJD(2) * qJD(3);
t93 = t90 * t99;
t85 = sin(qJ(4));
t86 = sin(qJ(3));
t89 = cos(qJ(4));
t77 = t85 * t86 * qJD(2) - t89 * t98;
t88 = cos(qJ(5));
t84 = sin(qJ(5));
t83 = qJD(3) + qJD(4);
t80 = -pkin(2) * t98 - t91 * qJD(1);
t79 = qJD(3) * pkin(2) - t86 * t99;
t78 = (t85 * t90 + t86 * t89) * qJD(2);
t76 = qJD(5) + t77;
t75 = t85 * t79 + t89 * t93;
t74 = -t89 * t79 + t85 * t93;
t73 = t88 * t78 + t84 * t83;
t72 = t84 * t78 - t88 * t83;
t1 = [qJD(1) ^ 2 / 0.2e1, t100, t91 * t95, -t87 * t95, t86 ^ 2 * t100, t86 * t92 * t90, t86 * t94, t90 * t94, qJD(3) ^ 2 / 0.2e1, (-t86 * t96 + t90 * t97) * qJD(1), (-t86 * t97 - t90 * t96) * qJD(1), t78 ^ 2 / 0.2e1, -t78 * t77, t78 * t83, -t77 * t83, t83 ^ 2 / 0.2e1, -t74 * t83 + t80 * t77, -t75 * t83 + t80 * t78, t73 ^ 2 / 0.2e1, -t73 * t72, t73 * t76, -t72 * t76, t76 ^ 2 / 0.2e1, (-t84 * t75 + t88 * t80) * t76 + t74 * t72, -(t88 * t75 + t84 * t80) * t76 + t74 * t73;];
T_reg  = t1;
