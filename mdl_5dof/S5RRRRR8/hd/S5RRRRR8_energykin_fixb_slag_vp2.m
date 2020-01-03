% Calculate kinetic energy for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:19
% EndTime: 2019-12-31 22:24:19
% DurationCPUTime: 0.44s
% Computational Cost: add. (482->86), mult. (1002->137), div. (0->0), fcn. (690->8), ass. (0->37)
t109 = -pkin(7) - pkin(6);
t108 = pkin(6) * mrSges(3,3);
t101 = cos(qJ(4));
t102 = cos(qJ(3));
t103 = cos(qJ(2));
t106 = qJD(1) * t103;
t99 = sin(qJ(2));
t107 = qJD(1) * t99;
t98 = sin(qJ(3));
t87 = t102 * t106 - t98 * t107;
t88 = (t102 * t99 + t103 * t98) * qJD(1);
t93 = (-pkin(2) * t103 - pkin(1)) * qJD(1);
t76 = -pkin(3) * t87 - pkin(8) * t88 + t93;
t91 = qJD(2) * pkin(2) + t109 * t107;
t92 = t109 * t106;
t81 = -t102 * t92 + t98 * t91;
t95 = qJD(2) + qJD(3);
t79 = pkin(8) * t95 + t81;
t97 = sin(qJ(4));
t70 = t101 * t79 + t97 * t76;
t69 = t101 * t76 - t79 * t97;
t80 = t102 * t91 + t98 * t92;
t78 = -pkin(3) * t95 - t80;
t86 = qJD(4) - t87;
t100 = cos(qJ(5));
t96 = sin(qJ(5));
t84 = qJD(5) + t86;
t83 = t101 * t88 + t95 * t97;
t82 = t101 * t95 - t88 * t97;
t73 = t100 * t83 + t82 * t96;
t72 = t100 * t82 - t83 * t96;
t71 = -pkin(4) * t82 + t78;
t68 = pkin(9) * t82 + t70;
t67 = pkin(4) * t86 - pkin(9) * t83 + t69;
t66 = t100 * t68 + t67 * t96;
t65 = t100 * t67 - t68 * t96;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t80 ^ 2 + t81 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t69 ^ 2 + t70 ^ 2 + t78 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t71 ^ 2) / 0.2e1 + (t80 * mrSges(4,1) - t81 * mrSges(4,2) + Ifges(4,3) * t95 / 0.2e1) * t95 + (t69 * mrSges(5,1) - t70 * mrSges(5,2) + Ifges(5,3) * t86 / 0.2e1) * t86 + (t65 * mrSges(6,1) - t66 * mrSges(6,2) + Ifges(6,3) * t84 / 0.2e1) * t84 + (t93 * mrSges(4,2) - t80 * mrSges(4,3) + Ifges(4,5) * t95 + Ifges(4,1) * t88 / 0.2e1) * t88 + (t78 * mrSges(5,2) - t69 * mrSges(5,3) + Ifges(5,5) * t86 + Ifges(5,1) * t83 / 0.2e1) * t83 + (t71 * mrSges(6,2) - t65 * mrSges(6,3) + Ifges(6,5) * t84 + Ifges(6,1) * t73 / 0.2e1) * t73 + (-t93 * mrSges(4,1) + t81 * mrSges(4,3) + Ifges(4,4) * t88 + Ifges(4,6) * t95 + Ifges(4,2) * t87 / 0.2e1) * t87 + (-t78 * mrSges(5,1) + t70 * mrSges(5,3) + Ifges(5,4) * t83 + Ifges(5,6) * t86 + Ifges(5,2) * t82 / 0.2e1) * t82 + (-t71 * mrSges(6,1) + t66 * mrSges(6,3) + Ifges(6,4) * t73 + Ifges(6,6) * t84 + Ifges(6,2) * t72 / 0.2e1) * t72 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t103 ^ 2 + t99 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (-pkin(1) * mrSges(3,2) + (t108 + Ifges(3,1) / 0.2e1) * t99) * t99 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t99 + (t108 + Ifges(3,2) / 0.2e1) * t103) * t103) * qJD(1) + ((-pkin(6) * mrSges(3,1) + Ifges(3,5)) * t99 + (-pkin(6) * mrSges(3,2) + Ifges(3,6)) * t103) * qJD(2)) * qJD(1);
T = t1;
