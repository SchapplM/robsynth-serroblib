% Calculate kinetic energy for
% S5RRRRR9
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:39
% EndTime: 2019-12-31 22:27:39
% DurationCPUTime: 0.43s
% Computational Cost: add. (506->87), mult. (1076->138), div. (0->0), fcn. (740->8), ass. (0->36)
t107 = pkin(6) * mrSges(3,3);
t100 = cos(qJ(4));
t101 = cos(qJ(3));
t102 = cos(qJ(2));
t98 = sin(qJ(2));
t85 = (-pkin(2) * t102 - pkin(7) * t98 - pkin(1)) * qJD(1);
t105 = qJD(1) * t102;
t91 = pkin(6) * t105 + qJD(2) * pkin(7);
t97 = sin(qJ(3));
t80 = t101 * t85 - t91 * t97;
t106 = qJD(1) * t98;
t87 = qJD(2) * t97 + t101 * t106;
t93 = qJD(3) - t105;
t75 = pkin(3) * t93 - pkin(8) * t87 + t80;
t81 = t101 * t91 + t97 * t85;
t86 = qJD(2) * t101 - t97 * t106;
t77 = pkin(8) * t86 + t81;
t96 = sin(qJ(4));
t69 = t100 * t77 + t96 * t75;
t68 = t100 * t75 - t77 * t96;
t90 = -qJD(2) * pkin(2) + pkin(6) * t106;
t92 = qJD(4) + t93;
t82 = -pkin(3) * t86 + t90;
t99 = cos(qJ(5));
t95 = sin(qJ(5));
t89 = qJD(5) + t92;
t79 = t100 * t87 + t86 * t96;
t78 = t100 * t86 - t87 * t96;
t72 = -pkin(4) * t78 + t82;
t71 = t78 * t95 + t79 * t99;
t70 = t78 * t99 - t79 * t95;
t67 = pkin(9) * t78 + t69;
t66 = pkin(4) * t92 - pkin(9) * t79 + t68;
t65 = t66 * t95 + t67 * t99;
t64 = t66 * t99 - t67 * t95;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t80 ^ 2 + t81 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t72 ^ 2) / 0.2e1 + (t80 * mrSges(4,1) - t81 * mrSges(4,2) + Ifges(4,3) * t93 / 0.2e1) * t93 + (t68 * mrSges(5,1) - t69 * mrSges(5,2) + Ifges(5,3) * t92 / 0.2e1) * t92 + (t64 * mrSges(6,1) - t65 * mrSges(6,2) + Ifges(6,3) * t89 / 0.2e1) * t89 + (t90 * mrSges(4,2) - t80 * mrSges(4,3) + Ifges(4,5) * t93 + Ifges(4,1) * t87 / 0.2e1) * t87 + (t82 * mrSges(5,2) - t68 * mrSges(5,3) + Ifges(5,5) * t92 + Ifges(5,1) * t79 / 0.2e1) * t79 + (t72 * mrSges(6,2) - t64 * mrSges(6,3) + Ifges(6,5) * t89 + Ifges(6,1) * t71 / 0.2e1) * t71 + (-t90 * mrSges(4,1) + t81 * mrSges(4,3) + Ifges(4,4) * t87 + Ifges(4,6) * t93 + Ifges(4,2) * t86 / 0.2e1) * t86 + (-t82 * mrSges(5,1) + t69 * mrSges(5,3) + Ifges(5,4) * t79 + Ifges(5,6) * t92 + Ifges(5,2) * t78 / 0.2e1) * t78 + (-t72 * mrSges(6,1) + t65 * mrSges(6,3) + Ifges(6,4) * t71 + Ifges(6,6) * t89 + Ifges(6,2) * t70 / 0.2e1) * t70 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t102 ^ 2 + t98 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (-pkin(1) * mrSges(3,2) + (t107 + Ifges(3,1) / 0.2e1) * t98) * t98 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t98 + (t107 + Ifges(3,2) / 0.2e1) * t102) * t102) * qJD(1) + ((-pkin(6) * mrSges(3,1) + Ifges(3,5)) * t98 + (-pkin(6) * mrSges(3,2) + Ifges(3,6)) * t102) * qJD(2)) * qJD(1);
T = t1;
