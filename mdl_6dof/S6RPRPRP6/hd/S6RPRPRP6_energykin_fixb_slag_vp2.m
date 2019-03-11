% Calculate kinetic energy for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:17:53
% EndTime: 2019-03-09 03:17:53
% DurationCPUTime: 0.50s
% Computational Cost: add. (413->97), mult. (979->130), div. (0->0), fcn. (652->6), ass. (0->37)
t97 = cos(pkin(9));
t112 = t97 ^ 2;
t111 = m(3) / 0.2e1;
t110 = pkin(3) + pkin(8);
t109 = cos(qJ(3));
t108 = pkin(7) + qJ(2);
t100 = cos(qJ(5));
t96 = sin(pkin(9));
t99 = sin(qJ(3));
t87 = (t109 * t96 + t97 * t99) * qJD(1);
t90 = qJD(2) + (-pkin(2) * t97 - pkin(1)) * qJD(1);
t103 = -qJ(4) * t87 + t90;
t105 = t96 * qJD(1);
t106 = qJD(1) * t97;
t86 = t99 * t105 - t109 * t106;
t70 = t110 * t86 + t103;
t88 = t108 * t105;
t89 = t108 * t106;
t78 = -t109 * t88 - t99 * t89;
t104 = qJD(4) - t78;
t73 = t87 * pkin(4) - t110 * qJD(3) + t104;
t98 = sin(qJ(5));
t67 = t100 * t70 + t98 * t73;
t79 = t109 * t89 - t99 * t88;
t77 = -qJD(3) * qJ(4) - t79;
t66 = t100 * t73 - t70 * t98;
t74 = -pkin(4) * t86 - t77;
t92 = -qJD(1) * pkin(1) + qJD(2);
t82 = qJD(5) + t87;
t81 = qJD(3) * t100 + t86 * t98;
t80 = -qJD(3) * t98 + t100 * t86;
t76 = -qJD(3) * pkin(3) + t104;
t75 = pkin(3) * t86 + t103;
t68 = -pkin(5) * t80 + qJD(6) + t74;
t65 = qJ(6) * t80 + t67;
t64 = pkin(5) * t82 - qJ(6) * t81 + t66;
t1 = m(4) * (t78 ^ 2 + t79 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(7) * (t64 ^ 2 + t65 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t66 ^ 2 + t67 ^ 2 + t74 ^ 2) / 0.2e1 + t92 ^ 2 * t111 + (t76 * mrSges(5,1) + t90 * mrSges(4,2) - t78 * mrSges(4,3) - t75 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t87) * t87 + (t66 * mrSges(6,1) + t64 * mrSges(7,1) - t67 * mrSges(6,2) - t65 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t82) * t82 + (t90 * mrSges(4,1) + t77 * mrSges(5,1) - t75 * mrSges(5,2) - t79 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t86 + (-Ifges(4,4) - Ifges(5,6)) * t87) * t86 + (t74 * mrSges(6,2) + t68 * mrSges(7,2) - t66 * mrSges(6,3) - t64 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t81 + (Ifges(6,5) + Ifges(7,5)) * t82) * t81 + (-t74 * mrSges(6,1) - t68 * mrSges(7,1) + t67 * mrSges(6,3) + t65 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t80 + (Ifges(6,6) + Ifges(7,6)) * t82 + (Ifges(6,4) + Ifges(7,4)) * t81) * t80 + (t78 * mrSges(4,1) - t79 * mrSges(4,2) + t76 * mrSges(5,2) - t77 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(3) + (-Ifges(5,4) + Ifges(4,5)) * t87 + (Ifges(5,5) - Ifges(4,6)) * t86) * qJD(3) + (t92 * (-mrSges(3,1) * t97 + mrSges(3,2) * t96) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t111 + mrSges(3,3)) * (t96 ^ 2 + t112) * qJ(2) + Ifges(3,2) * t112 / 0.2e1 + (Ifges(3,4) * t97 + Ifges(3,1) * t96 / 0.2e1) * t96) * qJD(1)) * qJD(1);
T  = t1;
