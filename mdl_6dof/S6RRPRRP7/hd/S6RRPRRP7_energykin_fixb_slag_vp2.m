% Calculate kinetic energy for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:40
% EndTime: 2019-03-09 12:16:40
% DurationCPUTime: 0.51s
% Computational Cost: add. (450->104), mult. (883->142), div. (0->0), fcn. (522->6), ass. (0->33)
t103 = sin(qJ(4));
t104 = sin(qJ(2));
t105 = cos(qJ(4));
t106 = cos(qJ(2));
t87 = (t103 * t104 + t105 * t106) * qJD(1);
t114 = pkin(7) * mrSges(3,3);
t113 = cos(qJ(5));
t102 = sin(qJ(5));
t110 = qJD(1) * t106;
t111 = qJD(1) * t104;
t90 = -qJD(1) * pkin(1) - pkin(2) * t110 - qJ(3) * t111;
t82 = pkin(3) * t110 - t90;
t88 = (-t103 * t106 + t104 * t105) * qJD(1);
t73 = pkin(4) * t87 - pkin(9) * t88 + t82;
t112 = pkin(7) * t111 + qJD(3);
t83 = -pkin(8) * t111 + (-pkin(2) - pkin(3)) * qJD(2) + t112;
t92 = pkin(7) * t110 + qJD(2) * qJ(3);
t89 = -pkin(8) * t110 + t92;
t78 = t103 * t83 + t105 * t89;
t99 = -qJD(2) + qJD(4);
t76 = pkin(9) * t99 + t78;
t70 = t102 * t73 + t113 * t76;
t77 = -t103 * t89 + t105 * t83;
t75 = -pkin(4) * t99 - t77;
t69 = -t102 * t76 + t113 * t73;
t91 = -qJD(2) * pkin(2) + t112;
t86 = qJD(5) + t87;
t80 = t102 * t99 + t113 * t88;
t79 = t102 * t88 - t113 * t99;
t71 = pkin(5) * t79 - qJ(6) * t80 + t75;
t68 = qJ(6) * t86 + t70;
t67 = -t86 * pkin(5) + qJD(6) - t69;
t1 = m(4) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t82 ^ 2) / 0.2e1 + m(7) * (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t75 ^ 2) / 0.2e1 + (t77 * mrSges(5,1) - t78 * mrSges(5,2) + Ifges(5,3) * t99 / 0.2e1) * t99 + (t82 * mrSges(5,2) - t77 * mrSges(5,3) + Ifges(5,5) * t99 + Ifges(5,1) * t88 / 0.2e1) * t88 - (-t82 * mrSges(5,1) + t78 * mrSges(5,3) + Ifges(5,4) * t88 + Ifges(5,6) * t99 - Ifges(5,2) * t87 / 0.2e1) * t87 + (-t91 * mrSges(4,1) + t92 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * qJD(2)) * qJD(2) + (t69 * mrSges(6,1) - t67 * mrSges(7,1) - t70 * mrSges(6,2) + t68 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t86) * t86 + (t75 * mrSges(6,2) + t67 * mrSges(7,2) - t69 * mrSges(6,3) - t71 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t80 + (Ifges(7,4) + Ifges(6,5)) * t86) * t80 + (t75 * mrSges(6,1) + t71 * mrSges(7,1) - t68 * mrSges(7,2) - t70 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t79 + (-Ifges(6,6) + Ifges(7,6)) * t86 + (-Ifges(6,4) + Ifges(7,5)) * t80) * t79 + ((-t90 * mrSges(4,1) + t92 * mrSges(4,2) + (-pkin(7) * mrSges(3,2) + Ifges(3,6) - Ifges(4,6)) * qJD(2)) * t106 + (t91 * mrSges(4,2) - t90 * mrSges(4,3) + (-pkin(7) * mrSges(3,1) + Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t104 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t104 ^ 2 + t106 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t114 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t106) * t106 + (-pkin(1) * mrSges(3,2) + (t114 + Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t104 + (Ifges(3,4) - Ifges(4,5)) * t106) * t104) * qJD(1)) * qJD(1);
T  = t1;
