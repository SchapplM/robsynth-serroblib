% Calculate kinetic energy for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:34
% EndTime: 2019-03-09 08:29:34
% DurationCPUTime: 0.46s
% Computational Cost: add. (438->103), mult. (1011->137), div. (0->0), fcn. (656->6), ass. (0->35)
t106 = pkin(3) + pkin(8);
t105 = pkin(7) * mrSges(3,3);
t104 = pkin(7) + qJ(3);
t94 = sin(qJ(2));
t101 = t94 * qJD(1);
t96 = cos(qJ(2));
t102 = qJD(1) * t96;
t103 = cos(pkin(9));
t92 = sin(pkin(9));
t83 = t92 * t101 - t103 * t102;
t84 = (t103 * t94 + t92 * t96) * qJD(1);
t89 = qJD(3) + (-pkin(2) * t96 - pkin(1)) * qJD(1);
t99 = -qJ(4) * t84 + t89;
t69 = t106 * t83 + t99;
t87 = qJD(2) * pkin(2) - t104 * t101;
t88 = t104 * t102;
t77 = t103 * t87 - t92 * t88;
t100 = qJD(4) - t77;
t72 = t84 * pkin(4) - t106 * qJD(2) + t100;
t93 = sin(qJ(5));
t95 = cos(qJ(5));
t66 = t95 * t69 + t93 * t72;
t78 = t103 * t88 + t92 * t87;
t76 = -qJD(2) * qJ(4) - t78;
t65 = -t69 * t93 + t95 * t72;
t73 = -pkin(4) * t83 - t76;
t82 = qJD(5) + t84;
t80 = qJD(2) * t95 + t83 * t93;
t79 = -qJD(2) * t93 + t83 * t95;
t75 = -qJD(2) * pkin(3) + t100;
t74 = pkin(3) * t83 + t99;
t67 = -pkin(5) * t79 + qJD(6) + t73;
t64 = qJ(6) * t79 + t66;
t63 = pkin(5) * t82 - qJ(6) * t80 + t65;
t1 = m(7) * (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t77 ^ 2 + t78 ^ 2 + t89 ^ 2) / 0.2e1 + (t75 * mrSges(5,1) + t89 * mrSges(4,2) - t77 * mrSges(4,3) - t74 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t84) * t84 + (t65 * mrSges(6,1) + t63 * mrSges(7,1) - t66 * mrSges(6,2) - t64 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t82) * t82 + (t89 * mrSges(4,1) + t76 * mrSges(5,1) - t74 * mrSges(5,2) - t78 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t83 + (-Ifges(4,4) - Ifges(5,6)) * t84) * t83 + (t73 * mrSges(6,2) + t67 * mrSges(7,2) - t65 * mrSges(6,3) - t63 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t80 + (Ifges(6,5) + Ifges(7,5)) * t82) * t80 + (-t73 * mrSges(6,1) - t67 * mrSges(7,1) + t66 * mrSges(6,3) + t64 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t79 + (Ifges(6,6) + Ifges(7,6)) * t82 + (Ifges(6,4) + Ifges(7,4)) * t80) * t79 + (t77 * mrSges(4,1) - t78 * mrSges(4,2) + t75 * mrSges(5,2) - t76 * mrSges(5,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(2) + (-Ifges(5,4) + Ifges(4,5)) * t84 + (Ifges(5,5) - Ifges(4,6)) * t83 + (Ifges(3,5) * t94 + Ifges(3,6) * t96 + (-mrSges(3,1) * t94 - mrSges(3,2) * t96) * pkin(7)) * qJD(1)) * qJD(2) + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t94 ^ 2 + t96 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t105) * t96) * t96 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t96 + (Ifges(3,1) / 0.2e1 + t105) * t94) * t94) * qJD(1) ^ 2;
T  = t1;
