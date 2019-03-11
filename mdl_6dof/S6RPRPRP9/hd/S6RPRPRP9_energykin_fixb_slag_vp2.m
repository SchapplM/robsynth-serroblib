% Calculate kinetic energy for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:29
% EndTime: 2019-03-09 03:27:30
% DurationCPUTime: 0.34s
% Computational Cost: add. (409->91), mult. (810->132), div. (0->0), fcn. (488->6), ass. (0->32)
t101 = cos(qJ(5));
t88 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t100 = t88 * mrSges(4,3);
t95 = sin(qJ(3));
t96 = cos(qJ(3));
t83 = (pkin(3) * t95 - qJ(4) * t96 + qJ(2)) * qJD(1);
t84 = qJD(3) * qJ(4) + t95 * t88;
t92 = sin(pkin(9));
t93 = cos(pkin(9));
t73 = t93 * t83 - t84 * t92;
t98 = t96 * qJD(1);
t86 = qJD(3) * t92 + t93 * t98;
t99 = t95 * qJD(1);
t70 = pkin(4) * t99 - pkin(8) * t86 + t73;
t74 = t92 * t83 + t93 * t84;
t85 = qJD(3) * t93 - t92 * t98;
t72 = pkin(8) * t85 + t74;
t94 = sin(qJ(5));
t67 = t101 * t72 + t94 * t70;
t66 = t101 * t70 - t94 * t72;
t82 = -qJD(3) * pkin(3) - t96 * t88 + qJD(4);
t77 = -pkin(4) * t85 + t82;
t97 = qJD(1) ^ 2;
t91 = t97 * qJ(2) ^ 2;
t90 = -qJD(1) * pkin(1) + qJD(2);
t89 = qJD(5) + t99;
t76 = t101 * t86 + t94 * t85;
t75 = -t101 * t85 + t86 * t94;
t68 = pkin(5) * t75 - qJ(6) * t76 + t77;
t65 = qJ(6) * t89 + t67;
t64 = -t89 * pkin(5) + qJD(6) - t66;
t1 = m(7) * (t64 ^ 2 + t65 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t66 ^ 2 + t67 ^ 2 + t77 ^ 2) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t82 ^ 2) / 0.2e1 + m(3) * (t90 ^ 2 + t91) / 0.2e1 + m(4) * (t91 + (t95 ^ 2 + t96 ^ 2) * t88 ^ 2) / 0.2e1 + (qJ(2) * mrSges(3,3) + Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * t97 + (t82 * mrSges(5,2) - t73 * mrSges(5,3) + Ifges(5,1) * t86 / 0.2e1) * t86 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t96 * mrSges(4,1) - t95 * mrSges(4,2)) * t88) * qJD(3) + (t90 * mrSges(3,2) + (qJ(2) * mrSges(4,2) * qJD(1) + Ifges(4,5) * qJD(3) + (-t100 + Ifges(4,1) * qJD(1) / 0.2e1) * t96) * t96 + (Ifges(5,5) * t86 - Ifges(4,6) * qJD(3) - t74 * mrSges(5,2) + t73 * mrSges(5,1) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t96) * qJD(1) + (-t100 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(1)) * t95) * t95) * qJD(1) + (Ifges(5,6) * t99 - t82 * mrSges(5,1) + t74 * mrSges(5,3) + Ifges(5,4) * t86 + Ifges(5,2) * t85 / 0.2e1) * t85 + (t66 * mrSges(6,1) - t64 * mrSges(7,1) - t67 * mrSges(6,2) + t65 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t89) * t89 + (t77 * mrSges(6,2) + t64 * mrSges(7,2) - t66 * mrSges(6,3) - t68 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t76 + (Ifges(7,4) + Ifges(6,5)) * t89) * t76 + (t77 * mrSges(6,1) + t68 * mrSges(7,1) - t65 * mrSges(7,2) - t67 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t75 + (-Ifges(6,6) + Ifges(7,6)) * t89 + (-Ifges(6,4) + Ifges(7,5)) * t76) * t75;
T  = t1;
