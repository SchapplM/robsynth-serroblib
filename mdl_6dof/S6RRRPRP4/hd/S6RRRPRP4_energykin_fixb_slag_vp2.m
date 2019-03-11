% Calculate kinetic energy for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:29
% EndTime: 2019-03-09 16:43:30
% DurationCPUTime: 0.49s
% Computational Cost: add. (490->103), mult. (1007->138), div. (0->0), fcn. (652->6), ass. (0->36)
t109 = pkin(3) + pkin(9);
t108 = -pkin(8) - pkin(7);
t107 = pkin(7) * mrSges(3,3);
t106 = cos(qJ(3));
t105 = cos(qJ(5));
t97 = sin(qJ(2));
t104 = t97 * qJD(1);
t89 = qJD(2) * pkin(2) + t108 * t104;
t98 = cos(qJ(2));
t103 = t98 * qJD(1);
t90 = t108 * t103;
t96 = sin(qJ(3));
t78 = t106 * t89 + t96 * t90;
t102 = qJD(4) - t78;
t86 = (t106 * t97 + t96 * t98) * qJD(1);
t94 = qJD(2) + qJD(3);
t72 = t86 * pkin(4) - t109 * t94 + t102;
t91 = (-pkin(2) * t98 - pkin(1)) * qJD(1);
t101 = -qJ(4) * t86 + t91;
t85 = -t106 * t103 + t96 * t104;
t73 = t109 * t85 + t101;
t95 = sin(qJ(5));
t68 = t105 * t73 + t95 * t72;
t79 = -t106 * t90 + t96 * t89;
t77 = -t94 * qJ(4) - t79;
t74 = -pkin(4) * t85 - t77;
t67 = t105 * t72 - t95 * t73;
t84 = qJD(5) + t86;
t81 = t105 * t94 + t95 * t85;
t80 = -t105 * t85 + t94 * t95;
t76 = -t94 * pkin(3) + t102;
t75 = pkin(3) * t85 + t101;
t69 = pkin(5) * t80 - qJ(6) * t81 + t74;
t66 = qJ(6) * t84 + t68;
t65 = -t84 * pkin(5) + qJD(6) - t67;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t78 ^ 2 + t79 ^ 2 + t91 ^ 2) / 0.2e1 + m(6) * (t67 ^ 2 + t68 ^ 2 + t74 ^ 2) / 0.2e1 + m(5) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(7) * (t65 ^ 2 + t66 ^ 2 + t69 ^ 2) / 0.2e1 + (t78 * mrSges(4,1) - t79 * mrSges(4,2) + t76 * mrSges(5,2) - t77 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t94) * t94 + (t67 * mrSges(6,1) - t65 * mrSges(7,1) - t68 * mrSges(6,2) + t66 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t84) * t84 + (t76 * mrSges(5,1) + t91 * mrSges(4,2) - t78 * mrSges(4,3) - t75 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t86 + (-Ifges(5,4) + Ifges(4,5)) * t94) * t86 + (t74 * mrSges(6,2) + t65 * mrSges(7,2) - t67 * mrSges(6,3) - t69 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t81 + (Ifges(7,4) + Ifges(6,5)) * t84) * t81 + (t91 * mrSges(4,1) + t77 * mrSges(5,1) - t75 * mrSges(5,2) - t79 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t85 + (Ifges(5,5) - Ifges(4,6)) * t94 + (-Ifges(4,4) - Ifges(5,6)) * t86) * t85 + (t74 * mrSges(6,1) + t69 * mrSges(7,1) - t66 * mrSges(7,2) - t68 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t80 + (-Ifges(6,6) + Ifges(7,6)) * t84 + (-Ifges(6,4) + Ifges(7,5)) * t81) * t80 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t97 ^ 2 + t98 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t107 + Ifges(3,2) / 0.2e1) * t98) * t98 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t98 + (t107 + Ifges(3,1) / 0.2e1) * t97) * t97) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t98 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t97) * qJD(2)) * qJD(1);
T  = t1;
