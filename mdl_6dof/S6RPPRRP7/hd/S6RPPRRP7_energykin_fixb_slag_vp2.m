% Calculate kinetic energy for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:31
% EndTime: 2019-03-09 02:12:31
% DurationCPUTime: 0.38s
% Computational Cost: add. (376->83), mult. (774->122), div. (0->0), fcn. (482->6), ass. (0->32)
t96 = cos(pkin(9));
t106 = t96 ^ 2;
t100 = cos(qJ(4));
t95 = sin(pkin(9));
t98 = sin(qJ(4));
t84 = (t100 * t95 + t96 * t98) * qJD(1);
t105 = m(3) / 0.2e1;
t88 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t103 = -pkin(7) * qJD(1) + t88;
t81 = t103 * t95;
t82 = t103 * t96;
t75 = t100 * t81 + t98 * t82;
t72 = qJD(4) * pkin(8) + t75;
t85 = (t100 * t96 - t95 * t98) * qJD(1);
t90 = qJD(1) * qJ(2) + qJD(3);
t86 = t95 * qJD(1) * pkin(3) + t90;
t73 = pkin(4) * t84 - pkin(8) * t85 + t86;
t97 = sin(qJ(5));
t99 = cos(qJ(5));
t66 = t99 * t72 + t97 * t73;
t104 = t95 ^ 2 + t106;
t65 = -t72 * t97 + t99 * t73;
t74 = t100 * t82 - t98 * t81;
t71 = -qJD(4) * pkin(4) - t74;
t91 = -qJD(1) * pkin(1) + qJD(2);
t83 = qJD(5) + t84;
t77 = qJD(4) * t97 + t85 * t99;
t76 = qJD(4) * t99 - t85 * t97;
t67 = -pkin(5) * t76 + qJD(6) + t71;
t64 = qJ(6) * t76 + t66;
t63 = pkin(5) * t83 - qJ(6) * t77 + t65;
t1 = t91 ^ 2 * t105 + m(5) * (t74 ^ 2 + t75 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t104 * t88 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t71 ^ 2) / 0.2e1 + m(7) * (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) / 0.2e1 + (t86 * mrSges(5,2) - t74 * mrSges(5,3) + Ifges(5,1) * t85 / 0.2e1) * t85 - (-t86 * mrSges(5,1) + t75 * mrSges(5,3) + Ifges(5,4) * t85 - Ifges(5,2) * t84 / 0.2e1) * t84 + (t74 * mrSges(5,1) - t75 * mrSges(5,2) + Ifges(5,5) * t85 - Ifges(5,6) * t84 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t65 * mrSges(6,1) + t63 * mrSges(7,1) - t66 * mrSges(6,2) - t64 * mrSges(7,2) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t83) * t83 + (t71 * mrSges(6,2) + t67 * mrSges(7,2) - t65 * mrSges(6,3) - t63 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t77 + (Ifges(6,5) + Ifges(7,5)) * t83) * t77 + (-t71 * mrSges(6,1) - t67 * mrSges(7,1) + t66 * mrSges(6,3) + t64 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t76 + (Ifges(6,6) + Ifges(7,6)) * t83 + (Ifges(6,4) + Ifges(7,4)) * t77) * t76 + (t90 * (mrSges(4,1) * t95 + mrSges(4,2) * t96) + t91 * mrSges(3,2) - t104 * t88 * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t105 + mrSges(3,3)) * qJ(2) + Ifges(4,1) * t106 / 0.2e1 + (-Ifges(4,4) * t96 + Ifges(4,2) * t95 / 0.2e1) * t95) * qJD(1)) * qJD(1);
T  = t1;
