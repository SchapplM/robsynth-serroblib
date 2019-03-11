% Calculate kinetic energy for
% S6RPPRRP8
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:11
% EndTime: 2019-03-09 02:15:12
% DurationCPUTime: 0.40s
% Computational Cost: add. (376->83), mult. (770->122), div. (0->0), fcn. (478->6), ass. (0->32)
t96 = cos(pkin(9));
t106 = t96 ^ 2;
t95 = sin(pkin(9));
t98 = sin(qJ(4));
t99 = cos(qJ(4));
t83 = (t95 * t99 + t96 * t98) * qJD(1);
t105 = m(3) / 0.2e1;
t104 = cos(qJ(5));
t87 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t102 = -pkin(7) * qJD(1) + t87;
t80 = t102 * t95;
t81 = t102 * t96;
t74 = t99 * t80 + t98 * t81;
t71 = qJD(4) * pkin(8) + t74;
t84 = (-t95 * t98 + t96 * t99) * qJD(1);
t89 = qJD(1) * qJ(2) + qJD(3);
t85 = t95 * qJD(1) * pkin(3) + t89;
t72 = pkin(4) * t83 - pkin(8) * t84 + t85;
t97 = sin(qJ(5));
t66 = t104 * t71 + t97 * t72;
t103 = t95 ^ 2 + t106;
t73 = -t98 * t80 + t81 * t99;
t70 = -qJD(4) * pkin(4) - t73;
t65 = t104 * t72 - t97 * t71;
t90 = -qJD(1) * pkin(1) + qJD(2);
t82 = qJD(5) + t83;
t76 = t97 * qJD(4) + t104 * t84;
t75 = -t104 * qJD(4) + t84 * t97;
t67 = pkin(5) * t75 - qJ(6) * t76 + t70;
t64 = qJ(6) * t82 + t66;
t63 = -t82 * pkin(5) + qJD(6) - t65;
t1 = t90 ^ 2 * t105 + m(5) * (t73 ^ 2 + t74 ^ 2 + t85 ^ 2) / 0.2e1 + m(4) * (t103 * t87 ^ 2 + t89 ^ 2) / 0.2e1 + m(7) * (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t70 ^ 2) / 0.2e1 + (t85 * mrSges(5,2) - t73 * mrSges(5,3) + Ifges(5,1) * t84 / 0.2e1) * t84 - (-t85 * mrSges(5,1) + t74 * mrSges(5,3) + Ifges(5,4) * t84 - Ifges(5,2) * t83 / 0.2e1) * t83 + (t73 * mrSges(5,1) - t74 * mrSges(5,2) + Ifges(5,5) * t84 - Ifges(5,6) * t83 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t65 * mrSges(6,1) - t63 * mrSges(7,1) - t66 * mrSges(6,2) + t64 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t82) * t82 + (t70 * mrSges(6,2) + t63 * mrSges(7,2) - t65 * mrSges(6,3) - t67 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t76 + (Ifges(7,4) + Ifges(6,5)) * t82) * t76 + (t70 * mrSges(6,1) + t67 * mrSges(7,1) - t64 * mrSges(7,2) - t66 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t75 + (-Ifges(6,6) + Ifges(7,6)) * t82 + (-Ifges(6,4) + Ifges(7,5)) * t76) * t75 + (t89 * (mrSges(4,1) * t95 + mrSges(4,2) * t96) + t90 * mrSges(3,2) - t103 * t87 * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t105 + mrSges(3,3)) * qJ(2) + Ifges(4,1) * t106 / 0.2e1 + (-Ifges(4,4) * t96 + Ifges(4,2) * t95 / 0.2e1) * t95) * qJD(1)) * qJD(1);
T  = t1;
