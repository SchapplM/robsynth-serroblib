% Calculate kinetic energy for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:37
% EndTime: 2019-03-09 02:07:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (214->75), mult. (377->104), div. (0->0), fcn. (164->4), ass. (0->26)
t85 = -pkin(1) - qJ(3);
t73 = -t85 * qJD(1) - qJD(2);
t89 = t73 ^ 2;
t88 = m(3) / 0.2e1;
t76 = qJD(1) * qJ(2) + qJD(3);
t72 = -qJD(1) * pkin(7) + t76;
t87 = t72 * mrSges(5,3);
t86 = qJD(1) / 0.2e1;
t80 = sin(qJ(4));
t82 = cos(qJ(4));
t65 = -qJD(2) + (pkin(4) * t80 - pkin(8) * t82 - t85) * qJD(1);
t67 = qJD(4) * pkin(8) + t80 * t72;
t79 = sin(qJ(5));
t81 = cos(qJ(5));
t61 = t79 * t65 + t81 * t67;
t84 = t82 * qJD(1);
t60 = t81 * t65 - t67 * t79;
t68 = -qJD(4) * pkin(4) - t82 * t72;
t77 = -qJD(1) * pkin(1) + qJD(2);
t75 = t80 * qJD(1) + qJD(5);
t70 = qJD(4) * t79 + t81 * t84;
t69 = qJD(4) * t81 - t79 * t84;
t62 = -pkin(5) * t69 + qJD(6) + t68;
t59 = qJ(6) * t69 + t61;
t58 = pkin(5) * t75 - qJ(6) * t70 + t60;
t1 = m(5) * (t89 + (t80 ^ 2 + t82 ^ 2) * t72 ^ 2) / 0.2e1 + t77 ^ 2 * t88 + m(4) * (t76 ^ 2 + t89) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t68 ^ 2) / 0.2e1 + m(7) * (t58 ^ 2 + t59 ^ 2 + t62 ^ 2) / 0.2e1 + (Ifges(5,3) * qJD(4) / 0.2e1 + (t82 * mrSges(5,1) - t80 * mrSges(5,2)) * t72) * qJD(4) + (t60 * mrSges(6,1) + t58 * mrSges(7,1) - t61 * mrSges(6,2) - t59 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t75) * t75 + (t68 * mrSges(6,2) + t62 * mrSges(7,2) - t60 * mrSges(6,3) - t58 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t70 + (Ifges(6,5) + Ifges(7,5)) * t75) * t70 + (-t68 * mrSges(6,1) - t62 * mrSges(7,1) + t61 * mrSges(6,3) + t59 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t69 + (Ifges(6,6) + Ifges(7,6)) * t75 + (Ifges(6,4) + Ifges(7,4)) * t70) * t69 + (t77 * mrSges(3,2) + t76 * mrSges(4,2) + t73 * mrSges(4,3) + (t73 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t86 - t87) * t82) * t82 + (-Ifges(5,4) * t84 + t73 * mrSges(5,1) - Ifges(5,6) * qJD(4) + (Ifges(5,2) * t86 - t87) * t80) * t80 + (Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t88 + mrSges(3,3)) * qJ(2)) * qJD(1)) * qJD(1);
T  = t1;
