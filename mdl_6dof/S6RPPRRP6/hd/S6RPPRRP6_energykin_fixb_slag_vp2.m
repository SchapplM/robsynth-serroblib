% Calculate kinetic energy for
% S6RPPRRP6
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:08
% EndTime: 2019-03-09 02:10:08
% DurationCPUTime: 0.26s
% Computational Cost: add. (214->75), mult. (373->104), div. (0->0), fcn. (160->4), ass. (0->26)
t84 = -pkin(1) - qJ(3);
t72 = -t84 * qJD(1) - qJD(2);
t89 = t72 ^ 2;
t88 = m(3) / 0.2e1;
t87 = cos(qJ(5));
t75 = qJD(1) * qJ(2) + qJD(3);
t71 = -qJD(1) * pkin(7) + t75;
t86 = t71 * mrSges(5,3);
t85 = qJD(1) / 0.2e1;
t80 = sin(qJ(4));
t81 = cos(qJ(4));
t64 = -qJD(2) + (pkin(4) * t80 - pkin(8) * t81 - t84) * qJD(1);
t66 = qJD(4) * pkin(8) + t80 * t71;
t79 = sin(qJ(5));
t62 = t79 * t64 + t87 * t66;
t83 = qJD(1) * t81;
t67 = -qJD(4) * pkin(4) - t81 * t71;
t61 = t87 * t64 - t79 * t66;
t76 = -qJD(1) * pkin(1) + qJD(2);
t74 = qJD(1) * t80 + qJD(5);
t69 = t79 * qJD(4) + t87 * t83;
t68 = -t87 * qJD(4) + t79 * t83;
t60 = pkin(5) * t68 - qJ(6) * t69 + t67;
t59 = qJ(6) * t74 + t62;
t58 = -t74 * pkin(5) + qJD(6) - t61;
t1 = m(5) * (t89 + (t80 ^ 2 + t81 ^ 2) * t71 ^ 2) / 0.2e1 + t76 ^ 2 * t88 + m(4) * (t75 ^ 2 + t89) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t67 ^ 2) / 0.2e1 + m(7) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + (Ifges(5,3) * qJD(4) / 0.2e1 + (t81 * mrSges(5,1) - t80 * mrSges(5,2)) * t71) * qJD(4) + (t61 * mrSges(6,1) - t58 * mrSges(7,1) - t62 * mrSges(6,2) + t59 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t74) * t74 + (t67 * mrSges(6,2) + t58 * mrSges(7,2) - t61 * mrSges(6,3) - t60 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t69 + (Ifges(7,4) + Ifges(6,5)) * t74) * t69 + (t67 * mrSges(6,1) + t60 * mrSges(7,1) - t59 * mrSges(7,2) - t62 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t68 + (-Ifges(6,6) + Ifges(7,6)) * t74 + (-Ifges(6,4) + Ifges(7,5)) * t69) * t68 + (t76 * mrSges(3,2) + t75 * mrSges(4,2) + t72 * mrSges(4,3) + (t72 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t85 - t86) * t81) * t81 + (-Ifges(5,4) * t83 + t72 * mrSges(5,1) - Ifges(5,6) * qJD(4) + (Ifges(5,2) * t85 - t86) * t80) * t80 + (Ifges(2,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t88 + mrSges(3,3)) * qJ(2)) * qJD(1)) * qJD(1);
T  = t1;
