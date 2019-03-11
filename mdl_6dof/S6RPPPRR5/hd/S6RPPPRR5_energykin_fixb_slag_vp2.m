% Calculate kinetic energy for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:10
% EndTime: 2019-03-09 01:37:10
% DurationCPUTime: 0.18s
% Computational Cost: add. (222->67), mult. (346->101), div. (0->0), fcn. (130->6), ass. (0->28)
t95 = m(3) / 0.2e1;
t82 = qJD(1) * qJ(2) + qJD(3);
t79 = qJD(1) * pkin(3) + t82;
t80 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t86 = sin(pkin(9));
t87 = cos(pkin(9));
t73 = t86 * t79 + t87 * t80;
t71 = qJD(1) * pkin(7) + t73;
t89 = sin(qJ(5));
t91 = cos(qJ(5));
t68 = t89 * qJD(4) + t91 * t71;
t94 = qJD(1) * t89;
t93 = qJD(1) * t91;
t72 = t79 * t87 - t86 * t80;
t67 = qJD(4) * t91 - t71 * t89;
t90 = cos(qJ(6));
t88 = sin(qJ(6));
t83 = -qJD(1) * pkin(1) + qJD(2);
t81 = qJD(6) - t93;
t75 = qJD(5) * t88 + t90 * t94;
t74 = qJD(5) * t90 - t88 * t94;
t70 = -qJD(1) * pkin(4) - t72;
t66 = (-pkin(5) * t91 - pkin(8) * t89 - pkin(4)) * qJD(1) - t72;
t65 = qJD(5) * pkin(8) + t68;
t64 = -qJD(5) * pkin(5) - t67;
t63 = t65 * t90 + t66 * t88;
t62 = -t65 * t88 + t90 * t66;
t1 = t83 ^ 2 * t95 + m(5) * (qJD(4) ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(4) * (t80 ^ 2 + t82 ^ 2) / 0.2e1 + m(7) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t67 ^ 2 + t68 ^ 2 + t70 ^ 2) / 0.2e1 + (t62 * mrSges(7,1) - t63 * mrSges(7,2) + Ifges(7,3) * t81 / 0.2e1) * t81 + (t67 * mrSges(6,1) - t68 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t64 * mrSges(7,2) - t62 * mrSges(7,3) + Ifges(7,5) * t81 + Ifges(7,1) * t75 / 0.2e1) * t75 + (-t64 * mrSges(7,1) + t63 * mrSges(7,3) + Ifges(7,4) * t75 + Ifges(7,6) * t81 + Ifges(7,2) * t74 / 0.2e1) * t74 + (t82 * mrSges(4,1) + t72 * mrSges(5,1) + t83 * mrSges(3,2) - t73 * mrSges(5,2) - t80 * mrSges(4,3) + (-t70 * mrSges(6,1) + t68 * mrSges(6,3) + Ifges(6,6) * qJD(5) + Ifges(6,2) * t93 / 0.2e1) * t91 + (t70 * mrSges(6,2) - t67 * mrSges(6,3) + Ifges(6,5) * qJD(5)) * t89 + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + (qJ(2) * t95 + mrSges(3,3)) * qJ(2) + (Ifges(6,4) * t91 + Ifges(6,1) * t89 / 0.2e1) * t89) * qJD(1)) * qJD(1);
T  = t1;
