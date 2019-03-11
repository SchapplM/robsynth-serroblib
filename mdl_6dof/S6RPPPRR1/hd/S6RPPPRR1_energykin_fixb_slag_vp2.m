% Calculate kinetic energy for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:48
% EndTime: 2019-03-09 01:29:48
% DurationCPUTime: 0.25s
% Computational Cost: add. (164->67), mult. (303->99), div. (0->0), fcn. (118->6), ass. (0->28)
t84 = cos(pkin(9));
t93 = -pkin(1) * t84 - pkin(2);
t92 = -qJ(4) + t93;
t72 = -t92 * qJD(1) - qJD(3);
t96 = t72 ^ 2;
t95 = m(3) / 0.2e1;
t83 = sin(pkin(9));
t78 = (-pkin(1) * t83 - qJ(3)) * qJD(1);
t76 = qJD(4) - t78;
t71 = -qJD(1) * pkin(7) + t76;
t86 = sin(qJ(5));
t88 = cos(qJ(5));
t69 = t88 * qJD(2) + t86 * t71;
t94 = qJD(1) * t88;
t68 = -qJD(2) * t86 + t71 * t88;
t89 = qJD(2) ^ 2;
t87 = cos(qJ(6));
t85 = sin(qJ(6));
t79 = qJD(1) * t86 + qJD(6);
t77 = t93 * qJD(1) + qJD(3);
t75 = qJD(5) * t85 + t87 * t94;
t74 = qJD(5) * t87 - t85 * t94;
t67 = -qJD(3) + (pkin(5) * t86 - pkin(8) * t88 - t92) * qJD(1);
t66 = qJD(5) * pkin(8) + t69;
t65 = -qJD(5) * pkin(5) - t68;
t64 = t66 * t87 + t67 * t85;
t63 = -t66 * t85 + t87 * t67;
t1 = m(5) * (t76 ^ 2 + t89 + t96) / 0.2e1 + m(4) * (t77 ^ 2 + t78 ^ 2 + t89) / 0.2e1 + t89 * t95 + m(7) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t68 ^ 2 + t69 ^ 2 + t96) / 0.2e1 + (t63 * mrSges(7,1) - t64 * mrSges(7,2) + Ifges(7,3) * t79 / 0.2e1) * t79 + (t68 * mrSges(6,1) - t69 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t65 * mrSges(7,2) - t63 * mrSges(7,3) + Ifges(7,5) * t79 + Ifges(7,1) * t75 / 0.2e1) * t75 + (-t65 * mrSges(7,1) + t64 * mrSges(7,3) + Ifges(7,4) * t75 + Ifges(7,6) * t79 + Ifges(7,2) * t74 / 0.2e1) * t74 + (t77 * mrSges(4,2) + t76 * mrSges(5,2) - t78 * mrSges(4,3) + t72 * mrSges(5,3) + (t72 * mrSges(6,2) - t68 * mrSges(6,3) + Ifges(6,5) * qJD(5) + Ifges(6,1) * t94 / 0.2e1) * t88 + (t72 * mrSges(6,1) - t69 * mrSges(6,3) - Ifges(6,6) * qJD(5)) * t86 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + (t84 * mrSges(3,1) - t83 * mrSges(3,2) + (t83 ^ 2 + t84 ^ 2) * t95 * pkin(1)) * pkin(1) + (-Ifges(6,4) * t88 + Ifges(6,2) * t86 / 0.2e1) * t86) * qJD(1)) * qJD(1);
T  = t1;
