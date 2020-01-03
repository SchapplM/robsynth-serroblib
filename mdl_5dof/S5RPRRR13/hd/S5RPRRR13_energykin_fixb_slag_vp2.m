% Calculate kinetic energy for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:28
% EndTime: 2019-12-31 19:14:28
% DurationCPUTime: 0.26s
% Computational Cost: add. (251->71), mult. (483->115), div. (0->0), fcn. (266->6), ass. (0->30)
t78 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t92 = t78 * mrSges(4,3);
t91 = qJD(1) / 0.2e1;
t85 = sin(qJ(3));
t88 = cos(qJ(3));
t71 = (pkin(3) * t85 - pkin(7) * t88 + qJ(2)) * qJD(1);
t72 = qJD(3) * pkin(7) + t85 * t78;
t84 = sin(qJ(4));
t87 = cos(qJ(4));
t64 = t84 * t71 + t87 * t72;
t90 = t88 * qJD(1);
t79 = t85 * qJD(1) + qJD(4);
t63 = t87 * t71 - t72 * t84;
t73 = -qJD(3) * pkin(3) - t88 * t78;
t89 = qJD(1) ^ 2;
t86 = cos(qJ(5));
t83 = sin(qJ(5));
t82 = t89 * qJ(2) ^ 2;
t80 = -qJD(1) * pkin(1) + qJD(2);
t77 = qJD(5) + t79;
t75 = qJD(3) * t84 + t87 * t90;
t74 = qJD(3) * t87 - t84 * t90;
t67 = -pkin(4) * t74 + t73;
t66 = t74 * t83 + t75 * t86;
t65 = t74 * t86 - t75 * t83;
t62 = pkin(8) * t74 + t64;
t61 = pkin(4) * t79 - pkin(8) * t75 + t63;
t60 = t61 * t83 + t62 * t86;
t59 = t61 * t86 - t62 * t83;
t1 = m(3) * (t80 ^ 2 + t82) / 0.2e1 + m(4) * (t82 + (t85 ^ 2 + t88 ^ 2) * t78 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t67 ^ 2) / 0.2e1 + (qJ(2) * mrSges(3,3) + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t89 + (t63 * mrSges(5,1) - t64 * mrSges(5,2) + Ifges(5,3) * t79 / 0.2e1) * t79 + (t59 * mrSges(6,1) - t60 * mrSges(6,2) + Ifges(6,3) * t77 / 0.2e1) * t77 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t88 * mrSges(4,1) - t85 * mrSges(4,2)) * t78) * qJD(3) + (t80 * mrSges(3,2) + (qJ(2) * mrSges(4,2) * qJD(1) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t91 - t92) * t88) * t88 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t88) * qJD(1) + (Ifges(4,2) * t91 - t92) * t85) * t85) * qJD(1) + (t73 * mrSges(5,2) - t63 * mrSges(5,3) + Ifges(5,5) * t79 + Ifges(5,1) * t75 / 0.2e1) * t75 + (t67 * mrSges(6,2) - t59 * mrSges(6,3) + Ifges(6,5) * t77 + Ifges(6,1) * t66 / 0.2e1) * t66 + (-t73 * mrSges(5,1) + t64 * mrSges(5,3) + Ifges(5,4) * t75 + Ifges(5,6) * t79 + Ifges(5,2) * t74 / 0.2e1) * t74 + (-t67 * mrSges(6,1) + t60 * mrSges(6,3) + Ifges(6,4) * t66 + Ifges(6,6) * t77 + Ifges(6,2) * t65 / 0.2e1) * t65;
T = t1;
