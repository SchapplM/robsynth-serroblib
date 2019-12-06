% Calculate kinetic energy for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:42
% EndTime: 2019-12-05 17:11:42
% DurationCPUTime: 0.26s
% Computational Cost: add. (243->70), mult. (539->117), div. (0->0), fcn. (346->8), ass. (0->31)
t88 = sin(qJ(2));
t81 = qJD(2) * pkin(6) + t88 * qJD(1);
t97 = t81 * mrSges(4,3);
t96 = qJD(2) / 0.2e1;
t87 = sin(qJ(3));
t94 = pkin(7) * qJD(2) + t81;
t75 = qJD(3) * pkin(3) - t94 * t87;
t91 = cos(qJ(3));
t76 = t94 * t91;
t86 = sin(qJ(4));
t90 = cos(qJ(4));
t68 = t86 * t75 + t90 * t76;
t92 = cos(qJ(2));
t95 = t92 * qJD(1);
t84 = qJD(3) + qJD(4);
t67 = t90 * t75 - t86 * t76;
t79 = -t95 + (-pkin(3) * t91 - pkin(2)) * qJD(2);
t89 = cos(qJ(5));
t85 = sin(qJ(5));
t83 = qJD(5) + t84;
t82 = -qJD(2) * pkin(2) - t95;
t78 = (t86 * t91 + t87 * t90) * qJD(2);
t77 = (-t86 * t87 + t90 * t91) * qJD(2);
t71 = -t77 * pkin(4) + t79;
t70 = t85 * t77 + t89 * t78;
t69 = t89 * t77 - t85 * t78;
t66 = t77 * pkin(8) + t68;
t65 = t84 * pkin(4) - t78 * pkin(8) + t67;
t64 = t85 * t65 + t89 * t66;
t63 = t89 * t65 - t85 * t66;
t1 = m(4) * (t82 ^ 2 + (t87 ^ 2 + t91 ^ 2) * t81 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t71 ^ 2) / 0.2e1 + (m(3) * (t88 ^ 2 + t92 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * t84 / 0.2e1) * t84 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t83 / 0.2e1) * t83 + (Ifges(4,3) * qJD(3) / 0.2e1 + (-mrSges(4,1) * t87 - mrSges(4,2) * t91) * t81) * qJD(3) + (t79 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * t84 + Ifges(5,1) * t78 / 0.2e1) * t78 + (t71 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t83 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t79 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t78 + Ifges(5,6) * t84 + Ifges(5,2) * t77 / 0.2e1) * t77 + (-t71 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t83 + Ifges(6,2) * t69 / 0.2e1) * t69 + (Ifges(3,3) * t96 + (mrSges(3,1) * t92 - mrSges(3,2) * t88) * qJD(1) + (-t82 * mrSges(4,1) + Ifges(4,6) * qJD(3) + (Ifges(4,2) * t96 + t97) * t91) * t91 + (Ifges(4,4) * t91 * qJD(2) + t82 * mrSges(4,2) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t96 + t97) * t87) * t87) * qJD(2);
T = t1;
