% Calculate kinetic energy for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:14
% EndTime: 2019-12-05 18:11:14
% DurationCPUTime: 0.46s
% Computational Cost: add. (445->80), mult. (1162->129), div. (0->0), fcn. (862->8), ass. (0->36)
t96 = cos(pkin(9));
t109 = t96 ^ 2;
t108 = qJD(1) * (pkin(6) + qJ(2));
t107 = m(3) / 0.2e1;
t101 = cos(qJ(4));
t102 = cos(qJ(3));
t95 = sin(pkin(9));
t87 = t95 * t108;
t88 = t96 * t108;
t99 = sin(qJ(3));
t79 = -t102 * t87 - t88 * t99;
t86 = (t102 * t95 + t96 * t99) * qJD(1);
t75 = qJD(3) * pkin(3) - pkin(7) * t86 + t79;
t80 = t102 * t88 - t99 * t87;
t85 = (t102 * t96 - t95 * t99) * qJD(1);
t76 = pkin(7) * t85 + t80;
t98 = sin(qJ(4));
t68 = t101 * t76 + t98 * t75;
t94 = qJD(3) + qJD(4);
t67 = t101 * t75 - t76 * t98;
t89 = qJD(2) + (-pkin(2) * t96 - pkin(1)) * qJD(1);
t81 = -pkin(3) * t85 + t89;
t100 = cos(qJ(5));
t97 = sin(qJ(5));
t91 = qJD(5) + t94;
t90 = -qJD(1) * pkin(1) + qJD(2);
t78 = t101 * t86 + t85 * t98;
t77 = t101 * t85 - t86 * t98;
t71 = -pkin(4) * t77 + t81;
t70 = t100 * t78 + t77 * t97;
t69 = t100 * t77 - t78 * t97;
t66 = pkin(8) * t77 + t68;
t65 = pkin(4) * t94 - pkin(8) * t78 + t67;
t64 = t100 * t66 + t65 * t97;
t63 = t100 * t65 - t66 * t97;
t1 = t90 ^ 2 * t107 + m(5) * (t67 ^ 2 + t68 ^ 2 + t81 ^ 2) / 0.2e1 + m(4) * (t79 ^ 2 + t80 ^ 2 + t89 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t71 ^ 2) / 0.2e1 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * t94 / 0.2e1) * t94 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t91 / 0.2e1) * t91 + (t89 * mrSges(4,2) - t79 * mrSges(4,3) + Ifges(4,1) * t86 / 0.2e1) * t86 + (-t89 * mrSges(4,1) + t80 * mrSges(4,3) + Ifges(4,4) * t86 + Ifges(4,2) * t85 / 0.2e1) * t85 + (t81 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * t94 + Ifges(5,1) * t78 / 0.2e1) * t78 + (t71 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t91 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t81 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t78 + Ifges(5,6) * t94 + Ifges(5,2) * t77 / 0.2e1) * t77 + (-t71 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t91 + Ifges(6,2) * t69 / 0.2e1) * t69 + (t79 * mrSges(4,1) - t80 * mrSges(4,2) + Ifges(4,5) * t86 + Ifges(4,6) * t85 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t90 * (-mrSges(3,1) * t96 + mrSges(3,2) * t95) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t107 + mrSges(3,3)) * (t95 ^ 2 + t109) * qJ(2) + Ifges(3,2) * t109 / 0.2e1 + (Ifges(3,4) * t96 + Ifges(3,1) * t95 / 0.2e1) * t95) * qJD(1)) * qJD(1);
T = t1;
