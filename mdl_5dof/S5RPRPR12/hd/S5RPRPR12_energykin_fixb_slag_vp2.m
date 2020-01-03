% Calculate kinetic energy for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:22
% EndTime: 2019-12-31 18:29:22
% DurationCPUTime: 0.44s
% Computational Cost: add. (391->80), mult. (970->127), div. (0->0), fcn. (686->8), ass. (0->37)
t98 = cos(pkin(8));
t110 = t98 ^ 2;
t109 = m(3) / 0.2e1;
t108 = cos(qJ(3));
t107 = pkin(6) + qJ(2);
t100 = sin(qJ(3));
t104 = t98 * qJD(1);
t96 = sin(pkin(8));
t105 = qJD(1) * t96;
t86 = t100 * t105 - t108 * t104;
t87 = (t100 * t98 + t108 * t96) * qJD(1);
t90 = qJD(2) + (-pkin(2) * t98 - pkin(1)) * qJD(1);
t74 = t86 * pkin(3) - t87 * qJ(4) + t90;
t88 = t107 * t105;
t89 = t107 * t104;
t79 = -t100 * t88 + t108 * t89;
t77 = qJD(3) * qJ(4) + t79;
t95 = sin(pkin(9));
t97 = cos(pkin(9));
t68 = t95 * t74 + t97 * t77;
t67 = t97 * t74 - t95 * t77;
t78 = -t100 * t89 - t108 * t88;
t76 = -qJD(3) * pkin(3) + qJD(4) - t78;
t101 = cos(qJ(5));
t99 = sin(qJ(5));
t92 = -qJD(1) * pkin(1) + qJD(2);
t82 = qJD(5) + t86;
t81 = t95 * qJD(3) + t97 * t87;
t80 = t97 * qJD(3) - t95 * t87;
t71 = t101 * t81 + t99 * t80;
t70 = t101 * t80 - t99 * t81;
t69 = -t80 * pkin(4) + t76;
t66 = t80 * pkin(7) + t68;
t65 = t86 * pkin(4) - t81 * pkin(7) + t67;
t64 = t101 * t66 + t99 * t65;
t63 = t101 * t65 - t99 * t66;
t1 = m(5) * (t67 ^ 2 + t68 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t69 ^ 2) / 0.2e1 + m(4) * (t78 ^ 2 + t79 ^ 2 + t90 ^ 2) / 0.2e1 + t92 ^ 2 * t109 + (t90 * mrSges(4,2) - t78 * mrSges(4,3) + Ifges(4,1) * t87 / 0.2e1) * t87 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t82 / 0.2e1) * t82 + (t76 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,1) * t81 / 0.2e1) * t81 + (-t76 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t81 + Ifges(5,2) * t80 / 0.2e1) * t80 + (t69 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t82 + Ifges(6,1) * t71 / 0.2e1) * t71 + (t78 * mrSges(4,1) - t79 * mrSges(4,2) + Ifges(4,5) * t87 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t69 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t71 + Ifges(6,6) * t82 + Ifges(6,2) * t70 / 0.2e1) * t70 + (t90 * mrSges(4,1) + t67 * mrSges(5,1) - t68 * mrSges(5,2) - t79 * mrSges(4,3) - Ifges(4,4) * t87 + Ifges(5,5) * t81 - Ifges(4,6) * qJD(3) + Ifges(5,6) * t80 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t86) * t86 + (t92 * (-mrSges(3,1) * t98 + mrSges(3,2) * t96) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t109 + mrSges(3,3)) * (t96 ^ 2 + t110) * qJ(2) + Ifges(3,2) * t110 / 0.2e1 + (Ifges(3,4) * t98 + Ifges(3,1) * t96 / 0.2e1) * t96) * qJD(1)) * qJD(1);
T = t1;
