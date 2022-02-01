% Calculate kinetic energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:11
% EndTime: 2022-01-23 09:16:12
% DurationCPUTime: 0.46s
% Computational Cost: add. (335->80), mult. (898->132), div. (0->0), fcn. (604->8), ass. (0->37)
t110 = m(3) / 0.2e1;
t101 = cos(qJ(4));
t95 = sin(pkin(8));
t97 = cos(pkin(8));
t83 = qJD(2) + (-pkin(2) * t97 - qJ(3) * t95 - pkin(1)) * qJD(1);
t96 = cos(pkin(9));
t82 = t96 * t83;
t94 = sin(pkin(9));
t73 = t82 + (-pkin(6) * t95 * t96 + (-qJ(2) * t94 - pkin(3)) * t97) * qJD(1);
t108 = t95 * qJD(1);
t105 = t94 * t108;
t106 = qJ(2) * qJD(1);
t104 = t97 * t106;
t78 = t96 * t104 + t94 * t83;
t76 = -pkin(6) * t105 + t78;
t99 = sin(qJ(4));
t68 = t101 * t76 + t99 * t73;
t107 = t97 * qJD(1);
t88 = t95 * t106 + qJD(3);
t84 = pkin(3) * t105 + t88;
t67 = t101 * t73 - t99 * t76;
t89 = qJD(4) - t107;
t100 = cos(qJ(5));
t98 = sin(qJ(5));
t91 = -qJD(1) * pkin(1) + qJD(2);
t87 = qJD(5) + t89;
t80 = (t101 * t96 - t94 * t99) * t108;
t79 = (-t101 * t94 - t96 * t99) * t108;
t77 = -t104 * t94 + t82;
t75 = -t79 * pkin(4) + t84;
t70 = t100 * t80 + t98 * t79;
t69 = t100 * t79 - t98 * t80;
t66 = t79 * pkin(7) + t68;
t65 = t89 * pkin(4) - t80 * pkin(7) + t67;
t64 = t100 * t66 + t98 * t65;
t63 = t100 * t65 - t98 * t66;
t1 = m(6) * (t63 ^ 2 + t64 ^ 2 + t75 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t84 ^ 2) / 0.2e1 + m(4) * (t77 ^ 2 + t78 ^ 2 + t88 ^ 2) / 0.2e1 + t91 ^ 2 * t110 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * t89 / 0.2e1) * t89 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t87 / 0.2e1) * t87 + (t84 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * t89 + Ifges(5,1) * t80 / 0.2e1) * t80 + (t75 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t87 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t84 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t80 + Ifges(5,6) * t89 + Ifges(5,2) * t79 / 0.2e1) * t79 + (-t75 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t87 + Ifges(6,2) * t69 / 0.2e1) * t69 + ((t88 * (mrSges(4,1) * t94 + mrSges(4,2) * t96) + t91 * mrSges(3,2) + (Ifges(4,1) * t96 ^ 2 / 0.2e1 + Ifges(3,1) / 0.2e1 + (-Ifges(4,4) * t96 + Ifges(4,2) * t94 / 0.2e1) * t94) * t108 + (-t77 * t96 - t78 * t94) * mrSges(4,3)) * t95 + (-t91 * mrSges(3,1) + t78 * mrSges(4,2) - t77 * mrSges(4,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t107 + (-Ifges(4,5) * t96 + Ifges(4,6) * t94 + Ifges(3,4)) * t108) * t97 + (Ifges(2,3) / 0.2e1 + (qJ(2) * t110 + mrSges(3,3)) * (t95 ^ 2 + t97 ^ 2) * qJ(2)) * qJD(1)) * qJD(1);
T = t1;
