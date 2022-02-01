% Calculate kinetic energy for
% S5RPRPR4
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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:21
% EndTime: 2022-01-23 09:22:21
% DurationCPUTime: 0.35s
% Computational Cost: add. (260->74), mult. (606->120), div. (0->0), fcn. (374->8), ass. (0->32)
t99 = m(3) / 0.2e1;
t88 = sin(pkin(8));
t82 = (pkin(1) * t88 + pkin(6)) * qJD(1);
t94 = cos(qJ(3));
t85 = t94 * qJD(2);
t92 = sin(qJ(3));
t98 = qJ(4) * qJD(1);
t74 = qJD(3) * pkin(3) + t85 + (-t82 - t98) * t92;
t77 = t92 * qJD(2) + t94 * t82;
t75 = t94 * t98 + t77;
t87 = sin(pkin(9));
t89 = cos(pkin(9));
t67 = t87 * t74 + t89 * t75;
t90 = cos(pkin(8));
t97 = -pkin(1) * t90 - pkin(2);
t66 = t89 * t74 - t87 * t75;
t78 = qJD(4) + (-pkin(3) * t94 + t97) * qJD(1);
t93 = cos(qJ(5));
t91 = sin(qJ(5));
t86 = qJD(3) + qJD(5);
t83 = t97 * qJD(1);
t80 = (t87 * t94 + t89 * t92) * qJD(1);
t79 = (-t87 * t92 + t89 * t94) * qJD(1);
t76 = -t92 * t82 + t85;
t70 = -t79 * pkin(4) + t78;
t69 = t91 * t79 + t93 * t80;
t68 = t93 * t79 - t91 * t80;
t65 = t79 * pkin(7) + t67;
t64 = qJD(3) * pkin(4) - t80 * pkin(7) + t66;
t63 = t91 * t64 + t93 * t65;
t62 = t93 * t64 - t91 * t65;
t1 = qJD(2) ^ 2 * t99 + m(5) * (t66 ^ 2 + t67 ^ 2 + t78 ^ 2) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t70 ^ 2) / 0.2e1 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t86 / 0.2e1) * t86 + (t78 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,1) * t80 / 0.2e1) * t80 + (-t78 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,4) * t80 + Ifges(5,2) * t79 / 0.2e1) * t79 + (t70 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t86 + Ifges(6,1) * t69 / 0.2e1) * t69 + (-t70 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t69 + Ifges(6,6) * t86 + Ifges(6,2) * t68 / 0.2e1) * t68 + (t76 * mrSges(4,1) + t66 * mrSges(5,1) - t77 * mrSges(4,2) - t67 * mrSges(5,2) + Ifges(5,5) * t80 + Ifges(5,6) * t79 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3)) * qJD(3) + (t83 * (-mrSges(4,1) * t94 + mrSges(4,2) * t92) + (-t76 * t92 + t77 * t94) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t92 + Ifges(4,6) * t94) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t90 * mrSges(3,1) - t88 * mrSges(3,2) + (t88 ^ 2 + t90 ^ 2) * t99 * pkin(1)) * pkin(1) + Ifges(4,2) * t94 ^ 2 / 0.2e1 + (Ifges(4,4) * t94 + Ifges(4,1) * t92 / 0.2e1) * t92) * qJD(1)) * qJD(1);
T = t1;
