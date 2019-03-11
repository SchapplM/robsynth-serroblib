% Calculate kinetic energy for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:39
% EndTime: 2019-03-09 02:02:39
% DurationCPUTime: 0.33s
% Computational Cost: add. (236->78), mult. (444->109), div. (0->0), fcn. (210->6), ass. (0->29)
t86 = sin(pkin(9));
t94 = -pkin(1) * t86 - qJ(3);
t80 = t94 * qJD(1);
t99 = t80 ^ 2;
t98 = m(3) / 0.2e1;
t97 = cos(qJ(5));
t87 = cos(pkin(9));
t95 = -pkin(1) * t87 - pkin(2);
t76 = qJD(3) + (-pkin(7) + t95) * qJD(1);
t89 = sin(qJ(4));
t90 = cos(qJ(4));
t73 = qJD(2) * t90 + t76 * t89;
t70 = qJD(4) * pkin(8) + t73;
t74 = (pkin(4) * t89 - pkin(8) * t90 - t94) * qJD(1);
t88 = sin(qJ(5));
t67 = t70 * t97 + t74 * t88;
t96 = qJD(1) * t90;
t72 = -qJD(2) * t89 + t76 * t90;
t69 = -qJD(4) * pkin(4) - t72;
t66 = -t70 * t88 + t74 * t97;
t91 = qJD(2) ^ 2;
t82 = qJD(1) * t89 + qJD(5);
t79 = qJD(1) * t95 + qJD(3);
t78 = qJD(4) * t88 + t96 * t97;
t77 = -qJD(4) * t97 + t88 * t96;
t65 = pkin(5) * t77 - qJ(6) * t78 + t69;
t64 = qJ(6) * t82 + t67;
t63 = -pkin(5) * t82 + qJD(6) - t66;
t1 = m(4) * (t79 ^ 2 + t91 + t99) / 0.2e1 + t91 * t98 + m(5) * (t72 ^ 2 + t73 ^ 2 + t99) / 0.2e1 + m(6) * (t66 ^ 2 + t67 ^ 2 + t69 ^ 2) / 0.2e1 + m(7) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + (t72 * mrSges(5,1) - t73 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t66 * mrSges(6,1) - t63 * mrSges(7,1) - t67 * mrSges(6,2) + t64 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t82) * t82 + (t69 * mrSges(6,2) + t63 * mrSges(7,2) - t66 * mrSges(6,3) - t65 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t78 + (Ifges(7,4) + Ifges(6,5)) * t82) * t78 + (t69 * mrSges(6,1) + t65 * mrSges(7,1) - t64 * mrSges(7,2) - t67 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t77 + (-Ifges(6,6) + Ifges(7,6)) * t82 + (-Ifges(6,4) + Ifges(7,5)) * t78) * t77 + (t79 * mrSges(4,2) - t80 * mrSges(4,3) + (-t80 * mrSges(5,2) - t72 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t96 / 0.2e1) * t90 + (-t80 * mrSges(5,1) - t73 * mrSges(5,3) - Ifges(5,6) * qJD(4)) * t89 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t87 * mrSges(3,1) - t86 * mrSges(3,2) + (t86 ^ 2 + t87 ^ 2) * t98 * pkin(1)) * pkin(1) + (-Ifges(5,4) * t90 + Ifges(5,2) * t89 / 0.2e1) * t89) * qJD(1)) * qJD(1);
T  = t1;
