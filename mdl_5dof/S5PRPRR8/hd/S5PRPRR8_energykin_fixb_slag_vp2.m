% Calculate kinetic energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:20
% EndTime: 2019-12-05 16:02:20
% DurationCPUTime: 0.20s
% Computational Cost: add. (136->59), mult. (285->96), div. (0->0), fcn. (150->8), ass. (0->30)
t88 = sin(qJ(2));
t84 = sin(pkin(5));
t98 = qJD(1) * t84;
t95 = t88 * t98;
t78 = qJD(2) * qJ(3) + t95;
t99 = t78 ^ 2;
t91 = cos(qJ(2));
t94 = -t91 * t98 + qJD(3);
t74 = (-pkin(2) - pkin(7)) * qJD(2) + t94;
t87 = sin(qJ(4));
t90 = cos(qJ(4));
t85 = cos(pkin(5));
t97 = qJD(1) * t85;
t71 = t87 * t74 + t90 * t97;
t96 = qJD(2) * t90;
t70 = t74 * t90 - t87 * t97;
t93 = qJD(1) ^ 2;
t89 = cos(qJ(5));
t86 = sin(qJ(5));
t82 = t85 ^ 2 * t93;
t81 = qJD(2) * t87 + qJD(5);
t77 = qJD(4) * t86 + t89 * t96;
t76 = qJD(4) * t89 - t86 * t96;
t75 = -qJD(2) * pkin(2) + t94;
t72 = t95 + (pkin(4) * t87 - pkin(8) * t90 + qJ(3)) * qJD(2);
t69 = qJD(4) * pkin(8) + t71;
t68 = -qJD(4) * pkin(4) - t70;
t67 = t69 * t89 + t72 * t86;
t66 = -t69 * t86 + t72 * t89;
t1 = m(6) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t70 ^ 2 + t71 ^ 2 + t99) / 0.2e1 + m(4) * (t75 ^ 2 + t82 + t99) / 0.2e1 + m(3) * (t82 + (t88 ^ 2 + t91 ^ 2) * t93 * t84 ^ 2) / 0.2e1 + m(2) * t93 / 0.2e1 + (t66 * mrSges(6,1) - t67 * mrSges(6,2) + Ifges(6,3) * t81 / 0.2e1) * t81 + (t70 * mrSges(5,1) - t71 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t68 * mrSges(6,2) - t66 * mrSges(6,3) + Ifges(6,5) * t81 + Ifges(6,1) * t77 / 0.2e1) * t77 + (-t68 * mrSges(6,1) + t67 * mrSges(6,3) + Ifges(6,4) * t77 + Ifges(6,6) * t81 + Ifges(6,2) * t76 / 0.2e1) * t76 + (t75 * mrSges(4,2) + t78 * mrSges(4,3) + (mrSges(3,1) * t91 - mrSges(3,2) * t88) * t98 + (t78 * mrSges(5,2) - t70 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t96 / 0.2e1) * t90 + (t78 * mrSges(5,1) - t71 * mrSges(5,3) - Ifges(5,6) * qJD(4)) * t87 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + (-Ifges(5,4) * t90 + Ifges(5,2) * t87 / 0.2e1) * t87) * qJD(2)) * qJD(2);
T = t1;
