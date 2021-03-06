% Calculate kinetic energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:13:47
% EndTime: 2019-03-08 20:13:47
% DurationCPUTime: 0.27s
% Computational Cost: add. (236->79), mult. (468->114), div. (0->0), fcn. (266->8), ass. (0->33)
t102 = sin(qJ(2));
t98 = sin(pkin(6));
t112 = qJD(1) * t98;
t109 = t102 * t112;
t91 = qJD(2) * qJ(3) + t109;
t113 = t91 ^ 2;
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t99 = cos(pkin(6));
t111 = qJD(1) * t99;
t105 = cos(qJ(2));
t108 = -t105 * t112 + qJD(3);
t87 = (-pkin(2) - pkin(8)) * qJD(2) + t108;
t82 = t101 * t87 + t104 * t111;
t80 = qJD(4) * pkin(9) + t82;
t85 = t109 + (pkin(4) * t101 - pkin(9) * t104 + qJ(3)) * qJD(2);
t76 = t100 * t85 + t103 * t80;
t110 = qJD(2) * t104;
t75 = -t100 * t80 + t103 * t85;
t81 = -t101 * t111 + t104 * t87;
t79 = -qJD(4) * pkin(4) - t81;
t107 = qJD(1) ^ 2;
t96 = t99 ^ 2 * t107;
t95 = t101 * qJD(2) + qJD(5);
t90 = t100 * qJD(4) + t103 * t110;
t89 = t103 * qJD(4) - t100 * t110;
t88 = -qJD(2) * pkin(2) + t108;
t77 = -t89 * pkin(5) + qJD(6) + t79;
t74 = t89 * qJ(6) + t76;
t73 = t95 * pkin(5) - t90 * qJ(6) + t75;
t1 = m(3) * (t96 + (t102 ^ 2 + t105 ^ 2) * t98 ^ 2 * t107) / 0.2e1 + m(2) * t107 / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t113) / 0.2e1 + m(4) * (t88 ^ 2 + t113 + t96) / 0.2e1 + m(7) * (t73 ^ 2 + t74 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t79 ^ 2) / 0.2e1 + (t81 * mrSges(5,1) - t82 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t75 * mrSges(6,1) + t73 * mrSges(7,1) - t76 * mrSges(6,2) - t74 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t95) * t95 + (t79 * mrSges(6,2) + t77 * mrSges(7,2) - t75 * mrSges(6,3) - t73 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t90 + (Ifges(6,5) + Ifges(7,5)) * t95) * t90 + (-t79 * mrSges(6,1) - t77 * mrSges(7,1) + t76 * mrSges(6,3) + t74 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t89 + (Ifges(6,6) + Ifges(7,6)) * t95 + (Ifges(6,4) + Ifges(7,4)) * t90) * t89 + (t88 * mrSges(4,2) + t91 * mrSges(4,3) + (mrSges(3,1) * t105 - mrSges(3,2) * t102) * t112 + (t91 * mrSges(5,2) - t81 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t110 / 0.2e1) * t104 + (t91 * mrSges(5,1) - t82 * mrSges(5,3) - Ifges(5,6) * qJD(4)) * t101 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + (-Ifges(5,4) * t104 + Ifges(5,2) * t101 / 0.2e1) * t101) * qJD(2)) * qJD(2);
T  = t1;
