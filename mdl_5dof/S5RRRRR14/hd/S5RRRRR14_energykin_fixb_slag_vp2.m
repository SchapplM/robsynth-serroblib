% Calculate kinetic energy for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR14_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:17
% EndTime: 2024-09-27 18:42:18
% DurationCPUTime: 0.00s
% Computational Cost: add. (564->75), mult. (918->133), div. (0->0), fcn. (630->10), ass. (0->38)
t107 = sin(qJ(4));
t111 = cos(qJ(4));
t103 = qJD(1) + qJD(2);
t105 = cos(pkin(5));
t101 = t105 * t103 + qJD(3);
t108 = sin(qJ(3));
t104 = sin(pkin(5));
t117 = t103 * t104;
t116 = pkin(9) * t117;
t112 = cos(qJ(3));
t113 = cos(qJ(2));
t119 = pkin(1) * qJD(1);
t99 = pkin(2) * t103 + t113 * t119;
t118 = t105 * t99;
t95 = t112 * t118;
t109 = sin(qJ(2));
t97 = pkin(8) * t117 + t109 * t119;
t85 = pkin(3) * t101 + t95 + (-t97 - t116) * t108;
t90 = t108 * t118 + t112 * t97;
t88 = t112 * t116 + t90;
t80 = t107 * t85 + t111 * t88;
t100 = qJD(4) + t101;
t79 = -t107 * t88 + t111 * t85;
t91 = (-pkin(3) * t103 * t112 - t99) * t104;
t110 = cos(qJ(5));
t106 = sin(qJ(5));
t98 = qJD(5) + t100;
t93 = (t107 * t112 + t108 * t111) * t117;
t92 = (-t107 * t108 + t111 * t112) * t117;
t89 = -t108 * t97 + t95;
t86 = -pkin(4) * t92 + t91;
t84 = t106 * t92 + t110 * t93;
t83 = -t106 * t93 + t110 * t92;
t78 = pkin(10) * t92 + t80;
t77 = pkin(4) * t100 - pkin(10) * t93 + t79;
t76 = t106 * t77 + t110 * t78;
t75 = -t106 * t78 + t110 * t77;
t1 = m(4) * (t104 ^ 2 * t99 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t79 ^ 2 + t80 ^ 2 + t91 ^ 2) / 0.2e1 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,3) * t98 / 0.2e1) * t98 + (t91 * mrSges(5,2) - t79 * mrSges(5,3) + Ifges(5,1) * t93 / 0.2e1) * t93 + (Ifges(2,3) / 0.2e1 + m(3) * (t109 ^ 2 + t113 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t89 * mrSges(4,1) - t90 * mrSges(4,2) + Ifges(4,3) * t101 / 0.2e1) * t101 + (-t91 * mrSges(5,1) + t80 * mrSges(5,3) + Ifges(5,4) * t93 + Ifges(5,2) * t92 / 0.2e1) * t92 + (t86 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,5) * t98 + Ifges(6,1) * t84 / 0.2e1) * t84 + (-t86 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t84 + Ifges(6,6) * t98 + Ifges(6,2) * t83 / 0.2e1) * t83 + (t79 * mrSges(5,1) - t80 * mrSges(5,2) + Ifges(5,5) * t93 + Ifges(5,6) * t92 + Ifges(5,3) * t100 / 0.2e1) * t100 + (Ifges(3,3) * t103 / 0.2e1 + (mrSges(3,1) * t113 - mrSges(3,2) * t109) * t119 + ((-t99 * (-mrSges(4,1) * t112 + mrSges(4,2) * t108) + (Ifges(4,2) * t112 ^ 2 / 0.2e1 + (Ifges(4,4) * t112 + Ifges(4,1) * t108 / 0.2e1) * t108) * t103) * t104 + (-t89 * t108 + t90 * t112) * mrSges(4,3) + t101 * (Ifges(4,5) * t108 + Ifges(4,6) * t112)) * t104) * t103;
T = t1;
