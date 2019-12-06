% Calculate kinetic energy for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:13
% EndTime: 2019-12-05 16:25:13
% DurationCPUTime: 0.31s
% Computational Cost: add. (242->75), mult. (576->124), div. (0->0), fcn. (387->10), ass. (0->35)
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t97 = sin(pkin(10));
t99 = cos(pkin(10));
t89 = (t102 * t97 - t105 * t99) * qJD(2);
t111 = qJ(4) * qJD(2);
t103 = sin(qJ(2));
t98 = sin(pkin(5));
t112 = qJD(1) * t98;
t92 = qJD(2) * pkin(7) + t103 * t112;
t100 = cos(pkin(5));
t110 = qJD(1) * t100;
t95 = t105 * t110;
t81 = qJD(3) * pkin(3) + t95 + (-t92 - t111) * t102;
t86 = t102 * t110 + t105 * t92;
t82 = t105 * t111 + t86;
t77 = t97 * t81 + t99 * t82;
t106 = cos(qJ(2));
t109 = t106 * t112;
t76 = t99 * t81 - t97 * t82;
t87 = -t109 + qJD(4) + (-pkin(3) * t105 - pkin(2)) * qJD(2);
t104 = cos(qJ(5));
t101 = sin(qJ(5));
t93 = -qJD(2) * pkin(2) - t109;
t90 = (t102 * t99 + t105 * t97) * qJD(2);
t88 = qJD(5) + t89;
t85 = -t102 * t92 + t95;
t84 = t101 * qJD(3) + t104 * t90;
t83 = t104 * qJD(3) - t101 * t90;
t78 = t89 * pkin(4) - t90 * pkin(8) + t87;
t75 = qJD(3) * pkin(8) + t77;
t74 = -qJD(3) * pkin(4) - t76;
t73 = t101 * t78 + t104 * t75;
t72 = -t101 * t75 + t104 * t78;
t1 = m(4) * (t85 ^ 2 + t86 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t76 ^ 2 + t77 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + (t87 * mrSges(5,2) - t76 * mrSges(5,3) + Ifges(5,1) * t90 / 0.2e1) * t90 + (t72 * mrSges(6,1) - t73 * mrSges(6,2) + Ifges(6,3) * t88 / 0.2e1) * t88 - (-t87 * mrSges(5,1) + t77 * mrSges(5,3) + Ifges(5,4) * t90 - Ifges(5,2) * t89 / 0.2e1) * t89 + (t74 * mrSges(6,2) - t72 * mrSges(6,3) + Ifges(6,5) * t88 + Ifges(6,1) * t84 / 0.2e1) * t84 + (-t74 * mrSges(6,1) + t73 * mrSges(6,3) + Ifges(6,4) * t84 + Ifges(6,6) * t88 + Ifges(6,2) * t83 / 0.2e1) * t83 + (m(3) * (t100 ^ 2 + (t103 ^ 2 + t106 ^ 2) * t98 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t93 * (-mrSges(4,1) * t105 + mrSges(4,2) * t102) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) * t105 ^ 2 / 0.2e1 + (Ifges(4,4) * t105 + Ifges(4,1) * t102 / 0.2e1) * t102) * qJD(2) + (mrSges(3,1) * t106 - mrSges(3,2) * t103) * t112 + (-t85 * t102 + t86 * t105) * mrSges(4,3)) * qJD(2) + (t85 * mrSges(4,1) + t76 * mrSges(5,1) - t86 * mrSges(4,2) - t77 * mrSges(5,2) + Ifges(5,5) * t90 - Ifges(5,6) * t89 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (Ifges(4,5) * t102 + Ifges(4,6) * t105) * qJD(2)) * qJD(3);
T = t1;
