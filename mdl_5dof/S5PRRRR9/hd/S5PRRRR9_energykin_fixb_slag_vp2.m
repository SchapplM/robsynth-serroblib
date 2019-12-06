% Calculate kinetic energy for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:18:58
% EndTime: 2019-12-05 17:18:58
% DurationCPUTime: 0.30s
% Computational Cost: add. (264->74), mult. (578->123), div. (0->0), fcn. (383->10), ass. (0->36)
t105 = sin(qJ(4));
t109 = cos(qJ(4));
t106 = sin(qJ(3));
t110 = cos(qJ(3));
t103 = cos(pkin(5));
t116 = qJD(1) * t103;
t107 = sin(qJ(2));
t102 = sin(pkin(5));
t117 = qJD(1) * t102;
t96 = qJD(2) * pkin(7) + t107 * t117;
t90 = t106 * t116 + t110 * t96;
t86 = qJD(3) * pkin(8) + t90;
t111 = cos(qJ(2));
t113 = t111 * t117;
t91 = -t113 + (-pkin(3) * t110 - pkin(8) * t106 - pkin(2)) * qJD(2);
t80 = t105 * t91 + t109 * t86;
t115 = qJD(2) * t106;
t114 = t110 * qJD(2);
t79 = -t105 * t86 + t109 * t91;
t100 = qJD(4) - t114;
t89 = -t106 * t96 + t110 * t116;
t85 = -qJD(3) * pkin(3) - t89;
t108 = cos(qJ(5));
t104 = sin(qJ(5));
t98 = qJD(5) + t100;
t97 = -qJD(2) * pkin(2) - t113;
t95 = t105 * qJD(3) + t109 * t115;
t94 = t109 * qJD(3) - t105 * t115;
t83 = t104 * t94 + t108 * t95;
t82 = -t104 * t95 + t108 * t94;
t81 = -t94 * pkin(4) + t85;
t78 = t94 * pkin(9) + t80;
t77 = t100 * pkin(4) - t95 * pkin(9) + t79;
t76 = t104 * t77 + t108 * t78;
t75 = -t104 * t78 + t108 * t77;
t1 = m(4) * (t89 ^ 2 + t90 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t79 ^ 2 + t80 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) / 0.2e1 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,3) * t98 / 0.2e1) * t98 + (t85 * mrSges(5,2) - t79 * mrSges(5,3) + Ifges(5,1) * t95 / 0.2e1) * t95 + (t89 * mrSges(4,1) - t90 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t85 * mrSges(5,1) + t80 * mrSges(5,3) + Ifges(5,4) * t95 + Ifges(5,2) * t94 / 0.2e1) * t94 + (t81 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,5) * t98 + Ifges(6,1) * t83 / 0.2e1) * t83 + (-t81 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t83 + Ifges(6,6) * t98 + Ifges(6,2) * t82 / 0.2e1) * t82 + (m(3) * (t103 ^ 2 + (t107 ^ 2 + t111 ^ 2) * t102 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t79 * mrSges(5,1) - t80 * mrSges(5,2) + Ifges(5,5) * t95 + Ifges(5,6) * t94 + Ifges(5,3) * t100 / 0.2e1) * t100 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t111 - mrSges(3,2) * t107) * t117 + (-t97 * mrSges(4,1) + t90 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t114 / 0.2e1) * t110 + (t97 * mrSges(4,2) - t89 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t110 + Ifges(4,1) * t106 / 0.2e1) * qJD(2)) * t106) * qJD(2);
T = t1;
