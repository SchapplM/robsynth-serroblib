% Calculate kinetic energy for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:32
% EndTime: 2019-12-05 15:49:33
% DurationCPUTime: 0.18s
% Computational Cost: add. (161->60), mult. (370->101), div. (0->0), fcn. (230->10), ass. (0->30)
t101 = cos(qJ(4));
t94 = sin(pkin(5));
t107 = qJD(1) * t94;
t99 = sin(qJ(2));
t105 = t99 * t107;
t102 = cos(qJ(2));
t86 = qJD(2) * pkin(2) + t102 * t107;
t93 = sin(pkin(10));
t95 = cos(pkin(10));
t82 = t95 * t105 + t93 * t86;
t80 = qJD(2) * pkin(7) + t82;
t96 = cos(pkin(5));
t90 = t96 * qJD(1) + qJD(3);
t98 = sin(qJ(4));
t76 = t101 * t80 + t98 * t90;
t106 = qJD(2) * t98;
t81 = -t93 * t105 + t95 * t86;
t75 = t101 * t90 - t98 * t80;
t100 = cos(qJ(5));
t97 = sin(qJ(5));
t91 = -t101 * qJD(2) + qJD(5);
t85 = t97 * qJD(4) + t100 * t106;
t84 = t100 * qJD(4) - t97 * t106;
t79 = -qJD(2) * pkin(3) - t81;
t77 = (-pkin(4) * t101 - pkin(8) * t98 - pkin(3)) * qJD(2) - t81;
t74 = qJD(4) * pkin(8) + t76;
t73 = -qJD(4) * pkin(4) - t75;
t72 = t100 * t74 + t97 * t77;
t71 = t100 * t77 - t97 * t74;
t1 = m(6) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(4) * (t81 ^ 2 + t82 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t75 ^ 2 + t76 ^ 2 + t79 ^ 2) / 0.2e1 + (t71 * mrSges(6,1) - t72 * mrSges(6,2) + Ifges(6,3) * t91 / 0.2e1) * t91 + (t75 * mrSges(5,1) - t76 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t73 * mrSges(6,2) - t71 * mrSges(6,3) + Ifges(6,5) * t91 + Ifges(6,1) * t85 / 0.2e1) * t85 + (-t73 * mrSges(6,1) + t72 * mrSges(6,3) + Ifges(6,4) * t85 + Ifges(6,6) * t91 + Ifges(6,2) * t84 / 0.2e1) * t84 + (m(3) * (t96 ^ 2 + (t102 ^ 2 + t99 ^ 2) * t94 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t81 * mrSges(4,1) - t82 * mrSges(4,2) + (mrSges(3,1) * t102 - mrSges(3,2) * t99) * t107 + (t79 * mrSges(5,2) - t75 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t106 / 0.2e1) * t98 + (-t79 * mrSges(5,1) + t76 * mrSges(5,3) + Ifges(5,6) * qJD(4)) * t101 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (Ifges(5,4) * t98 + Ifges(5,2) * t101 / 0.2e1) * t101) * qJD(2)) * qJD(2);
T = t1;
