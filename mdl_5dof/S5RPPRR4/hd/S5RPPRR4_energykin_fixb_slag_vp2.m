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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:43:59
% EndTime: 2019-12-05 17:43:59
% DurationCPUTime: 0.47s
% Computational Cost: add. (335->80), mult. (898->132), div. (0->0), fcn. (604->8), ass. (0->37)
t109 = m(3) / 0.2e1;
t100 = cos(qJ(4));
t94 = sin(pkin(8));
t96 = cos(pkin(8));
t82 = qJD(2) + (-pkin(2) * t96 - qJ(3) * t94 - pkin(1)) * qJD(1);
t95 = cos(pkin(9));
t81 = t95 * t82;
t93 = sin(pkin(9));
t72 = t81 + (-pkin(6) * t94 * t95 + (-qJ(2) * t93 - pkin(3)) * t96) * qJD(1);
t107 = t94 * qJD(1);
t104 = t93 * t107;
t105 = qJD(1) * qJ(2);
t103 = t96 * t105;
t77 = t95 * t103 + t93 * t82;
t75 = -pkin(6) * t104 + t77;
t98 = sin(qJ(4));
t67 = t100 * t75 + t98 * t72;
t106 = t96 * qJD(1);
t87 = t94 * t105 + qJD(3);
t83 = pkin(3) * t104 + t87;
t66 = t100 * t72 - t98 * t75;
t88 = qJD(4) - t106;
t99 = cos(qJ(5));
t97 = sin(qJ(5));
t90 = -qJD(1) * pkin(1) + qJD(2);
t86 = qJD(5) + t88;
t79 = (t100 * t95 - t93 * t98) * t107;
t78 = (-t100 * t93 - t95 * t98) * t107;
t76 = -t103 * t93 + t81;
t74 = -t78 * pkin(4) + t83;
t69 = t97 * t78 + t99 * t79;
t68 = t99 * t78 - t97 * t79;
t65 = t78 * pkin(7) + t67;
t64 = t88 * pkin(4) - t79 * pkin(7) + t66;
t63 = t97 * t64 + t99 * t65;
t62 = t99 * t64 - t97 * t65;
t1 = m(6) * (t62 ^ 2 + t63 ^ 2 + t74 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t83 ^ 2) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t87 ^ 2) / 0.2e1 + t90 ^ 2 * t109 + (t66 * mrSges(5,1) - t67 * mrSges(5,2) + Ifges(5,3) * t88 / 0.2e1) * t88 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t86 / 0.2e1) * t86 + (t83 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,5) * t88 + Ifges(5,1) * t79 / 0.2e1) * t79 + (t74 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t86 + Ifges(6,1) * t69 / 0.2e1) * t69 + (-t83 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,4) * t79 + Ifges(5,6) * t88 + Ifges(5,2) * t78 / 0.2e1) * t78 + (-t74 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t69 + Ifges(6,6) * t86 + Ifges(6,2) * t68 / 0.2e1) * t68 + ((t87 * (mrSges(4,1) * t93 + mrSges(4,2) * t95) + t90 * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,1) * t95 ^ 2 / 0.2e1 + (-Ifges(4,4) * t95 + Ifges(4,2) * t93 / 0.2e1) * t93) * t107 + (-t76 * t95 - t77 * t93) * mrSges(4,3)) * t94 + (-t90 * mrSges(3,1) + t77 * mrSges(4,2) - t76 * mrSges(4,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t106 + (-Ifges(4,5) * t95 + Ifges(4,6) * t93 + Ifges(3,4)) * t107) * t96 + (Ifges(2,3) / 0.2e1 + (qJ(2) * t109 + mrSges(3,3)) * (t94 ^ 2 + t96 ^ 2) * qJ(2)) * qJD(1)) * qJD(1);
T = t1;
