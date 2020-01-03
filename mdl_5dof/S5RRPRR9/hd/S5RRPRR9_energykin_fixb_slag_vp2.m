% Calculate kinetic energy for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:19:58
% EndTime: 2019-12-31 20:19:58
% DurationCPUTime: 0.41s
% Computational Cost: add. (428->86), mult. (1002->136), div. (0->0), fcn. (690->8), ass. (0->36)
t106 = pkin(6) * mrSges(3,3);
t105 = pkin(6) + qJ(3);
t100 = cos(qJ(2));
t103 = qJD(1) * t100;
t97 = sin(qJ(2));
t104 = qJD(1) * t97;
t93 = sin(pkin(9));
t94 = cos(pkin(9));
t85 = t94 * t103 - t93 * t104;
t86 = (t100 * t93 + t94 * t97) * qJD(1);
t91 = qJD(3) + (-pkin(2) * t100 - pkin(1)) * qJD(1);
t74 = -pkin(3) * t85 - pkin(7) * t86 + t91;
t89 = qJD(2) * pkin(2) - t105 * t104;
t90 = t105 * t103;
t79 = t93 * t89 + t94 * t90;
t77 = qJD(2) * pkin(7) + t79;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t68 = t96 * t74 + t99 * t77;
t67 = t99 * t74 - t77 * t96;
t78 = t89 * t94 - t93 * t90;
t76 = -qJD(2) * pkin(3) - t78;
t84 = qJD(4) - t85;
t98 = cos(qJ(5));
t95 = sin(qJ(5));
t82 = qJD(5) + t84;
t81 = qJD(2) * t96 + t86 * t99;
t80 = qJD(2) * t99 - t86 * t96;
t71 = t80 * t95 + t81 * t98;
t70 = t80 * t98 - t81 * t95;
t69 = -pkin(4) * t80 + t76;
t66 = pkin(8) * t80 + t68;
t65 = pkin(4) * t84 - pkin(8) * t81 + t67;
t64 = t65 * t95 + t66 * t98;
t63 = t65 * t98 - t66 * t95;
t1 = m(4) * (t78 ^ 2 + t79 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t69 ^ 2) / 0.2e1 + (t91 * mrSges(4,2) - t78 * mrSges(4,3) + Ifges(4,1) * t86 / 0.2e1) * t86 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * t84 / 0.2e1) * t84 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t82 / 0.2e1) * t82 + (-t91 * mrSges(4,1) + t79 * mrSges(4,3) + Ifges(4,4) * t86 + Ifges(4,2) * t85 / 0.2e1) * t85 + (t76 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * t84 + Ifges(5,1) * t81 / 0.2e1) * t81 + (t69 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t82 + Ifges(6,1) * t71 / 0.2e1) * t71 + (-t76 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t81 + Ifges(5,6) * t84 + Ifges(5,2) * t80 / 0.2e1) * t80 + (-t69 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t71 + Ifges(6,6) * t82 + Ifges(6,2) * t70 / 0.2e1) * t70 + (t78 * mrSges(4,1) - t79 * mrSges(4,2) + Ifges(4,5) * t86 + Ifges(4,6) * t85 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t97 + Ifges(3,6) * t100 + (-mrSges(3,1) * t97 - mrSges(3,2) * t100) * pkin(6)) * qJD(1)) * qJD(2) + (m(3) * (pkin(1) ^ 2 + (t100 ^ 2 + t97 ^ 2) * pkin(6) ^ 2) / 0.2e1 + Ifges(2,3) / 0.2e1 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + t106) * t97) * t97 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t97 + (Ifges(3,2) / 0.2e1 + t106) * t100) * t100) * qJD(1) ^ 2;
T = t1;
