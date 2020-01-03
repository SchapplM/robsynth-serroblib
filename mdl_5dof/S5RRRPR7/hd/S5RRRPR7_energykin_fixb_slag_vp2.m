% Calculate kinetic energy for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:02
% EndTime: 2019-12-31 21:16:03
% DurationCPUTime: 0.40s
% Computational Cost: add. (470->86), mult. (1002->135), div. (0->0), fcn. (690->8), ass. (0->36)
t107 = -pkin(7) - pkin(6);
t106 = pkin(6) * mrSges(3,3);
t105 = cos(qJ(3));
t100 = cos(qJ(2));
t103 = qJD(1) * t100;
t98 = sin(qJ(2));
t104 = qJD(1) * t98;
t97 = sin(qJ(3));
t85 = -t105 * t103 + t97 * t104;
t86 = (t100 * t97 + t105 * t98) * qJD(1);
t91 = (-pkin(2) * t100 - pkin(1)) * qJD(1);
t75 = pkin(3) * t85 - qJ(4) * t86 + t91;
t89 = qJD(2) * pkin(2) + t107 * t104;
t90 = t107 * t103;
t80 = -t105 * t90 + t97 * t89;
t93 = qJD(2) + qJD(3);
t78 = qJ(4) * t93 + t80;
t94 = sin(pkin(9));
t95 = cos(pkin(9));
t69 = t94 * t75 + t95 * t78;
t68 = t95 * t75 - t78 * t94;
t79 = t105 * t89 + t97 * t90;
t77 = -t93 * pkin(3) + qJD(4) - t79;
t99 = cos(qJ(5));
t96 = sin(qJ(5));
t84 = qJD(5) + t85;
t82 = t86 * t95 + t93 * t94;
t81 = -t86 * t94 + t93 * t95;
t72 = t81 * t96 + t82 * t99;
t71 = t81 * t99 - t82 * t96;
t70 = -t81 * pkin(4) + t77;
t67 = pkin(8) * t81 + t69;
t66 = pkin(4) * t85 - pkin(8) * t82 + t68;
t65 = t66 * t96 + t67 * t99;
t64 = t66 * t99 - t67 * t96;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t79 ^ 2 + t80 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t70 ^ 2) / 0.2e1 + (t79 * mrSges(4,1) - t80 * mrSges(4,2) + Ifges(4,3) * t93 / 0.2e1) * t93 + (t64 * mrSges(6,1) - t65 * mrSges(6,2) + Ifges(6,3) * t84 / 0.2e1) * t84 + (t77 * mrSges(5,2) - t68 * mrSges(5,3) + Ifges(5,1) * t82 / 0.2e1) * t82 + (t91 * mrSges(4,2) - t79 * mrSges(4,3) + Ifges(4,5) * t93 + Ifges(4,1) * t86 / 0.2e1) * t86 + (-t77 * mrSges(5,1) + t69 * mrSges(5,3) + Ifges(5,4) * t82 + Ifges(5,2) * t81 / 0.2e1) * t81 + (t70 * mrSges(6,2) - t64 * mrSges(6,3) + Ifges(6,5) * t84 + Ifges(6,1) * t72 / 0.2e1) * t72 + (-t70 * mrSges(6,1) + t65 * mrSges(6,3) + Ifges(6,4) * t72 + Ifges(6,6) * t84 + Ifges(6,2) * t71 / 0.2e1) * t71 + (t91 * mrSges(4,1) + t68 * mrSges(5,1) - t69 * mrSges(5,2) - t80 * mrSges(4,3) - Ifges(4,4) * t86 + Ifges(5,5) * t82 - Ifges(4,6) * t93 + Ifges(5,6) * t81 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t85) * t85 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t100 ^ 2 + t98 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (-pkin(1) * mrSges(3,2) + (t106 + Ifges(3,1) / 0.2e1) * t98) * t98 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t98 + (t106 + Ifges(3,2) / 0.2e1) * t100) * t100) * qJD(1) + ((-pkin(6) * mrSges(3,1) + Ifges(3,5)) * t98 + (-pkin(6) * mrSges(3,2) + Ifges(3,6)) * t100) * qJD(2)) * qJD(1);
T = t1;
