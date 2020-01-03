% Calculate kinetic energy for
% S5RPPRR3
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:24
% EndTime: 2020-01-03 11:27:24
% DurationCPUTime: 0.32s
% Computational Cost: add. (243->69), mult. (583->113), div. (0->0), fcn. (372->8), ass. (0->32)
t99 = m(3) / 0.2e1;
t88 = sin(pkin(8));
t83 = (pkin(1) * t88 + qJ(3)) * qJD(1);
t89 = cos(pkin(9));
t85 = t89 * qJD(2);
t87 = sin(pkin(9));
t98 = pkin(6) * qJD(1);
t74 = t85 + (-t83 - t98) * t87;
t77 = t87 * qJD(2) + t89 * t83;
t75 = t89 * t98 + t77;
t92 = sin(qJ(4));
t94 = cos(qJ(4));
t67 = t92 * t74 + t94 * t75;
t90 = cos(pkin(8));
t97 = -pkin(1) * t90 - pkin(2);
t66 = t94 * t74 - t92 * t75;
t78 = qJD(3) + (-pkin(3) * t89 + t97) * qJD(1);
t93 = cos(qJ(5));
t91 = sin(qJ(5));
t86 = qJD(4) + qJD(5);
t82 = t97 * qJD(1) + qJD(3);
t80 = (t87 * t94 + t89 * t92) * qJD(1);
t79 = (-t87 * t92 + t89 * t94) * qJD(1);
t76 = -t87 * t83 + t85;
t70 = -t79 * pkin(4) + t78;
t69 = t91 * t79 + t93 * t80;
t68 = t93 * t79 - t91 * t80;
t65 = t79 * pkin(7) + t67;
t64 = qJD(4) * pkin(4) - t80 * pkin(7) + t66;
t63 = t91 * t64 + t93 * t65;
t62 = t93 * t64 - t91 * t65;
t1 = m(4) * (t76 ^ 2 + t77 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t70 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t78 ^ 2) / 0.2e1 + qJD(2) ^ 2 * t99 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t86 / 0.2e1) * t86 + (t78 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,1) * t80 / 0.2e1) * t80 + (-t78 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,4) * t80 + Ifges(5,2) * t79 / 0.2e1) * t79 + (t70 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t86 + Ifges(6,1) * t69 / 0.2e1) * t69 + (-t70 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t69 + Ifges(6,6) * t86 + Ifges(6,2) * t68 / 0.2e1) * t68 + (t66 * mrSges(5,1) - t67 * mrSges(5,2) + Ifges(5,5) * t80 + Ifges(5,6) * t79 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t82 * (-mrSges(4,1) * t89 + mrSges(4,2) * t87) + (-t76 * t87 + t77 * t89) * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t90 * mrSges(3,1) - t88 * mrSges(3,2) + (t88 ^ 2 + t90 ^ 2) * t99 * pkin(1)) * pkin(1) + Ifges(4,2) * t89 ^ 2 / 0.2e1 + (Ifges(4,4) * t89 + Ifges(4,1) * t87 / 0.2e1) * t87) * qJD(1)) * qJD(1);
T = t1;
