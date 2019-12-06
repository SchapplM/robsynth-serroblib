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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:41:34
% EndTime: 2019-12-05 17:41:35
% DurationCPUTime: 0.33s
% Computational Cost: add. (243->69), mult. (583->113), div. (0->0), fcn. (372->8), ass. (0->32)
t97 = m(3) / 0.2e1;
t86 = sin(pkin(8));
t81 = (pkin(1) * t86 + qJ(3)) * qJD(1);
t87 = cos(pkin(9));
t83 = t87 * qJD(2);
t85 = sin(pkin(9));
t96 = pkin(6) * qJD(1);
t72 = t83 + (-t81 - t96) * t85;
t75 = t85 * qJD(2) + t87 * t81;
t73 = t87 * t96 + t75;
t90 = sin(qJ(4));
t92 = cos(qJ(4));
t65 = t90 * t72 + t92 * t73;
t88 = cos(pkin(8));
t95 = -pkin(1) * t88 - pkin(2);
t64 = t92 * t72 - t90 * t73;
t76 = qJD(3) + (-pkin(3) * t87 + t95) * qJD(1);
t91 = cos(qJ(5));
t89 = sin(qJ(5));
t84 = qJD(4) + qJD(5);
t80 = t95 * qJD(1) + qJD(3);
t78 = (t85 * t92 + t87 * t90) * qJD(1);
t77 = (-t85 * t90 + t87 * t92) * qJD(1);
t74 = -t85 * t81 + t83;
t68 = -t77 * pkin(4) + t76;
t67 = t89 * t77 + t91 * t78;
t66 = t91 * t77 - t89 * t78;
t63 = t77 * pkin(7) + t65;
t62 = qJD(4) * pkin(4) - t78 * pkin(7) + t64;
t61 = t89 * t62 + t91 * t63;
t60 = t91 * t62 - t89 * t63;
t1 = m(6) * (t60 ^ 2 + t61 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t74 ^ 2 + t75 ^ 2 + t80 ^ 2) / 0.2e1 + qJD(2) ^ 2 * t97 + (t60 * mrSges(6,1) - t61 * mrSges(6,2) + Ifges(6,3) * t84 / 0.2e1) * t84 + (t76 * mrSges(5,2) - t64 * mrSges(5,3) + Ifges(5,1) * t78 / 0.2e1) * t78 + (-t76 * mrSges(5,1) + t65 * mrSges(5,3) + Ifges(5,4) * t78 + Ifges(5,2) * t77 / 0.2e1) * t77 + (t68 * mrSges(6,2) - t60 * mrSges(6,3) + Ifges(6,5) * t84 + Ifges(6,1) * t67 / 0.2e1) * t67 + (-t68 * mrSges(6,1) + t61 * mrSges(6,3) + Ifges(6,4) * t67 + Ifges(6,6) * t84 + Ifges(6,2) * t66 / 0.2e1) * t66 + (t64 * mrSges(5,1) - t65 * mrSges(5,2) + Ifges(5,5) * t78 + Ifges(5,6) * t77 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t80 * (-mrSges(4,1) * t87 + mrSges(4,2) * t85) + (-t74 * t85 + t75 * t87) * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t88 * mrSges(3,1) - t86 * mrSges(3,2) + (t86 ^ 2 + t88 ^ 2) * t97 * pkin(1)) * pkin(1) + Ifges(4,2) * t87 ^ 2 / 0.2e1 + (Ifges(4,4) * t87 + Ifges(4,1) * t85 / 0.2e1) * t85) * qJD(1)) * qJD(1);
T = t1;
