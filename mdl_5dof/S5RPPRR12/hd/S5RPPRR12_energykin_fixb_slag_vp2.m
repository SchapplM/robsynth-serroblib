% Calculate kinetic energy for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:54
% EndTime: 2019-12-31 18:06:54
% DurationCPUTime: 0.29s
% Computational Cost: add. (220->63), mult. (459->104), div. (0->0), fcn. (262->6), ass. (0->29)
t82 = cos(pkin(8));
t92 = t82 ^ 2;
t81 = sin(pkin(8));
t84 = sin(qJ(4));
t86 = cos(qJ(4));
t70 = (t81 * t86 + t82 * t84) * qJD(1);
t91 = m(3) / 0.2e1;
t74 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t89 = -pkin(6) * qJD(1) + t74;
t67 = t89 * t81;
t68 = t89 * t82;
t62 = t86 * t67 + t84 * t68;
t90 = t81 ^ 2 + t92;
t76 = qJD(1) * qJ(2) + qJD(3);
t72 = t81 * qJD(1) * pkin(3) + t76;
t61 = -t84 * t67 + t86 * t68;
t85 = cos(qJ(5));
t83 = sin(qJ(5));
t77 = -qJD(1) * pkin(1) + qJD(2);
t71 = (-t81 * t84 + t82 * t86) * qJD(1);
t69 = qJD(5) + t70;
t64 = t83 * qJD(4) + t85 * t71;
t63 = t85 * qJD(4) - t83 * t71;
t60 = t70 * pkin(4) - t71 * pkin(7) + t72;
t59 = qJD(4) * pkin(7) + t62;
t58 = -qJD(4) * pkin(4) - t61;
t57 = t85 * t59 + t83 * t60;
t56 = -t83 * t59 + t85 * t60;
t1 = t77 ^ 2 * t91 + m(5) * (t61 ^ 2 + t62 ^ 2 + t72 ^ 2) / 0.2e1 + m(4) * (t90 * t74 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + (t72 * mrSges(5,2) - t61 * mrSges(5,3) + Ifges(5,1) * t71 / 0.2e1) * t71 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * t69 / 0.2e1) * t69 - (-t72 * mrSges(5,1) + t62 * mrSges(5,3) + Ifges(5,4) * t71 - Ifges(5,2) * t70 / 0.2e1) * t70 + (t58 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * t69 + Ifges(6,1) * t64 / 0.2e1) * t64 + (-t58 * mrSges(6,1) + t57 * mrSges(6,3) + Ifges(6,4) * t64 + Ifges(6,6) * t69 + Ifges(6,2) * t63 / 0.2e1) * t63 + (t61 * mrSges(5,1) - t62 * mrSges(5,2) + Ifges(5,5) * t71 - Ifges(5,6) * t70 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t76 * (mrSges(4,1) * t81 + mrSges(4,2) * t82) + t77 * mrSges(3,2) - t90 * t74 * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t91 + mrSges(3,3)) * qJ(2) + Ifges(4,1) * t92 / 0.2e1 + (-Ifges(4,4) * t82 + Ifges(4,2) * t81 / 0.2e1) * t81) * qJD(1)) * qJD(1);
T = t1;
