% Calculate kinetic energy for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:20
% EndTime: 2019-12-31 20:49:21
% DurationCPUTime: 0.20s
% Computational Cost: add. (280->67), mult. (409->102), div. (0->0), fcn. (210->6), ass. (0->26)
t71 = qJD(1) + qJD(2);
t85 = t71 / 0.2e1;
t74 = sin(qJ(2));
t82 = pkin(1) * qJD(1);
t68 = t71 * pkin(7) + t74 * t82;
t84 = t68 * mrSges(4,3);
t75 = cos(qJ(3));
t83 = t71 * t75;
t73 = sin(qJ(3));
t79 = qJ(4) * t71 + t68;
t62 = qJD(3) * pkin(3) - t79 * t73;
t63 = t79 * t75;
t72 = sin(pkin(8));
t81 = cos(pkin(8));
t59 = t72 * t62 + t81 * t63;
t76 = cos(qJ(2));
t80 = t76 * t82;
t58 = t81 * t62 - t72 * t63;
t64 = -t80 + qJD(4) + (-pkin(3) * t75 - pkin(2)) * t71;
t69 = -t71 * pkin(2) - t80;
t66 = (t72 * t75 + t81 * t73) * t71;
t65 = t72 * t73 * t71 - t81 * t83;
t57 = qJD(3) * qJ(5) + t59;
t56 = -qJD(3) * pkin(4) + qJD(5) - t58;
t55 = t65 * pkin(4) - t66 * qJ(5) + t64;
t1 = m(4) * (t69 ^ 2 + (t73 ^ 2 + t75 ^ 2) * t68 ^ 2) / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t74 ^ 2 + t76 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t64 * mrSges(5,2) + t56 * mrSges(6,2) - t58 * mrSges(5,3) - t55 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t66) * t66 + (t64 * mrSges(5,1) + t55 * mrSges(6,1) - t57 * mrSges(6,2) - t59 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t65 + (-Ifges(5,4) + Ifges(6,5)) * t66) * t65 + (Ifges(3,3) * t85 + (mrSges(3,1) * t76 - mrSges(3,2) * t74) * t82 + (-t69 * mrSges(4,1) + (Ifges(4,2) * t85 + t84) * t75) * t75 + (Ifges(4,4) * t83 + t69 * mrSges(4,2) + (Ifges(4,1) * t85 + t84) * t73) * t73) * t71 + (t58 * mrSges(5,1) - t56 * mrSges(6,1) - t59 * mrSges(5,2) + t57 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + (Ifges(4,5) * t73 + Ifges(4,6) * t75) * t71 + (-t73 * mrSges(4,1) - t75 * mrSges(4,2)) * t68 + (Ifges(6,4) + Ifges(5,5)) * t66 + (-Ifges(5,6) + Ifges(6,6)) * t65) * qJD(3);
T = t1;
