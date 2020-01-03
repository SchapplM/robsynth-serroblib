% Calculate kinetic energy for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:38
% EndTime: 2019-12-31 19:27:38
% DurationCPUTime: 0.10s
% Computational Cost: add. (161->43), mult. (206->71), div. (0->0), fcn. (68->6), ass. (0->20)
t69 = qJD(1) + qJD(2);
t81 = t69 / 0.2e1;
t76 = cos(qJ(2));
t80 = pkin(1) * qJD(1);
t79 = -t76 * t80 + qJD(3);
t64 = (-pkin(2) - pkin(3)) * t69 + t79;
t74 = sin(qJ(2));
t67 = t69 * qJ(3) + t74 * t80;
t71 = sin(pkin(8));
t72 = cos(pkin(8));
t62 = t71 * t64 + t72 * t67;
t61 = t72 * t64 - t71 * t67;
t75 = cos(qJ(5));
t73 = sin(qJ(5));
t66 = -t69 * pkin(2) + t79;
t60 = -t69 * pkin(7) + t62;
t59 = t69 * pkin(4) - t61;
t58 = t73 * qJD(4) + t75 * t60;
t57 = t75 * qJD(4) - t73 * t60;
t1 = m(5) * (qJD(4) ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t74 ^ 2 + t76 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t57 * mrSges(6,1) - t58 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t62 * mrSges(5,2) - t61 * mrSges(5,1) - t66 * mrSges(4,1) + t67 * mrSges(4,3) + Ifges(5,3) * t81 + (mrSges(3,1) * t76 - mrSges(3,2) * t74) * t80 + (Ifges(6,2) * t75 * t81 + t59 * mrSges(6,1) - t58 * mrSges(6,3) - Ifges(6,6) * qJD(5)) * t75 + (-t59 * mrSges(6,2) + t57 * mrSges(6,3) - Ifges(6,5) * qJD(5)) * t73 + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(6,4) * t75 + Ifges(6,1) * t73 / 0.2e1) * t73) * t69) * t69;
T = t1;
