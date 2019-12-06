% Calculate kinetic energy for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:49
% EndTime: 2019-12-05 15:00:49
% DurationCPUTime: 0.15s
% Computational Cost: add. (116->50), mult. (279->86), div. (0->0), fcn. (182->8), ass. (0->27)
t82 = sin(pkin(8));
t84 = cos(pkin(8));
t92 = qJD(1) * cos(qJ(3));
t93 = qJD(1) * sin(qJ(3));
t75 = t82 * t92 + t84 * t93;
t71 = qJD(3) * qJ(4) + t75;
t81 = sin(pkin(9));
t83 = cos(pkin(9));
t67 = t81 * qJD(2) + t83 * t71;
t94 = pkin(6) * qJD(3);
t73 = -t82 * t93 + t84 * t92;
t91 = qJD(4) - t73;
t90 = qJD(1) ^ 2;
t89 = qJD(2) ^ 2;
t87 = cos(qJ(5));
t85 = sin(qJ(5));
t80 = t83 * qJD(2);
t74 = (t81 * t87 + t83 * t85) * qJD(3);
t72 = (-t81 * t85 + t83 * t87) * qJD(3);
t70 = -qJD(3) * pkin(3) + t91;
t68 = (-pkin(4) * t83 - pkin(3)) * qJD(3) + t91;
t66 = -t81 * t71 + t80;
t65 = t83 * t94 + t67;
t64 = t80 + (-t71 - t94) * t81;
t63 = t85 * t64 + t87 * t65;
t62 = t87 * t64 - t85 * t65;
t1 = m(3) * (t89 + (t82 ^ 2 + t84 ^ 2) * t90) / 0.2e1 + m(2) * t90 / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t73 ^ 2 + t75 ^ 2 + t89) / 0.2e1 + (t68 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,1) * t74 / 0.2e1) * t74 + (-t68 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t74 + Ifges(6,2) * t72 / 0.2e1) * t72 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,5) * t74 + Ifges(6,6) * t72 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t70 * (-mrSges(5,1) * t83 + mrSges(5,2) * t81) - t75 * mrSges(4,2) + t73 * mrSges(4,1) + (Ifges(5,2) * t83 ^ 2 / 0.2e1 + Ifges(4,3) / 0.2e1 + (Ifges(5,4) * t83 + Ifges(5,1) * t81 / 0.2e1) * t81) * qJD(3) + (-t66 * t81 + t67 * t83) * mrSges(5,3)) * qJD(3);
T = t1;
