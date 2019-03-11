% Return the minimum parameter vector for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [29x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPPR5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t143 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t131 = sin(pkin(10));
t133 = cos(pkin(10));
t142 = t131 * t133;
t141 = Ifges(6,4) * t142;
t140 = -pkin(8) * m(7) - mrSges(7,3);
t135 = pkin(8) ^ 2;
t137 = (pkin(5) ^ 2);
t122 = Ifges(6,2) + (t135 + t137) * m(7) + t143;
t125 = m(7) * t135 + Ifges(6,1) + t143;
t127 = t131 ^ 2;
t128 = t133 ^ 2;
t139 = t128 * t122 + t127 * t125 + Ifges(4,2) + Ifges(5,3);
t138 = (2 * pkin(7) * mrSges(4,3)) + t139 + 0.2e1 * t141;
t136 = pkin(7) ^ 2;
t134 = cos(pkin(9));
t132 = sin(pkin(9));
t126 = t140 * pkin(5) + Ifges(6,5);
t124 = -0.2e1 * t141;
t1 = [Ifges(2,3) + t134 ^ 2 * (Ifges(3,2) + ((pkin(2) ^ 2 + t136) * m(4)) + t138) + (0.2e1 * t134 * Ifges(3,4) + ((m(4) * t136) + Ifges(3,1) + t138) * t132) * t132; mrSges(2,1); mrSges(2,2); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2); pkin(7) * m(4) + mrSges(3,3) + mrSges(4,3); m(3) + m(4); (m(7) * t137) + Ifges(4,1) + Ifges(5,2) + Ifges(6,3) + t124 - t139; -Ifges(6,6) * t133 - t126 * t131 + Ifges(4,4) + Ifges(5,6); -Ifges(6,6) * t131 + t126 * t133 - Ifges(5,4) + Ifges(4,5); Ifges(4,6) - Ifges(5,5) - (t128 - t127) * Ifges(6,4) + (t122 - t125) * t142; t122 * t127 + t125 * t128 + Ifges(5,1) + Ifges(4,3) + t124; mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t140; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
