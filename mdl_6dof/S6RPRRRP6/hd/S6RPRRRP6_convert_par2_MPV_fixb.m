% Return the minimum parameter vector for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% MPV [30x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRRP6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t124 = (m(5) + m(6));
t136 = (-Ifges(6,2) - Ifges(7,2));
t129 = (pkin(3) ^ 2);
t135 = (t129 * t124 + Ifges(4,2));
t125 = (pkin(9) ^ 2);
t128 = (pkin(4) ^ 2);
t134 = (Ifges(5,2) + (t125 + t128) * m(6));
t119 = (m(4) + t124);
t133 = 2 * pkin(7) * mrSges(4,3) + t135;
t132 = -pkin(9) * m(6) - mrSges(6,3);
t116 = (mrSges(5,3) - t132);
t131 = pkin(8) * t124 + t116;
t130 = 2 * pkin(9) * mrSges(6,3) + 2 * pkin(8) * t116 + t134 - t136;
t127 = pkin(7) ^ 2;
t126 = pkin(8) ^ 2;
t123 = cos(pkin(10));
t122 = sin(pkin(10));
t1 = [Ifges(2,3) + t123 ^ 2 * (Ifges(3,2) + (pkin(2) ^ 2 + t127) * t119 + t133) + (0.2e1 * t123 * Ifges(3,4) + (t127 * t119 + Ifges(3,1) + t133) * t122) * t122; mrSges(2,1); mrSges(2,2); pkin(2) * t119 + mrSges(3,1); mrSges(3,2); pkin(7) * t119 + mrSges(3,3) + mrSges(4,3); m(3) + t119; t126 * t124 + Ifges(4,1) + t130 - t135; t131 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t126 + t129) * t124 + t130; pkin(3) * t124 + mrSges(4,1); mrSges(4,2) - t131; t125 * m(6) + Ifges(5,1) - t134; Ifges(5,4); t132 * pkin(4) + Ifges(5,5); Ifges(5,6); t128 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) + Ifges(7,1) + t136; Ifges(6,4) + Ifges(7,4); Ifges(6,5) + Ifges(7,5); Ifges(6,6) + Ifges(7,6); 2 * pkin(5) * mrSges(7,1) + Ifges(6,3) + Ifges(7,3); mrSges(6,1) + mrSges(7,1); mrSges(6,2) + mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
