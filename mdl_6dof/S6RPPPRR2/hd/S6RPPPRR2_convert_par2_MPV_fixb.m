% Return the minimum parameter vector for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPPPRR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t99 = (m(6) + m(7));
t108 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t103 = (pkin(5) ^ 2);
t107 = (t103 * m(7) + Ifges(6,2));
t106 = 2 * pkin(7) * mrSges(6,3) + t107;
t105 = pkin(8) * m(7) + mrSges(7,3);
t96 = sin(pkin(9));
t98 = cos(pkin(9));
t104 = t98 * mrSges(3,1) - t96 * mrSges(3,2);
t102 = pkin(7) ^ 2;
t101 = pkin(8) ^ 2;
t97 = cos(pkin(10));
t95 = sin(pkin(10));
t1 = [Ifges(2,3) + Ifges(3,3) + Ifges(4,1) + t97 ^ 2 * (t102 * t99 + Ifges(5,1) + t106) + (-0.2e1 * t97 * Ifges(5,4) + (Ifges(5,2) + (pkin(4) ^ 2 + t102) * t99 + t106) * t95) * t95 + 0.2e1 * t104 * pkin(1); mrSges(2,1) + t104; t96 * mrSges(3,1) + t98 * mrSges(3,2) + mrSges(2,2); m(3); mrSges(4,2); mrSges(4,3); m(4); pkin(4) * t99 + mrSges(5,1); mrSges(5,2); pkin(7) * t99 + mrSges(5,3) + mrSges(6,3); m(5) + t99; m(7) * t101 + Ifges(6,1) - t107 + t108; t105 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t101 + t103) * m(7) + t108; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t105; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
