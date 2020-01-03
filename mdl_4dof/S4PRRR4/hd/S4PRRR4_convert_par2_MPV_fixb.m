% Return the minimum parameter vector for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [18x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4PRRR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_convert_par2_MPV_fixb: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t48 = (m(4) + m(5));
t52 = -pkin(6) * m(5) - mrSges(5,3);
t51 = (mrSges(4,3) - t52);
t50 = (pkin(3) ^ 2);
t49 = (pkin(6) ^ 2);
t47 = (t49 + t50);
t1 = [m(2) + m(3) + t48; Ifges(3,3) + Ifges(4,2) + Ifges(5,2) + 2 * pkin(6) * mrSges(5,3) + t47 * m(5) + 2 * pkin(5) * t51 + (pkin(2) ^ 2 + pkin(5) ^ 2) * t48; pkin(2) * t48 + mrSges(3,1); -pkin(5) * t48 + mrSges(3,2) - t51; Ifges(4,1) - Ifges(4,2) + (-t47 + t49) * m(5); Ifges(4,4); t52 * pkin(3) + Ifges(4,5); Ifges(4,6); t50 * m(5) + Ifges(4,3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); Ifges(5,1) - Ifges(5,2); Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3); mrSges(5,1); mrSges(5,2);];
MPV = t1;
