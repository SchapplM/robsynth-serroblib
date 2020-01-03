% Return the minimum parameter vector for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% MPV [21x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RRRP5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_convert_par2_MPV_fixb: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP5_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t67 = (-Ifges(4,2) - Ifges(5,3));
t66 = -pkin(6) * m(4) - mrSges(4,3);
t65 = (mrSges(3,3) - t66);
t64 = (pkin(2) ^ 2);
t63 = (pkin(6) ^ 2);
t62 = (m(3) + m(4));
t61 = (t63 + t64);
t1 = [Ifges(2,3) + Ifges(3,2) + 2 * pkin(6) * mrSges(4,3) + t61 * m(4) + 2 * pkin(5) * t65 + (pkin(1) ^ 2 + pkin(5) ^ 2) * t62 - t67; pkin(1) * t62 + mrSges(2,1); -pkin(5) * t62 + mrSges(2,2) - t65; Ifges(3,1) - Ifges(3,2) + (-t61 + t63) * m(4); Ifges(3,4); t66 * pkin(2) + Ifges(3,5); Ifges(3,6); t64 * m(4) + Ifges(3,3); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2); Ifges(4,1) + Ifges(5,1) + t67; Ifges(4,4) - Ifges(5,5); Ifges(4,5) + Ifges(5,4); Ifges(4,6) - Ifges(5,6); Ifges(4,3) + Ifges(5,2); mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5);];
MPV = t1;
