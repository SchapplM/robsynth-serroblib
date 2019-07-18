% Return the minimum parameter vector for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% MPV [10x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4PRRR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_convert_par2_MPV_fixb: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t26 = m(4) + m(5);
t1 = [m(2) + m(3) + t26; pkin(1) ^ 2 * t26 + Ifges(3,3); pkin(1) * t26 + mrSges(3,1); mrSges(3,2); pkin(2) ^ 2 * m(5) + Ifges(4,3); m(5) * pkin(2) + mrSges(4,1); mrSges(4,2); Ifges(5,3); mrSges(5,1); mrSges(5,2);];
MPV  = t1;
