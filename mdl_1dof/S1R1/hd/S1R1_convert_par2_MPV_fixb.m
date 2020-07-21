% Return the minimum parameter vector for
% S1R1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1]';
% m [2x1]
%   mass of all robot links (including the base)
% mrSges [2x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [2x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [3x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S1R1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(2,1),zeros(2,3),zeros(2,6)}
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1R1_convert_par2_MPV_fixb: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [2 1]), ...
  'S1R1_convert_par2_MPV_fixb: m has to be [2x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [2,3]), ...
  'S1R1_convert_par2_MPV_fixb: mrSges has to be [2x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [2 6]), ...
  'S1R1_convert_par2_MPV_fixb: Ifges has to be [2x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t1 = [Ifges(2,3); mrSges(2,1); mrSges(2,2);];
MPV = t1;