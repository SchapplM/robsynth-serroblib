% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:18
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:18:50
% EndTime: 2019-02-22 12:18:50
% DurationCPUTime: 0.05s
% Computational Cost: add. (49->6), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
t83 = cos(qJ(1));
t82 = sin(qJ(1));
t81 = qJ(2) + qJ(3) + qJ(4);
t80 = cos(t81);
t79 = sin(t81);
t78 = t83 * t80;
t77 = t83 * t79;
t76 = t82 * t80;
t75 = t82 * t79;
t1 = [t83, 0, 0, 0, 0, 0; t82, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t76, t77, t77, t77, 0, 0; -t78, t75, t75, t75, 0, 0; 0, -t80, -t80, -t80, 0, 0; -t75, t78, t78, t78, 0, 0; t77, t76, t76, t76, 0, 0; 0, t79, t79, t79, 0, 0;];
JR_rot  = t1;
