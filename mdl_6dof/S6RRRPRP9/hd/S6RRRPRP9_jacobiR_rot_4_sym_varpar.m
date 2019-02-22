% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:01
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP9_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_jacobiR_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:00:58
% EndTime: 2019-02-22 12:00:58
% DurationCPUTime: 0.10s
% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t76 = sin(qJ(2));
t77 = sin(qJ(1));
t87 = t77 * t76;
t78 = cos(qJ(3));
t86 = t77 * t78;
t75 = sin(qJ(3));
t79 = cos(qJ(2));
t85 = t79 * t75;
t84 = t79 * t78;
t80 = cos(qJ(1));
t83 = t80 * t76;
t82 = t80 * t78;
t81 = t80 * t79;
t74 = t77 * t75 + t78 * t81;
t73 = t75 * t81 - t86;
t72 = -t80 * t75 + t77 * t84;
t71 = -t77 * t85 - t82;
t1 = [-t72, -t76 * t82, -t73, 0, 0, 0; t74, -t76 * t86, t71, 0, 0, 0; 0, t84, -t76 * t75, 0, 0, 0; -t87, t81, 0, 0, 0, 0; t83, t77 * t79, 0, 0, 0, 0; 0, t76, 0, 0, 0, 0; t71, -t75 * t83, t74, 0, 0, 0; t73, -t75 * t87, t72, 0, 0, 0; 0, t85, t76 * t78, 0, 0, 0;];
JR_rot  = t1;
