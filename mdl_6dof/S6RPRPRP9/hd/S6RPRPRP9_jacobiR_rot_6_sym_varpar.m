% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:29
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:29:48
% EndTime: 2019-02-22 10:29:48
% DurationCPUTime: 0.03s
% Computational Cost: add. (40->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t82 = sin(qJ(3));
t83 = sin(qJ(1));
t90 = t83 * t82;
t81 = pkin(9) + qJ(5);
t79 = sin(t81);
t84 = cos(qJ(3));
t89 = t84 * t79;
t80 = cos(t81);
t88 = t84 * t80;
t85 = cos(qJ(1));
t87 = t85 * t82;
t86 = t85 * t84;
t78 = -t83 * t79 + t80 * t87;
t77 = t79 * t87 + t83 * t80;
t76 = t85 * t79 + t80 * t90;
t75 = t79 * t90 - t85 * t80;
t1 = [t78, 0, t83 * t88, 0, -t75, 0; t76, 0, -t80 * t86, 0, t77, 0; 0, 0, -t82 * t80, 0, -t89, 0; -t86, 0, t90, 0, 0, 0; -t83 * t84, 0, -t87, 0, 0, 0; 0, 0, t84, 0, 0, 0; t77, 0, t83 * t89, 0, t76, 0; t75, 0, -t79 * t86, 0, -t78, 0; 0, 0, -t82 * t79, 0, t88, 0;];
JR_rot  = t1;
