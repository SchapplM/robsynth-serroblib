% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:21
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPP4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:21:29
% EndTime: 2019-02-22 11:21:29
% DurationCPUTime: 0.03s
% Computational Cost: add. (40->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t84 = sin(qJ(2));
t85 = sin(qJ(1));
t92 = t85 * t84;
t83 = qJ(4) + pkin(9);
t81 = sin(t83);
t86 = cos(qJ(2));
t91 = t86 * t81;
t82 = cos(t83);
t90 = t86 * t82;
t87 = cos(qJ(1));
t89 = t87 * t84;
t88 = t87 * t86;
t80 = -t81 * t92 + t87 * t82;
t79 = t87 * t81 + t82 * t92;
t78 = t81 * t89 + t85 * t82;
t77 = t85 * t81 - t82 * t89;
t1 = [t80, t81 * t88, 0, -t77, 0, 0; t78, t85 * t91, 0, t79, 0, 0; 0, t84 * t81, 0, -t90, 0, 0; -t85 * t86, -t89, 0, 0, 0, 0; t88, -t92, 0, 0, 0, 0; 0, t86, 0, 0, 0, 0; t79, -t82 * t88, 0, t78, 0, 0; t77, -t85 * t90, 0, -t80, 0, 0; 0, -t84 * t82, 0, -t91, 0, 0;];
JR_rot  = t1;
