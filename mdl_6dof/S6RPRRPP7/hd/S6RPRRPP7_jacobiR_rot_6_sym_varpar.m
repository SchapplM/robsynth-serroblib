% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:42
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPP7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_jacobiR_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:42:10
% EndTime: 2019-02-22 10:42:10
% DurationCPUTime: 0.04s
% Computational Cost: add. (15->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t71 = sin(qJ(3));
t72 = sin(qJ(1));
t82 = t72 * t71;
t73 = cos(qJ(4));
t81 = t72 * t73;
t70 = sin(qJ(4));
t74 = cos(qJ(3));
t80 = t74 * t70;
t79 = t74 * t73;
t75 = cos(qJ(1));
t78 = t75 * t71;
t77 = t75 * t73;
t76 = t75 * t74;
t69 = -t72 * t70 + t71 * t77;
t68 = t70 * t78 + t81;
t67 = t75 * t70 + t71 * t81;
t66 = t70 * t82 - t77;
t1 = [t69, 0, t72 * t79, -t66, 0, 0; t67, 0, -t73 * t76, t68, 0, 0; 0, 0, -t71 * t73, -t80, 0, 0; t68, 0, t72 * t80, t67, 0, 0; t66, 0, -t70 * t76, -t69, 0, 0; 0, 0, -t71 * t70, t79, 0, 0; t76, 0, -t82, 0, 0, 0; t72 * t74, 0, t78, 0, 0, 0; 0, 0, -t74, 0, 0, 0;];
JR_rot  = t1;
