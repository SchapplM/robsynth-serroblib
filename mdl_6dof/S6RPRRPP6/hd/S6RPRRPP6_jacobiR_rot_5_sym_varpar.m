% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:41
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPP6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:41:34
% EndTime: 2019-02-22 10:41:34
% DurationCPUTime: 0.04s
% Computational Cost: add. (40->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t72 = sin(qJ(3));
t73 = sin(qJ(1));
t80 = t73 * t72;
t71 = qJ(4) + pkin(9);
t69 = sin(t71);
t74 = cos(qJ(3));
t79 = t74 * t69;
t70 = cos(t71);
t78 = t74 * t70;
t75 = cos(qJ(1));
t77 = t75 * t72;
t76 = t75 * t74;
t68 = -t73 * t69 + t70 * t77;
t67 = t69 * t77 + t73 * t70;
t66 = t75 * t69 + t70 * t80;
t65 = -t69 * t80 + t75 * t70;
t1 = [t68, 0, t73 * t78, t65, 0, 0; t66, 0, -t70 * t76, t67, 0, 0; 0, 0, -t72 * t70, -t79, 0, 0; -t67, 0, -t73 * t79, -t66, 0, 0; t65, 0, t69 * t76, t68, 0, 0; 0, 0, t72 * t69, -t78, 0, 0; -t76, 0, t80, 0, 0, 0; -t73 * t74, 0, -t77, 0, 0, 0; 0, 0, t74, 0, 0, 0;];
JR_rot  = t1;
